

import yaml
import os
from pprint import pprint

import numpy as np
from scipy import optimize
from scipy.stats import linregress


class KineticsSeries:
    """
    Organize a set of kinetics data, which are fundamentally sets of 
    1d timeseries with associated metadata:
    
        * E_0: the initial enzyme concentration
        * S_0: the initial substrate concentration
        * dt:  the time delay between each data point
        
    Data are "indexed" by E_0 and S_0 and accessing them in this
    way is likely to be useful.
    """
    
    def __init__(self, yaml_file, prefix='', corrections=None, sumfxn=np.median):
        
        # _datasets is the main "under the hood" data structure
        # it maps : (protein_conc, substrate_conc) -->
        #              [ ...,
        #                { timeseries, dt, cntrl_i, blank_i }, 
        #                ..., 
        #              ]

        self._datasets = {}
        self._sumfxn = sumfxn


        # >> load corrections
        if corrections:
            corr_file = open(corrections, 'r')
            self._corrections = yaml.safe_load(corr_file)
        else:
            self._corrections = {}


        # >> load kinetics data
        if type(yaml_file) is str:
            yaml_file = open(yaml_file, 'r')
        self._yaml = yaml.safe_load(yaml_file)
        
        for data_file in self._yaml.keys():
            print('Loading: %s...' % data_file)
            if self._yaml[data_file]['exclude'] == 'all':
                print(' ... excluded')
            else:
                self._load_data(os.path.join(prefix, data_file),
                                np.array(self._yaml[data_file]['p_conc_uM']),
                                np.array(self._yaml[data_file]['s_conc_uM']),
                                self._yaml[data_file]['gain'],
                                np.array(self._yaml[data_file]['blank']),
                                self._yaml[data_file]['dt_s'],
                                self._yaml[data_file]['exclude']
                               )
        
        return


    def _rfu_to_conc(self, s_conc_uM, blank_i, dt, timeseries):
        """
        """

        # subtract blank
        if self._corrections.get('zero_mode', None):
            zero_mode = True
            timeseries = timeseries - timeseries[0]
        else:
            timeseries = timeseries - blank_i

        # baseline correction
        if self._corrections.get('baseline_correction', None):
            m = self._corrections['baseline_correction'].get(s_conc_uM, 0.0)
            if m == 0.0:
                print(' ! Warning, baseline corr for [S]=%f not found!' % s_conc_uM)
            bl = m * dt * np.arange(len(timeseries))
            timeseries = timeseries + bl
        
        # convert to concentration
        if self._corrections.get('rfu_to_conc', None):
            
            mdl = self._corrections['rfu_to_conc']['model']

            if mdl == 'powerlaw':
                p_line = self._corrections['rfu_to_conc']['params']
                params = np.array([float(e) for e in p_line])
                if zero_mode:
                    params[4] = 0.0
                ts_uM = _powerlaw(timeseries, s_conc_uM, *params)
            else:
                raise ValueError('rfu_to_conc::model = `%s` not implemented' % mdl)

        else: # no RFU --> conc info
            ts_uM = timeseries

        return ts_uM
    

    def _load_data(self, file_path, p_conc_uM, s_conc_uM, gain, blank, dt_s, exclude):

        if gain != 1400.0:
            print(' !!! CRITICAL WARNING !!! ')
            print(' GAIN: %f != 1400 !!!!!!! ' % gain)

        if (len(p_conc_uM) != 12) or (len(s_conc_uM) != 8):
            raise ValueError('include all 12 cols and 8 rows in data entry!',
                             len(p_conc_uM), len(s_conc_uM))

        data_block = np.genfromtxt(file_path, delimiter=',').reshape(-1, 8, 12)
        
        for i,s in enumerate(s_conc_uM):
            for j,p in enumerate(p_conc_uM):
                
                k = (p, s)

                # -1 is no data
                if (s == -1) or (p == -1):
                    continue
               
                # >> load the data, and mark any flagged entries
                ts = data_block[:,i,j]
                if np.any(np.isnan(ts)):
                    print(p, s, ts)
                    raise ValueError('nan in data -- did you label them correctly?')
              
                if [i+1, j+1] in exclude:
                    print(' ... excluding E=%.2f / S=%.2f' % (p, s))
                    exclude_entry = True
                else:
                    exclude_entry = False

                # >> correct the data & convert from RFU to [P]
                blank_i = self._sumfxn(data_block[:,int(blank[0])-1,int(blank[1])-1])
                ts_in_uM = self._rfu_to_conc(s, blank_i, dt_s, ts)

                # >> build final dictionary to hold data
                entry = {
                    'timeseries': ts_in_uM,
                    'dt':         dt_s,
                    'exclude':    exclude_entry,
                }
                
                if k in self._datasets.keys():
                    self._datasets[k].append(entry)
                else:
                    self._datasets[k] = [entry]

        return
    
    
    @property
    def protein_concs(self):
        return np.unique([item[0] for item in self._datasets.keys()])
    
    
    @property
    def substrate_concs(self):
        return np.unique([item[1] for item in self._datasets.keys()])
    
    
    @property
    def all_conditions(self):
        """
        (protein, substrate)
        """
        return list(self._datasets.keys())


    def fit_v0(self, r2_threshold=0.0, **kwargs):

        for k in self.all_conditions:
            for e in self.get(*k):

                e.update(kwargs)

                v0, b, stderr_v0, r2 = fit_linear_v0(**e)
                e['v0']        = v0
                e['b']         = b
                e['stderr_v0'] = stderr_v0
                e['r2']        = r2

                if r2 < r2_threshold:
                    e['exclude'] = True

                if e['v0'] < 0.0:
                    e['exclude'] = True

        return

    
    def get(self, prot_c, subs_c):
        """
        Returns
        -------         
        timeseries : np.ndarray
            time series data
            
        dt : float
            the time spacing, in seconds
            
        cntrl_i : float
            the relevant control intensity
        """
        k = (prot_c, subs_c)
        if k in self._datasets.keys():
            ret = self._datasets[k]
        else:
            raise KeyError(k, 'not in dataset')
        return ret


    def get_v0s(self, prot_c, subs_c):
        
        v0s = []
        for e in self.get(prot_c, subs_c):
            if 'v0' in e.keys() and (not e['exclude']):
                v0s.append(e['v0'])

        return np.array(v0s)


    def get_set_v0s(self, prot_cs, subs_cs):
        """
        get all for one protein conc
        """

        ss  = []
        ps  = []
        v0s = []

        for i,s in enumerate(subs_cs):
            for j,p in enumerate(prot_cs):

                v = self.get_v0s(p, s)
                ss.extend( [s,] * len(v) )
                ps.extend( [p,] * len(v) )
                v0s.extend(v)

        return np.array(ss), np.array(ps), np.array(v0s)


def fit_linear_v0(timeseries, dt=1.0, regions=None,
                  **kwargs):
    """
    Given a 1-d numpy array that represents fluorescence-vs-time,
    fit a line to the initial (linear) part of the time series.
    
    Parameters
    ----------
    timeseries : np.ndarray
        time series data

    dt : float
        the time spacing, in seconds

    regions : int or list of ints
        the first N data points to fit (linear region), if None
        or list will scan a set an choose the best
        
    Returns
    -------
    v0 : float
        the initial velocity, in s-1
        
    b : float
        the intercept from the fit
        
    stderr_v0 : float
        the error on v0
        
    r2 : float 
        the R-squared of the linear fit
    """
    
    N = len(timeseries)
    x = np.arange(N) * dt

    if regions is None:
        regions = np.array([4096, 2048, 1024, 528, 256, 128, 64, 32, 16, 8])
    elif type(regions) is int:
        regions = [regions]
    elif type(regions) is list:
        pass
    else:
        raise TypeError()

    r2s = []
    for fit_region in regions:

        v0, b, r, p, stderr_v0 = linregress(x[:fit_region], timeseries[:fit_region])
        r2s.append( r ** 2 )

    best = np.argmax(r2s)
    fit_region = regions[best]
    v0, b, r, p, stderr_v0 = linregress(x[:fit_region], timeseries[:fit_region])
    r2 = r ** 2

    return v0, b, stderr_v0, r2
    

def fit_mm(v0s, substrate_concs, enzyme_conc):
    """
    Non-linear LSQ fit of the Michaelis-Mentin equation:
    
        V = (E_0 * k_cat * S) / (K_m + S)
    
    Parameters
    ----------
    v0s: np.ndarray
        1D array of the intial velocities
        
    substrate_concs: np.ndarray
        1D array of the initial substrate concentrations, same shape as v0s
        
    enzyme_conc: float
        the total enzyme concentration
    
    Returns
    -------    
    k_cat : float
    K_m : float
    """
    
    S = np.array(substrate_concs)
    
    def mm(tup):
        k_cat, K_m = tup
        V = (enzyme_conc * k_cat * S) / (K_m + S)
        return V - v0s
    
    res = optimize.least_squares(mm, x0=(1.0, 1.0))
    
    return res['x']


def mm(E_0, S, k_cat, K_m):
    return (E_0 * k_cat * S) / (K_m + S)


def fit_haldane(v0s, substrate_concs, enzyme_conc):
    """
    Non-linear LSQ fit of the Haldane equation:
    
        V = (E_0 * k_cat * S) / (K_m + S + S^2 / K_i)

	References
	----------
	.[1] Reed et al Bioessays (2010) 32, 422-429.

    Parameters
    ----------
    v0s: np.ndarray
        1D array of the intial velocities
        
    substrate_concs: np.ndarray
        1D array of the initial substrate concentrations, same shape as v0s
        
    enzyme_conc: float
        the total enzyme concentration
    
	k4 : bool
		if True, then the ES2 complex has a turnover rate (k4)

    Returns
    -------    
    k_cat : float
    K_m : float
    """

    S = np.array(substrate_concs)
    E = np.array(enzyme_conc)

    def h_simple(tup):
        k_cat, K_m, K_i = tup
        V = (E * k_cat * S) / (K_m + S + np.square(S) / K_i)
        return V - v0s

    res = optimize.least_squares(h_simple, x0=(0.01, 10.0, 1.0),
             bounds=([0.0,]*3, [np.inf,]*3) )

    return res['x']


def haldane(E_0, S, k_cat, K_m, K_i):
	V = (E_0 * k_cat * S) / (K_m + S + np.square(S) / K_i)
	return V


def _powerlaw(ts, S, a, b, c, d, e):
    """
    Power law model for the inner filter effect:

        f = (a*S^b + c) * P + d*S + e

    This function takes a timeseries `ts` in RFU ("f" in eqn above) and
    converts it into concentration of product ("P"), accounting for the
    linear filter effect of the substrate "S".
    """

    ts_uM = (ts - d*S - e) / (a*np.power(S,b) + c)

    return ts_uM


if __name__ == '__main__': 
	ks = KineticsSeries(yaml_string, prefix='./wt')

	print(ks.protein_concs)
	print(ks.substrate_concs)
	print(ks.all_conditions)
	print(ks.get(2.0, 160.0))

