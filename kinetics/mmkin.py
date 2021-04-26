

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
        * cntrl_i: the intensity of a well with no enzyme (S_0 intensity)
        * blank_i: the intensity of a well with nothing
        
    Data are "indexed" by E_0 and S_0 and accessing them in this
    way is likely to be useful.
    """
    
    def __init__(self, yaml_file, prefix='', sumfxn=np.median):
        
        # _datasets is the main "under the hood" data structure
        # it maps : (protein_conc, substrate_conc) -->
        #              [ ...,
        #                { timeseries, dt, cntrl_i, blank_i }, 
        #                ..., 
        #              ]
        self._datasets = {}
        self._sumfxn = sumfxn
   
        if type(yaml_file) is str:
            yaml_file = open(yaml_file, 'r')
        self._yaml = yaml.safe_load(yaml_file)

        #pprint(self._yaml)
        
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
    

    def _load_data(self, file_path, p_conc_uM, s_conc_uM, gain, blank, dt_s, exclude):

        if (len(p_conc_uM) != 12) or (len(s_conc_uM) != 8):
            raise ValueError('include all 12 cols and 8 rows in data entry!',
                             len(p_conc_uM), len(s_conc_uM))

        data_block = np.genfromtxt(file_path, delimiter=',').reshape(-1, 8, 12)
        
        for i,s in enumerate(s_conc_uM):
            for j,p in enumerate(p_conc_uM):
                
                k = (p, s)

                # -1 is no data, 0 means no protein (control)
                if (s == -1) or (p == -1) or (p == 0):
                    continue
               
                if 0 in p_conc_uM:
                    cntrl_idx = int(np.where(p_conc_uM == 0)[0])
                    cntrl_i   = self._sumfxn(data_block[:,i,cntrl_idx])
                else:
                    cntrl_i   = 0.0
                
                ts = data_block[:,i,j] / gain
                if np.any(np.isnan(ts)):
                    print(p, s, ts)
                    raise ValueError('nan in data -- did you label them correctly?')
              
                if [i+1, j+1] in exclude:
                    print(' ... excluding E=%.2f / S=%.2f' % (p, s))
                    exclude_entry = True
                else:
                    exclude_entry = False

                entry = {
                    'timeseries': ts,
                    'dt':         dt_s,
                    'cntrl_i':    cntrl_i,
                    'blank_i':    self._sumfxn(data_block[:,int(blank[0])-1,int(blank[1])-1]),
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


    def fit_v0(self):
        for k in self.all_conditions:
            for e in self.get(*k):
                v0, b, stderr_v0, r2 = fit_linear_v0(**e)
                e['v0'] = v0
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

        blank_i : float
            relevant blank intensity
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
            if 'v0' in e.keys():
                v0s.append(e['v0'])

        return np.array(v0s)


    def get_set_v0s(self, prot_c, subc_concentrations):
        """
        get all for one protein conc
        """

        ss  = []
        v0s = []

        for s in subc_concentrations:
            v = self.get_v0s(prot_c, s)
            ss.extend( [s,] * len(v) )
            v0s.extend(v)

        return np.array(ss), np.array(v0s)


def fit_linear_v0(timeseries, dt=1.0, blank_i=0.0, cntrl_i=0.0, fit_region=15, **kwargs):
    """
    Given a 1-d numpy array that represents fluorescence-vs-time,
    fit a line to the initial (linear) part of the time series.
    
    Parameters
    ----------
    timeseries : np.ndarray
        time series data

    dt : float
        the time spacing, in seconds

    cntrl_i : float
        the relevant control intensity

    blank_i : float
        relevant blank intensity
        
    fit_region : int
        the first N data points to fit (linear region)
        
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

    v0, b, r, p, stderr_v0 = linregress(x[:fit_region], timeseries[:fit_region])
    #if (r ** 2) < 0.9:
    #    print('WARNING: high R^2 =', r**2)
        
    snr = v0 / stderr_v0
    #if snr < 10.0:
    #    print('WARNING: high SNR =', snr)
    
    return v0, b, stderr_v0, r**2
    

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

    def h_simple(tup):
        k_cat, K_m, K_i = tup
        V = (enzyme_conc * k_cat * S) / (K_m + S + np.square(S) / K_i)
        return V - v0s

    res = optimize.least_squares(h_simple, x0=(0.01, 10.0, 1.0),
             bounds=([0.0,]*3, [np.inf,]*3) )

    return res['x']


def haldane(E_0, S, k_cat, K_m, K_i):
	V = (E_0 * k_cat * S) / (K_m + S + np.square(S) / K_i)
	return V


if __name__ == '__main__': 
	ks = KineticsSeries(yaml_string, prefix='./wt')

	print(ks.protein_concs)
	print(ks.substrate_concs)
	print(ks.all_conditions)
	print(ks.get(2.0, 160.0))

