

import numpy as np
import pandas as pd


def load_raw(file_path, cycles=620):
    
    record = False
    results = np.zeros([cycles, 96])
    
    with open(file_path, 'r') as f:
        for line in f:

            if line.startswith('[Raw Data]'):
                record = True

            elif line.startswith('Well'):
                # first line after [Raw Data]
                # Well    Well Position   Cycle   x4-m4
                pass

            elif line.strip() == '':
                record = False
                #break

            elif record:
                d = line.split()
                
                well   = int(d[0]) - 1
                cycle  = int(d[2]) - 1
                fluoro = float(d[3].replace(',', ''))
                
                results[cycle, well] = fluoro
            
    return results


def load_melt(file_path):
    
    # Well    Well Position   Reading Temperature     Fluorescence    Derivative      Target Name
    
    record = False
    series = {'temperature' : []}
    
    with open(file_path, 'r') as f:
        for line in f:
            
            if line.startswith('[Melt Curve Raw Data]'):
                record = True
                
            elif line.startswith('Well'):
                pass
            
            elif line.strip() == '':
                record = False
            
            elif record:
                
                l = line.split()

                name = l[1]
                read = int(l[2]) - 1
                temp = float(l[3])
                fluo = float(l[4].replace(',', ''))
                derv = float(l[5].replace(',', ''))
                
                if name not in series.keys():
                    series[name] = []    
                series[name].append(derv)

                if len(series.keys()) <= 2: # the first time through
                    series['temperature'].append(temp)
                else:
                    #print(series['temperature'][read], temp)
                    assert np.abs(series['temperature'][read] - temp) < 0.2
                    
    return series

