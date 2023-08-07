#### Functions for processing forced photometry from ZTF  
from astropy.coordinates import SkyCoord 
from astropy.table import Table, vstack, hstack
from astropy.io import ascii 
import numpy as np
import matplotlib.pylab as plt
from numpy import random
import pandas as pd

import random
import math
import scipy 
from scipy.optimize import curve_fit
from scipy.special import gamma
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from scipy.stats import median_abs_deviation

import sfdmap 
import extinction
from PyAstronomy import pyasl
from astropy.coordinates import Distance


from astropy import cosmology

from iminuit import Minuit



import pprint



############################################################################################################################################################################################################################


############################################################################################################################################################################################################################


############################################################################################################################################################################################################################

    
def to_mag(f):
    '''
    this funtion converts flux to magnitude
    
    parameters
    ----------
    f.  column or array?
    
    returns
    -------
    '''
    f = abs(f)
    return -2.5*np.log10(f)

def error_on_conv_to_mag(eflux,flux):
    '''
    this function computes the error on the conversion from flux to mag
    
    parameters
    ----------
    eflux   error on flux
    flux.   flux values
    returns
    -------
    '''
    eflux = abs(eflux)
    flux  = abs(flux)
    return (2.5*eflux)/(np.log(10)*flux)




def read_table_forced_phot(path_name):
        '''
        This function reads the table and sets the column names for the ZTF forced photometry
        
        
        '''


        colname = ['index', 'field', 'ccdid', 'qid', 'filter', 'pid', 'infobitssci', 'sciinpseeing', 'scibckgnd', 'scisigpix', 'zpmaginpsci', 'zpmaginpsciunc', 'zpmaginpscirms', 
                    'clrcoeff', 'clrcoeffunc', 'ncalmatches', 'exptime', 'adpctdif1', 'adpctdif2', 'diffmaglim', 'zpdiff', 'programid', 'jd', 'rfid', 'forcediffimflux', 'forcediffimfluxunc',
                    'forcediffimsnr', 'forcediffimchisq', 'forcediffimfluxap', 'forcediffimfluxuncap', 'forcediffimsnrap', 'aperturecorr', 'dnearestrefsrc', 'nearestrefmag', 
                    'nearestrefmagunc', 'nearestrefchi', 'nearestrefsharp', 'refjdstart', 'refjdend', 'procstatus']

    
        table            = ascii.read(path_name, names = colname)

        table.remove_row(0)

        return table
        
def remove_procstatus_61(table):
    '''
    This procstatus causes to have "null" measurements which are preventing the creation of the photometric table
    
    '''

    table = table[table['procstatus'] != '61']
        
    return table 


def procstatus0_only(table):
    '''
    This procstatus causes to have "null" measurements which are preventing the creation of the photometric table
    
    '''

    table = table[table['procstatus'] == '0']
        
    return table 



def filter_bad_seeing(table, seeing = 4):
    '''
    This function filters the FP table on the science image seeing value. 
    By default we filter out any measurement with a seeing > 4"

    parameters
    ----------
    table [astropy.table]
    seeing [float] seeing value in as 
    '''

    table = table[table['sciinpseeing']<=4]
    return table

    
def remove_procstatus_57(table):
    '''
    This procstatus causes to have "null" measurements which are preventing the creation of the photometric table
    
    '''
    zaple = []
    zaple = [x for x in table if x['procstatus'][:2] != '57' and x['procstatus'][3:5] != '57'  ]
    if len(zaple) > 0 :
                  
        pable = vstack(zaple)
    
        return pable
    else :
        pable = table[table['procstatus']=='0']
        return pable


def remove_null_measurements(table):
    '''
    This function removes null measurement, independantly of procstatus... 
    (note: in the future, I'd just rely on this and forget about these stupid procstatus which are inconsistent....! )
    '''
    #print(type(table['forcediffimflux'][3]))
    table = table[table['forcediffimflux'] != 'null']
    #table = table[table['forcediffimfluxunc'] != 'null']
    #table = table[table['forcediffimchisq'] != 'null']
    #print(table)
    
    return table

def remove_infobits_nonzero(table):
    '''
    Removes all the entries with infobits different from 0 (which are supposed to be baddly processed?)
    '''
    table = table[table['infobitssci'] == '0.0']
    return table


def filter_StN_ratio(table, threshold):
    '''
    This function filters out on the S/N provided for the difference flux estimation

    parameters
    ----------
    threshold [float] S/N ratio
    '''


    table = table[table['forcediffimsnr'] >= threshold]
    return table
        
def convert_table_float(table):
        '''
        All the entries in the forced photometry table are strings, so we need to convert the relevant ones into floats
        
        '''
            
        table['procstatus'] = [ x[0:2] for x in table['procstatus'] ]
            
        
        list_keys = list(table.colnames)


        #NOTE: sometimes the the measurements are not null, but the nearest reference source within 5" returns nothing so it still gives a null. 
        # Make sure to remove the null flux measurements before converting to nans and flaots


        for lacle in list_keys:
            for _ in range(len(table[lacle])):
                if table[lacle][_] == 'null':
                    table[lacle][_] = 'nan'

        for lacle in list_keys:
            # table[lacle] = ['nan' for x in table[lacle] if x == 'null'] 
            table[lacle] = [float(x) for x in table[lacle]]
        table['index'] = [int(x) for x in table['index']]
        #print(table)
        return table



def remove_big_err_point(table, threshold = 5):
    '''
    Removes data points with large error bars

    should implement a sigma clipping method? 
    '''

    thres  = threshold * np.median(table['extcorrforcediffimfluxunc'])

    _table = table[table['extcorrforcediffimfluxunc'] <= thres ]   

    return _table


def create_new_id(table):
    '''
    This function creates a unique identification key for the field/quadrant...etc 
    '''

    new_id = []
    for _ in range(len(table['field'])):
        bli = float(str(math.floor(table['field'][_]))+str(math.floor(table['ccdid'][_]))+str(math.floor(table['qid'][_])))
        new_id.append(bli)
    table['obs_id'] = new_id
    return table


def clean_baseline_phot(table):
    '''
    This function is used to clean the baseline of outliers. It looks for the baseline median and std and kicks out points which are 5*std away 
    from the median
    '''
    _keep       = table[table['tfromFD'] >= -2.5 ] # want to conserve the information about the rise if it happens earlier 

    _tempbase   = table[table['tfromFD'] < -2.5  ]

    if len(_tempbase) != 0:
        _med, _mad  = np.median(_tempbase['forcediffimflux']) , median_abs_deviation(_tempbase['forcediffimflux'])

        _tempbase   = _tempbase[(_tempbase['forcediffimflux'] < _med + 5 * _mad ) & (_tempbase['forcediffimflux'] > _med - 5 * _mad )]
    #print(_tempbase)

        return vstack([_tempbase,_keep])

    elif len(_tempbase) == 0:
        return table



def correct_4_baseline_offsets(table):
    '''
    This function corrects the baseline for any offset for each 

    NOTE: you need to create the unique ID before using this function. see "create_new_id" function

    parameter
    ---------
    table [astrppytable]

    '''
    
    obids       = list(np.unique(table['obs_id']))

    # a very disgusting way to create a new table with the same keys
    
    lacorrectab = []
    # lacorrectab.remove_rows(slice(0,len(table)))
    

    for obs_id in obids: 
        
        _namebase           = 'temp_'+str(math.floor(obs_id)) 
        locals()[_namebase] = table[table['obs_id'] == obs_id] 

        _temp       = locals().get(_namebase)

        _tempbase   = _temp[_temp['tfromFD'] < -2.5  ]

        if len(_tempbase) != 0: 
            _med                            = np.median(_tempbase['forcediffimflux']) 
        
            _temp['forcediffimflux'] = _temp['forcediffimflux'] - _med
            
            lacorrectab = vstack([lacorrectab,_temp])


        else: 
            lacorrectab = vstack([lacorrectab,_temp])

    if len(lacorrectab) !=0:
        # print(lacorrectab)
        lacorrectab.sort('jd')
        return lacorrectab

    elif len(lacorrectab) == 0: 
        print('I could not compute a baseline or correct it')
        return table


   

def correct_4_baseline_offset_withoutsep(table):
    '''
    This function corrects the baseline for any offset for each 

    NOTE: you need to create the unique ID before using this function. see "create_new_id" function

    parameter
    ---------
    table [astrppytable]

    '''
    
    
    
    _tempbase   = table[table['tfromFD'] < -2.5  ]


        

    if len(_tempbase) != 0: 
        _med                            = np.median(_tempbase['forcediffimflux']) 
    
        table['forcediffimflux'] = table['forcediffimflux'] - _med

        return table
        
        
    else: 
        print('Not enough data to compute baseline shift')
        return table
            
    
    





############################################################################################################################################################################################################################


############################################################################################################################################################################################################################


############################################################################################################################################################################################################################