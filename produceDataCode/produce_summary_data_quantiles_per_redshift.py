# from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import time
import sys
import copy
#Quick fudge to make import from ../Scripts work
import sys
sys.path.append('../../common_code')

# from time import sleep
# from IPython.display import clear_output, display
import pandas as pd 
import string
from PostProcessingScripts import * 






        

def calculateConfidenceIntervals(BPSmodelName='A', MSSFRnameslist=MSSFRnameslist, DCOtype='BHNS',  path_dir = '/Volumes/SimonsFoundation/DataDCO/'): 
    """ 
    Calculate confidence intervals  (distribution quantiles) summary 
    input:
    

    """



    nZbins = 6
    Metallicity_lower_bins, Metallicity_upper_bins = np.zeros(nZbins), np.zeros(nZbins)
    metallicity_list = np.concatenate((np.asarray(metallicities_list),np.asarray([0.03])))

    



    for n_rough_Z_bins in range(6): 
        lowerZ, upperZ = 0 + 9*n_rough_Z_bins, 8 + 9*n_rough_Z_bins
        Metallicity_lower_bins[n_rough_Z_bins], Metallicity_upper_bins[n_rough_Z_bins]  = metallicity_list[lowerZ], metallicity_list[upperZ]



    # prepare DataFrame 
    xvarHeaders = ['Mass1', 'Mass2', 'tc',\
                   'log10(tc)', 'TotMass', 'ChirpMass', 'q', 'metallicitySystems', 'log10metallicitySystems', 'tdelay',\
                   'log10(tdelay)']

    xvarUnits = ['Msun', 'Msun', 'Gyr',\
                   'log10(Gyr)', 'Msun', 'Msun', '#', '#', '#', 'Gyr', 'log10(Gyr)']

    # quantiles that I want to know
    y_quantiles  =          [0.005,   0.05,   0.16,   0.25,   0.5,   0.75,   0.84,   0.95,  0.995]
    indexnames   = ['unit', '0.005', '0.05', '0.16', '0.25', '0.5', '0.75', '0.84', '0.95', '0.995']
    # nr of rows and columns that will be used:
    ncol_var = len(xvarHeaders)   

    ncol_MSSFR = len(Metallicity_lower_bins)
    ncol_Rate_det = 1

    nrows = len(y_quantiles) + 1 # +1 for units (see below)
    # store variables, and Observed and intrinsic rates for all MSSFR variations:
    ncol = ncol_var * (ncol_MSSFR) # for each MSSFR I want to give quantiles for each xparam 

    # do the following for all channels:
    adjustedChannelList = ['classic', 'stable B no CEE', 'vii', 'immediate CE',  r'double-core CE', 'other']
    for nrC, Channel in enumerate(adjustedChannelList): 
        ind_wanted = dictFormationChannelIndex[Channel]


        df_placeholder = np.zeros((nrows, ncol)) # will be filled in loop: 

        headernames=[]
        units=[]
        for ind_s, ss in enumerate(xvarHeaders):
            for ind_Zbin,  mssfr in enumerate(Metallicity_lower_bins):
                sss = ss + '_' + str(Metallicity_lower_bins[ind_Zbin]) + '<Z<' + str(Metallicity_upper_bins[ind_Zbin])
                headernames.append(sss)
                units.append(xvarUnits[ind_s])

        # store dataFrame with zeros that we will fill on the go:
        dfw = pd.DataFrame(data=df_placeholder, columns=headernames, index=indexnames)   
        # add units at first row (index=0)
        dfw.iloc[0]=units        
            
                    
            
            
        print('now at m=', BPSmodelName)






        # constants
        Zsolar=0.0142


        #####

        DCOname = DCOname_dict[DCOtype]
        path_ = path_dir
        path_ = path_ + alphabetDirDict[BPSmodelName] +'/'
        path  = path_ + 'COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'

        print(path)
                



        # read in data 
        fdata = h5.File(path)
        # obtain BH and NS masses
        M1 = fdata['doubleCompactObjects']['M1'][...].squeeze()
        M2 = fdata['doubleCompactObjects']['M2'][...].squeeze()
        MBH, MNS = obtainM1BHandM2BHassymetric(M1, M2)
        del M1
        del M2
        seeds = fdata['doubleCompactObjects']['seed'][...].squeeze()
        channels = fdata['doubleCompactObjects']['formaton channel'][...].squeeze()

        tc = fdata['doubleCompactObjects']['tc'][...].squeeze() /1000.  # in Gyr
        metallicitySystems = fdata['doubleCompactObjects']['Metallicity1'][...].squeeze()
        tform = fdata['doubleCompactObjects']['tform'][...].squeeze() /1000.  # devide by 1000 to get units in Gyr 
        tdelay = tc + tform  # delay time 



        xvarlist = [MBH, MNS, tc, np.log10(tc), (MBH+MNS), chirpmass(MBH, MNS), (MBH/MNS), metallicitySystems, np.log10(metallicitySystems), tdelay, np.log10(tdelay)]
        # gives 90%, 99% and median

        # obtain weights (STROOPWAFEL)
        weights_  = fdata['doubleCompactObjects']['weight'][...].squeeze() 


        print('doing lower bins:', Metallicity_lower_bins)
        for ind_Zbin, lowerZbin in enumerate(Metallicity_lower_bins):
            upperZbin = Metallicity_upper_bins[ind_Zbin]    
            # print('running for Zbin: [',lowerZbin, ',', upperZbin, ']' )        
            mask_ = (channels==ind_wanted) & (metallicitySystems>=lowerZbin) & (metallicitySystems<=upperZbin)


            
            print('now at Zbin', str(Metallicity_lower_bins[ind_Zbin]) + '<Z<' + str(Metallicity_upper_bins[ind_Zbin]))
            # should do bootstrapping here at some point.


            # too few data points             
            if np.sum(mask_)<=10:
                print('%s data points for KDE, this is below the threshold 10, not drawing KDE'%np.sum(mask_))
            else:

                for ind_xvar, xvar in enumerate(xvarlist):



                    xvar_quantiles = weighted_quantile(values=np.asarray(xvar)[mask_], quantiles=y_quantiles, \
                             sample_weight=weights_[mask_])

                  # Nrepeats=1
                  # boot_xvar = np.asarray(xvar)
                  # boot_index = np.arange(len(boot_xvar))
                  # bootstrap_array = np.zeros((Nrepeats, len(y_quantiles)))
                  # for ind_r, repeat in enumerate(range(Nrepeats)):
                  #     # (bootstrap) re-sample a random set of samples with replacement from existing samples
                  #     # do this by drawing random sample indecis nrs, each nr corresponds to a sample  
                  #     boot_randindex = np.random.choice(boot_index, size=len(boot_index), replace=True, p=None)
                  #     boot_randweight   = np.asarray(weights)[boot_randindex] # get the corresponding weights 
                  #     boot_randsamples   = np.asarray(boot_xvar)[boot_randindex] # get the corresponding samples
                  #     xvar_quantiles = weighted_quantile(values=boot_randsamples[mask_], quantiles=y_quantiles, \
                  #            sample_weight=boot_randweight)
                  #     bootstrap_array[ind_r] = xvar_quantiles
                  # xvar_quantiles_bootstrapped = np.mean(bootstrap_array, axis=0)

                    dfw_key = xvarHeaders[ind_xvar] + '_' + str(Metallicity_lower_bins[ind_Zbin]) + '<Z<' + str(Metallicity_upper_bins[ind_Zbin])
                    dfw[dfw_key][1:] = xvar_quantiles

                
            
        dfCSVname = pathQuantiles + 'distribution_percentiles_model_' + DCOtype  + '_' + BPSmodelName + '_' + Channel + '.csv' 
        dfw.to_csv(dfCSVname) 
        
        
        print()
        print('finished')

        fdata.close()
        
    return







# pathQuantiles = '/Users/floorbroekgaarden/Projects/GitHub/Double-Compact-Object-Mergers/dataFiles/summary_data_Fig_4_5_6/'
# for DCOtype in ['BHNS', 'BBH', 'BNS']:
#     print('---------------------------')
#     print('---------------------------')
#     print('---------------------------')
#     print('at DCO model', DCOtype)
#     for ind_m, BPSmodelName in enumerate(['A', 'B', 'C',  'D', 'E', 'F', 'G', 'H',  'I' ,'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q']):
#         print()
#         print('------------------------------------------')
#         calculateConfidenceIntervals(BPSmodelName=BPSmodelName, MSSFRnameslist=MSSFRnameslist, DCOtype=DCOtype, whichWeight='det')


pathQuantiles = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/summary_quantiles/'
for DCOtype in [ 'BNS', 'BBH','BHNS']:
    print('---------------------------')
    print('---------------------------')
    print('---------------------------')
    print('at DCO model', DCOtype)
    for ind_m, BPSmodelName in enumerate(['A', 'B', 'C',  'D', 'E', 'F', 'G', 'H',  'I' ,'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T']):
        print()
        print('------------------------------------------')
        calculateConfidenceIntervals(BPSmodelName=BPSmodelName, MSSFRnameslist=MSSFRnameslist, DCOtype=DCOtype)
















