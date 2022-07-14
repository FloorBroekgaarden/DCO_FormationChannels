# from __future__ import print_function
# from __future__ import division # undo if in Python 2 
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import time
import sys
import copy
#Quick fudge to make import from ../Scripts work
import sys
sys.path.append('../Scripts')
import string

import ClassCosmicIntegrator  as CI #Given settings and redshifts returns rates (2D arrays) Loads the data
from PostProcessingScripts import * 
# from ClassFormationChannels_5mainchannels import * 

import pandas as pd
from astropy import units as u
from astropy import constants as const

from formation_channels import * 





MSSFRnameslist = []
MSSFRnameslist.append('000') # add phenomenological 

for ind_GSMF, GSMF in enumerate(GSMFs):
    ind_y = ind_GSMF + 1
    for ind_MZ, MZ in enumerate(MZs):
        ind_z = ind_MZ +1
        for ind_SFR, SFR in enumerate(SFRs):
            ind_x = ind_SFR+1

            MSSFRnameslist.append('%s%s%s'%(ind_x, ind_y, ind_z))


GSMFs = [1,2,3]
SFRs = [1,2,3]
MZs=[1,2,3]


MSSFRnameslistCSV = []
MSSFRnameslistCSV.append('.0.0.0') # add phenomenological 


for ind_GSMF, GSMF in enumerate(GSMFs):
    ind_y = ind_GSMF + 1
    for ind_MZ, MZ in enumerate(MZs):
        ind_z = ind_MZ +1

        for ind_SFR, SFR in enumerate(SFRs):
            ind_x = ind_SFR+1            




            MSSFRnameslistCSV.append('.%s.%s.%s'%(ind_x, ind_y, ind_z))








def initialize_CSV_files_GW_formation_channel(DCOname='BHBH',  GWname='GW150629'):

    # namesEMlist=[]

    iii=0
    

    # CREATE PANDAS FILE 
    nModels=26
    BPSnameslist = list(string.ascii_uppercase)[0:nModels]

    NAMES = []
    stringgg =  'FormationChannels_rate_' + GWname 


    nModels=26
    BPSnameslist = list(string.ascii_uppercase)[0:nModels]

    NAMES = []
    stringgg =  'FormationChannels_rate_' + GWname   
    writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv'


    for ind_l, BPSmodelName in enumerate(BPSnameslist):
        # str_z0 = str(L + ' intrinsic (z=0) [Gpc^{-3} yr^{-1}]')
        # str_obs = str(L + ' observed (design LVK) [yr^{-1}]')

        namez0 = BPSmodelName + '_' +'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
        nameObs = BPSmodelName + '_' +'All observed (design LVK) [yr^{-1}]'


        namez0_oth = BPSmodelName + '_' +'Other intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
        nameObs_oth = BPSmodelName +'_' + 'Other observed (design LVK) [yr^{-1}]'

        namez0_I = BPSmodelName + '_' +'channel I intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
        nameObs_I = BPSmodelName +'_' + 'channel I observed (design LVK) [yr^{-1}]'
        namez0_II = BPSmodelName + '_' +'channel II intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
        nameObs_II = BPSmodelName + '_' +'channel II observed (design LVK) [yr^{-1}]'
        namez0_III = BPSmodelName +'_' + 'channel III intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
        nameObs_III = BPSmodelName +'_' + 'channel III observed (design LVK) [yr^{-1}]'
        namez0_IV = BPSmodelName + '_' +'channel IV intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
        nameObs_IV = BPSmodelName + '_' +'channel IV observed (design LVK) [yr^{-1}]'
        namez0_V = BPSmodelName + '_' +'channel V intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
        nameObs_V = BPSmodelName + '_' +'channel V observed (design LVK) [yr^{-1}]'




        NAMES.append(namez0)
        NAMES.append(nameObs)

        NAMES.append(namez0_oth)
        NAMES.append(nameObs_oth)

        NAMES.append(namez0_I)
        NAMES.append(nameObs_I)
        NAMES.append(namez0_II)
        NAMES.append(nameObs_II)
        NAMES.append(namez0_III)
        NAMES.append(nameObs_III)
        NAMES.append(namez0_IV)
        NAMES.append(nameObs_IV)
        NAMES.append(namez0_V)
        NAMES.append(nameObs_V)






    datas=[]

    for i in range(len(BPSnameslist)):
        for ii in range(7):
            datas.append(np.zeros_like(np.zeros(len(MSSFRnameslist))))
            datas.append(np.zeros_like(np.zeros(len(MSSFRnameslist))))


    df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
    df.columns =   df.columns.map(str)
    df.index.names = ['xyz']
    df.columns.names = ['m']

    # print(df) 

    df.to_csv(writePath)



    # df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv')





def GW_credible_intervals(GW_name, mode):

    # GW_list = ['GW151226','GW170729', 'GW190517_055101', 'GW190412','GW191109_010717'  ,'GW191103_012549', 'GW191126_115259']

    
    dfCSVname= '/Users/floorbroekgaarden/Projects/GitHub/MRR_Project/dataFiles/GWTC_posterior_samples/' 
    dfCSVname_ = dfCSVname + 'posteriorSamples_' + GW_name  + '.csv' 


    df = pd.read_csv(dfCSVname_, index_col=0, skiprows=[1])
    mass_ratio = df['M2']/df['M1']
    total_mass = df['M2'] + df['M1']
    spin1 = df['spin1'] 
    spin2 = df['spin2'] 
    
    chi_eff = df['chi_eff']
    chirp_mass = chirpmass(df['M2'], df['M1'])
    y_quantiles = [0.05, 0.5, 0.95]

    if mode=='normal':
    
        total_mass_CI = weighted_quantile(values=total_mass, quantiles=y_quantiles)
        mass1_CI = weighted_quantile(values=df['M1'], quantiles=y_quantiles)
        mass2_CI = weighted_quantile(values=df['M2'], quantiles=y_quantiles)
        chirp_mass_CI = weighted_quantile(values=chirp_mass, quantiles=y_quantiles)
        mass_ratio_CI = weighted_quantile(values=mass_ratio, quantiles=y_quantiles)
        spin1_CI = weighted_quantile(values=spin1, quantiles=y_quantiles)    
        spin2_CI = weighted_quantile(values=spin2, quantiles=y_quantiles)    
        chi_eff_quantiles = weighted_quantile(values=chi_eff, quantiles=y_quantiles)    

    elif mode=='spin1_is_zero':
        mask_spin = (abs(spin1)<0.05) & (spin2>0.05) # non MRR, spin 2 is the spinning one, we want spin1 to be zero 
        total_mass_CI = weighted_quantile(values=total_mass[mask_spin], quantiles=y_quantiles)
        mass1_CI = weighted_quantile(values=df['M1'][mask_spin], quantiles=y_quantiles)
        mass2_CI = weighted_quantile(values=df['M2'][mask_spin], quantiles=y_quantiles)
        chirp_mass_CI = weighted_quantile(values=chirp_mass[mask_spin], quantiles=y_quantiles)
        mass_ratio_CI = weighted_quantile(values=mass_ratio[mask_spin], quantiles=y_quantiles)
        spin1_CI = weighted_quantile(values=spin1[mask_spin], quantiles=y_quantiles)    
        spin2_CI = weighted_quantile(values=spin2[mask_spin], quantiles=y_quantiles)    
        chi_eff_quantiles = weighted_quantile(values=chi_eff[mask_spin], quantiles=y_quantiles)   

    elif  mode=='spin2_is_zero':
        mask_spin = (abs(spin2)<0.05) & (spin1>0.05)  # MRR spin1 is the spinning one, we want the other one to be zero
        total_mass_CI = weighted_quantile(values=total_mass[mask_spin], quantiles=y_quantiles)
        mass1_CI = weighted_quantile(values=df['M1'][mask_spin], quantiles=y_quantiles)
        mass2_CI = weighted_quantile(values=df['M2'][mask_spin], quantiles=y_quantiles)
        chirp_mass_CI = weighted_quantile(values=chirp_mass[mask_spin], quantiles=y_quantiles)
        mass_ratio_CI = weighted_quantile(values=mass_ratio[mask_spin], quantiles=y_quantiles)
        spin1_CI = weighted_quantile(values=spin1[mask_spin], quantiles=y_quantiles)    
        spin2_CI = weighted_quantile(values=spin2[mask_spin], quantiles=y_quantiles)    
        chi_eff_quantiles = weighted_quantile(values=chi_eff[mask_spin], quantiles=y_quantiles)   

    return total_mass_CI, mass1_CI, mass2_CI, chirp_mass_CI, mass_ratio_CI, spin1_CI, spin2_CI, chi_eff_quantiles
    



def writeToRatesFile_GW_FormationCHannels(BPSmodelName='Z', DCOtype='BHNS',  GWnameList=['GW150629']):
    """writes NS-BH rate to CSV file for all models"""


    if DCOtype=='BHNS':
        DCOname='BHNS'
    elif DCOtype=='BBH':
        DCOname='BHBH'
    elif DCOtype=='BNS':
        DCOname='NSNS'


    # path for files 
    path_dir = '/Volumes/Andromeda2/DATA/AllDCO_bugfix/'
    path_ = path_dir
    path_ = path_ + alphabetDirDict[BPSmodelName] +'/'
    path  = path_ + 'COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'
            
    # read in data 
    fdata = h5.File(path, 'r')
    fDCO      = fdata['doubleCompactObjects'] # hdf5 file with the DCO information
    fSN       = fdata['supernovae']  # hdf5 file with the SN information
    #
    M1 = fDCO['M1'][...].squeeze()   # Compact object mass [Msun] of the initially more massive star
    M2 = fDCO['M2'][...].squeeze()  # Compact object mass [Msun] of the initially less massive star

    print('using indices')
    seedsDCO = fDCO['seed'][...].squeeze()  # get the seeds in the DCO file 
    seedsSN = fSN['randomSeed'][...].squeeze()    # get the seeds in the SN file 
    indices = np.sort(np.unique(seedsSN[1::2], return_index=True)[1])
    maskSNdco = np.in1d(seedsSN,  seedsDCO) # mask in the SNe files the SNe that correspond to our DCO
    whichSN = fSN['whichStar'][...].squeeze()[maskSNdco]  # this is 1 if the initially primary star goes SN and 2 if the secondary goes supernova     
    whichSN2 = whichSN[1::2][indices]

    # either SN2 = primary (1) and M1 is > M2, or SN2 = secondary & M1 < M2 
    # this takes into account (first term) rejuvenation 
    MRR_mask = ((whichSN2==1) & (M1>M2) ) | ((whichSN2==2) & (M1<M2)) 

    M1LVK, M2LVK = obtainM1BHandM2BHassymetric(M1, M2)
    chirp_mass = chirpmass(M1LVK, M2LVK)
    mass_ratio_LVK =  M2LVK/M1LVK

    del M1
    del M2
    del whichSN2
    del whichSN
    del maskSNdco
    del indices
    del seedsSN
    del seedsDCO
    del fDCO
    del fSN


    seeds = fdata['doubleCompactObjects']['seed'][...].squeeze()
    channels = identify_formation_channels(seeds=seeds, file=fdata)

    # fdata.close()



    spin = COspin(data_path=path, state='he_depletion')  # set class 
    spin.setCOMPASData() # reads in the COMPAS DCO parameters 
    spinMZAMS1, spinMZAMS2  = spin.BaveraSpin()


    spinLVKM1, spinLVKM2 = np.zeros_like(spinMZAMS1), np.zeros_like(spinMZAMS1)
    spinLVKM1[MRR_mask] = spinMZAMS2[MRR_mask]  # MRR so M1 comes from M2ZAMS, we assign it spin from M2ZAMS
    spinLVKM1[~MRR_mask] = spinMZAMS1[~MRR_mask]  # no MRR so M1 comes from M1ZAMS, we assign it spin from M1ZAMS
    spinLVKM2[MRR_mask] = spinMZAMS1[MRR_mask]   # MRR so M2 comes from M1ZAMS, we assign it spin from M1ZAMS
    spinLVKM2[~MRR_mask] = spinMZAMS2[~MRR_mask]   # no MRR so M2 comes from M2ZAMS, we assign it spin from M2ZAMS     

    chi_eff = ((spinLVKM1*M1LVK) + (spinLVKM2*M2LVK)) / (M1LVK + M2LVK)


    # spin_threshold = 0.05 # definition of "spinning BH"





    # get intrinsic weights
    fparam_intrinsic = 'weights_intrinsic'
    # get detected weights
    fparam_detected = 'weights_detected'


    headerDict_intrinsic = {6:'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  0:'Other intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  1:'channel I intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 2:'channel II intrinsic (z=0) [Gpc^{-3} yr^{-1}]',3:'channel III intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 4:'channel IV intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 5:'channel V intrinsic (z=0) [Gpc^{-3} yr^{-1}]'}
    
    headerDict_observed  = {6:'All observed (design LVK) [yr^{-1}]',  0:'Other observed (design LVK) [yr^{-1}]', 1:'channel I observed (design LVK) [yr^{-1}]', 2:'channel II observed (design LVK) [yr^{-1}]', 3:'channel III observed (design LVK) [yr^{-1}]', 4:'channel IV observed (design LVK) [yr^{-1}]',5:'channel V observed (design LVK) [yr^{-1}]',}    





    # fdata = h5.File(path)



    DCOSeeds = fdata['doubleCompactObjects']['seed'][...].squeeze()

     

    for ind_GW, GWname in enumerate(GWnameList):  

        stringgg =  'FormationChannels_rate_' + GWname 
        
        writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv'

        total_mass_CI, mass1_CI, mass2_CI, chirp_mass_CI, mass_ratio_CI, spin1_CI, spin2_CI, chi_eff_CI = GW_credible_intervals(GWname, mode='normal')
        mask_GW =  (total_mass_CI[0]<=(M1LVK+M2LVK))  & ((M1LVK+M2LVK)<=total_mass_CI[2])  &  (chirp_mass_CI[0]<=chirp_mass) & (chirp_mass<=chirp_mass_CI[2])  & (mass_ratio_CI[0]<= mass_ratio_LVK) & (mass_ratio_LVK<=mass_ratio_CI[2]) & (chi_eff_CI[0]<=chi_eff) & (chi_eff<=chi_eff_CI[2]) & (mass1_CI[0]<= M1LVK) & (M1LVK<=mass1_CI[2]) & (mass2_CI[0]<= M2LVK) & (M2LVK<=mass2_CI[2])
        # mask_GW_match = (mask_GW==1) 



        df = pd.read_csv(writePath, index_col=0)
    

        for nrC, Channel in enumerate(list(np.unique(channels))):    
            print('reaching inside channels ')          
        #           #Get the seeds that relate to sorted indices
            
            
            intrinsicRates = np.zeros(len(MSSFRnameslist))
            detectedRates = np.zeros(len(MSSFRnameslist))       


            for ind_mssfr, mssfr in enumerate(MSSFRnameslist):


                weightheader = 'w_' + mssfr
                w_int = fdata[fparam_intrinsic][weightheader][...].squeeze()
                w_det = fdata[fparam_detected][weightheader][...].squeeze()



                if nrC==6: 
                    # TOTAL RATE
                    mask_match = (mask_GW==1)
                    intrinsicRates[ind_mssfr] = np.sum(w_int[mask_match])
                    detectedRates[ind_mssfr]  = np.sum(w_det[mask_match])  
                else:
                    mask_C = (channels==Channel) & (mask_GW==1)
                    # CHANNEL RATE 
                    intrinsicRates[ind_mssfr] = np.sum(w_int[mask_C])
                    detectedRates[ind_mssfr]  = np.sum(w_det[mask_C])          

                

                
            
            namez0  = BPSmodelName + '_' + headerDict_intrinsic[Channel]
            nameObs = BPSmodelName + '_' + headerDict_observed[Channel]
            print('reaching just before assigning')
            df[namez0] = intrinsicRates
            df[nameObs] = detectedRates


        print(GWname)


        df.to_csv(writePath)


    fdata.close() 










GW_MRR_list = ['GW151226','GW170729', 'GW190517_055101', 'GW190412','GW191109_010717'  ,'GW191103_012549', 'GW191126_115259']


GWTC1 = ['GW150914', 'GW151012', 'GW151226', 'GW170104', 'GW170608', 'GW170729',  'GW170809', 'GW170814', 'GW170818', 'GW170823'] # 'GW170817', 
#     GWTC2 = ['GW190408_181802','GW190412','GW190413_052954','GW190413_134308','GW190421_213856',\
#     'GW190424_180648','GW190503_185404','GW190512_180714',\
#     'GW190513_205428','GW190514_065416','GW190517_055101','GW190519_153544','GW190521_074359',\
#     'GW190521','GW190527_092055','GW190602_175927','GW190620_030421','GW190630_185205','GW190701_203306',\
#     'GW190706_222641','GW190707_093326','GW190708_232457','GW190720_000836',\
#     'GW190727_060333','GW190728_064510','GW190731_140936','GW190803_022701','GW190828_063405',\
#     'GW190828_065509','GW190910_112807','GW190915_235702','GW190924_021846','GW190929_012149',\
#      'GW190930_133541', 'GW190425', 'GW190814', 'GW190426_152155']
    
GWTC2 = [ 'GW190408_181802', 'GW190412', 'GW190413_052954', 'GW190413_134308', 'GW190421_213856', 'GW190424_180648', 'GW190503_185404', 'GW190512_180714', 'GW190513_205428', 'GW190514_065416', 'GW190517_055101', 'GW190519_153544', 'GW190521', 'GW190521_074359', 'GW190527_092055', 'GW190602_175927', 'GW190620_030421', 'GW190630_185205', 'GW190701_203306', 'GW190706_222641', 'GW190707_093326', 'GW190708_232457', 'GW190719_215514', 'GW190720_000836', 'GW190727_060333', 'GW190728_064510', 'GW190731_140936', 'GW190803_022701', 'GW190814', 'GW190828_063405', 'GW190828_065509', 'GW190910_112807', 'GW190915_235702', 'GW190924_021846', 'GW190929_012149', 'GW190930_133541'] #, 'GW191103_012549', 'GW191105_143521', 'GW191109_010717', 'GW191113_071753', 'GW191126_115259', 'GW191127_050227', 'GW191129_134029', 'GW191204_110529', 'GW191204_171526', 'GW191215_223052', 'GW191216_213338', 'GW191222_033537', 'GW191230_180458', 'GW200112_155838', 'GW200128_022011', 'GW200129_065458', 'GW200202_154313', 'GW200208_130117', 'GW200208_222617', 'GW200209_085452', 'GW200210_092255', 'GW200216_220804', 'GW200219_094415', 'GW200220_061928', 'GW200220_124850', 'GW200224_222234', 'GW200225_060421', 'GW200302_015811', 'GW200306_093714', 'GW200308_173609', 'GW200311_115853', 'GW200316_215756', 'GW200322_091133']
#     GWTC3 = ['GW191103_012549', 'GW191126_115259','GW191109_010717' ]    
GWTC3 = ['GW191103_012549', 'GW191105_143521', 'GW191109_010717', 'GW191113_071753', 'GW191126_115259', 'GW191127_050227', 'GW191129_134029', 'GW191204_110529', 'GW191204_171526', 'GW191215_223052', 'GW191216_213338', 'GW191222_033537', 'GW191230_180458', 'GW200112_155838', 'GW200128_022011', 'GW200129_065458', 'GW200202_154313', 'GW200208_130117', 'GW200208_222617', 'GW200209_085452', 'GW200210_092254', 'GW200216_220804', 'GW200219_094415', 'GW200220_061928', 'GW200220_124850', 'GW200224_222234', 'GW200225_060421', 'GW200302_015811', 'GW200306_093714', 'GW200308_173609', 'GW200311_115853', 'GW200316_215756', 'GW200322_091133']             
    


BBHlist_GWTC_123 = ['GW150914', 'GW151012', 'GW151226', 'GW170104', 'GW170608', 'GW170729', 'GW170809', \
                    'GW170814', 'GW170818', 'GW170823', 'GW190408_181802', 'GW190412', 'GW190413_052954',\
                    'GW190413_134308', 'GW190421_213856', 'GW190424_180648', 'GW190503_185404', 'GW190512_180714', \
                    'GW190513_205428', 'GW190514_065416', 'GW190517_055101', 'GW190519_153544', 'GW190521', \
                    'GW190521_074359', 'GW190527_092055', 'GW190602_175927', 'GW190620_030421', 'GW190630_185205', \
                    'GW190701_203306', 'GW190706_222641', 'GW190707_093326', 'GW190708_232457', 'GW190719_215514',\
                    'GW190720_000836', 'GW190727_060333', 'GW190728_064510', 'GW190731_140936', 'GW190803_022701',\
                    'GW190814', 'GW190828_063405', 'GW190828_065509', 'GW190910_112807', 'GW190915_235702', \
                    'GW190924_021846', 'GW190929_012149', 'GW190930_133541', 'GW191103_012549', 'GW191105_143521',\
                    'GW191109_010717', 'GW191113_071753', 'GW191126_115259', 'GW191127_050227', 'GW191129_134029', \
                    'GW191204_110529', 'GW191204_171526', 'GW191215_223052', 'GW191216_213338', 'GW191222_033537', \
                    'GW191230_180458', 'GW200112_155838', 'GW200128_022011', 'GW200129_065458', 'GW200202_154313', \
                    'GW200208_130117', 'GW200208_222617', 'GW200209_085452', 'GW200210_092254', 'GW200216_220804', \
                    'GW200219_094415', 'GW200220_061928', 'GW200220_124850', 'GW200224_222234', 'GW200225_060421',\
                    'GW200302_015811', 'GW200306_093714', 'GW200308_173609', 'GW200311_115853', 'GW200316_215756', 'GW200322_091133']



# ['GW150914', 'GW151012', 'GW151226', 'GW170104', 'GW170608', 'GW170729', 'GW170809', 'GW170814', 'GW170818', 'GW170823', 'GW190408_181802', 'GW190412', 'GW190413_052954', 'GW190413_134308', 'GW190421_213856', 'GW190424_180648', 'GW190503_185404', 'GW190512_180714', 'GW190513_205428', 'GW190514_065416', 'GW190517_055101', 'GW190519_153544', 'GW190521', 'GW190521_074359', 'GW190527_092055', 'GW190602_175927', 'GW190620_030421', 'GW190630_185205', 'GW190701_203306', 'GW190706_222641', 'GW190707_093326', 'GW190708_232457', 'GW190719_215514', 'GW190720_000836', 'GW190727_060333', 'GW190728_064510', 'GW190731_140936', 'GW190803_022701', 'GW190814', 'GW190828_063405', 'GW190828_065509', 'GW190910_112807', 'GW190915_235702', 'GW190924_021846', 'GW190929_012149', 'GW190930_133541', 'GW191103_012549', 'GW191105_143521', 'GW191109_010717', 'GW191113_071753', 'GW191126_115259', 'GW191127_050227', 'GW191129_134029', 'GW191204_110529', 'GW191204_171526', 'GW191215_223052', 'GW191216_213338', 'GW191222_033537', 'GW191230_180458', 'GW200112_155838', 'GW200128_022011', 'GW200129_065458', 'GW200202_154313', 'GW200208_130117', 'GW200208_222617', 'GW200209_085452', 'GW200210_092255', 'GW200216_220804', 'GW200219_094415', 'GW200220_061928', 'GW200220_124850', 'GW200224_222234', 'GW200225_060421', 'GW200302_015811', 'GW200306_093714', 'GW200308_173609', 'GW200311_115853', 'GW200316_215756', 'GW200322_091133']







# INITIALIZE
INITIALIZE_GENERAL = False # True #False #True #False#True #False
INITIALIZE_lightestBHfirst = False #True
INITIALIZE_MRR_FormationChannels = False
INITIALIZE_MRR_Spins = False
INITIALIZE_runMRR_nonMRR_ratio = False
INITIALIZE_FormationChannels_perGW = True
# 
spin_threshold=0.05





if INITIALIZE_FormationChannels_perGW==True:
    for GWname in BBHlist_GWTC_123:
        print('initializing ', GWname)
        # initialize_CSV_files_general(DCOname='BHNS')
        initialize_CSV_files_GW_formation_channel(DCOname='BHBH',  GWname=GWname)











if INITIALIZE_GENERAL==True:
    # initialize_CSV_files_general(DCOname='BHNS')
    initialize_CSV_files_general(DCOname='BHBH')
    # initialize_CSV_files_general(DCOname='NSNS')


if INITIALIZE_runMRR_nonMRR_ratio==True:
    for GWname in BBHlist_GWTC_123:
        print('initializing ', GWname)
        # initialize_CSV_files_general(DCOname='BHNS')
        initialize_CSV_files_MRR_nonMRR_ratio(DCOname='BHBH', spin_threshold=spin_threshold, GWname=GWname)



if INITIALIZE_lightestBHfirst==True:
    # initialize_CSV_files_general(DCOname='BHNS')
    initialize_CSV_files_lightestBHfirst(DCOname='BHBH')
    # initialize_CSV_files_general(DCOname='NSNS')

if INITIALIZE_MRR_FormationChannels==True:
    # initialize MRR formation channel file 
    initialize_CSV_files_MRRformationChannels(DCOname='BHBH')


if INITIALIZE_MRR_Spins==True:
    # initialize MRR formation channel file 
    initialize_CSV_files_MRRspins(DCOname='BHBH', spin_threshold=spin_threshold)

#### RUN different simulation summaries : 
runMejecta = False 
runFormationChannels =False 
runNSBH = False
runGeneralBHNS = False
runGeneralBHBH = False
runGeneralNSNS = False


runLightestFormsFirst=False
runMRR_FormationChannels = False
runMRR_Spins = False
runMRR_nonMRR_ratio = False

runFormationChannelsPerGW = True



if runFormationChannelsPerGW==True:
    for BPS in [ 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T' ]:
    # for BPS in ['C']:
    # for BPS in ['N', 'O', 'P', 'Q', 'R', 'S', 'T' ]:
        print(BPS)
        for DCOtype in ['BBH']:
            print('at DCOtype =', DCOtype)
            # for GWname in GW_MRR_list:
            
            # print('at GW =', GWname)
            writeToRatesFile_GW_FormationCHannels(BPSmodelName=BPS, DCOtype=DCOtype, GWnameList=BBHlist_GWTC_123)

            print('done with ', BPS)
            print('----------------')
            print()






if runMRR_nonMRR_ratio==True:
    for BPS in [ 'C',  'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T' ]:
    # for BPS in ['C']:
    # for BPS in ['N', 'O', 'P', 'Q', 'R', 'S', 'T' ]:
        print(BPS)
        for DCOtype in ['BBH']:
            print('at DCOtype =', DCOtype)
            # for GWname in GW_MRR_list:
            
            # print('at GW =', GWname)
            writeToRatesFile_MRR_nonMRR_ratio(BPSmodelName=BPS, DCOtype=DCOtype, spin_threshold=spin_threshold, GWnameList=BBHlist_GWTC_123)

            print('done with ', BPS)
            print('----------------')
            print()





if runMRR_Spins==True:
    for BPS in ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T' ]:
        print(BPS)
        for DCOtype in ['BBH']:
            writeToRatesFile_MRR_Spins(BPSmodelName=BPS, DCOtype=DCOtype, spin_threshold=spin_threshold)



if runMRR_FormationChannels==True:
    for BPS in ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T' ]:
        print(BPS)
        for DCOtype in ['BBH']:
            print('at DCOtype =', DCOtype)
            writeToRatesFile_MRR_FormationChannels(BPSmodelName=BPS, DCOtype=DCOtype)
            print('done with ', BPS)




# if runMRR_FormationChannels==True:
#   for BPS in ['F', 'H', 'K']:
#       print(BPS)
#       for DCOtype in ['BBH']:
#           print('at DCOtype =', DCOtype)
#           writeToRatesFile_MRR_FormationChannels(BPSmodelName=BPS, DCOtype=DCOtype)
#           print('done with ', BPS)



if runLightestFormsFirst==True:
    for BPS in ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T' ]:
        print(BPS)
        for DCOtype in ['BBH']:
            print('at DCOtype =', DCOtype)
            writeToRatesFile_lightestFormsFirst(BPSmodelName=BPS, DCOtype=DCOtype)
            print('done with ', BPS)



if runGeneralBHBH ==True:
    for BPS in  ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T' ]:
        print(BPS)
        for DCOtype in ['BBH']:
            print('at DCOtype =', DCOtype)
            writeToRatesFile_GENERAL(BPSmodelName=BPS, DCOtype=DCOtype)
            print('done with ', BPS)

    # # INITIALIZE FILE 
    # namesEMlist=[]

    # DCOname ='BHBH'
    # iii=0
    

    # # CREATE PANDAS FILE 
    # nModels=26
    # BPSnameslist = list(string.ascii_uppercase)[0:nModels]

    # NAMES = []
    # stringgg = 'lightestFormsFirst'

    # for ind_l, L in enumerate(BPSnameslist):
    #   str_z0 = str(L + ' intrinsic (z=0) [Gpc^{-3} yr^{-1}]')
    #   str_obs = str(L + ' observed (design LVK) [yr^{-1}]')
    #   NAMES.append(str_z0)
    #   NAMES.append(str_obs)
        
        


    # datas=[]

    # for i in range(len(BPSnameslist)):
    #   datas.append(np.zeros_like(MSSFRnameslist))
    #   datas.append(np.zeros_like(MSSFRnameslist))
        
        
    # df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
    # df.columns =   df.columns.map(str)
    # df.index.names = ['.x.y.z']
    # df.columns.names = ['m']

    # # print(df) 

    # df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/MRR_Project/dataFiles/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv')

    # # run calculation 
    # for BPS in   ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q' ]:
    #   print(BPS)
    #   for DCOtype in ['BBH']:
    #       print('at DCOtype =', DCOtype)
    #       writeToRatesFile_lightestFormsFirst(BPSmodelName=BPS, DCOtype=DCOtype)
    #       print('done with ', BPS)














# if runMejecta ==True:
#   for BPS in  ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q' ]:
#       print(BPS)
#       for DCOtype in ['BHNS']:
#           print('at DCOtype =', DCOtype)
#           writeToRatesFile_Mejecta(BPSmodelName=BPS)
#           print('done with ', BPS)




# if runFormationChannels ==True:
#   # for BPS in  ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O' ]:
#   #   print(BPS)
#   #   for DCOtype in ['BNS']:
#   #       print('at DCOtype =', DCOtype)
#   #       writeToRatesFile_FormationChannels(BPSmodelName=BPS, DCOtype='BNS')
#   #       print('done with ', BPS)


#   for BPS in  ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q' ]:
#       print(BPS)
#       for DCOtype in ['BBH']:
#           print('at DCOtype =', DCOtype)
#           writeToRatesFile_FormationChannels(BPSmodelName=BPS, DCOtype='BBH')
#           print('done with ', BPS)

# if runNSBH ==True:
#   for BPS in ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q' ]:
#       print(BPS)
#       for DCOtype in ['BHNS']:
#           print('at DCOtype =', DCOtype)
#           writeToRatesFile_NSBH(BPSmodelName=BPS)
#           print('done with ', BPS)




# if runGeneralBHNS ==True:
#   for BPS in ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q' , 'R', 'S', 'T']:
#       print(BPS)
#       for DCOtype in ['BHNS']:
#           print('at DCOtype =', DCOtype)
#           writeToRatesFile_GENERAL(BPSmodelName=BPS, DCOtype=DCOtype)
#           print('done with ', BPS)


# if runGeneralNSNS ==True:
#   for BPS in  ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T' ]:
#       print(BPS)
#       for DCOtype in ['BBH']:
#           print('at DCOtype =', DCOtype)
#           writeToRatesFile_GENERAL(BPSmodelName=BPS, DCOtype=DCOtype)
#           print('done with ', BPS)

# if runGeneralBHBH ==True:
#   for BPS in  ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T' ]:
#       print(BPS)
#       for DCOtype in ['BNS']:
#           print('at DCOtype =', DCOtype)
#           writeToRatesFile_GENERAL(BPSmodelName=BPS, DCOtype=DCOtype)
#           print('done with ', BPS)














# Models to RUN 

# May 20: I am updating my data with the AllDCO focused runs :-) 

# this is an overwrite with better data (old ones are in BHNS copy)


# for DCOtype in ['BHNS', 'BBH', 'BNS']:
#   print('at DCOtype =', DCOtype)
#   pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/fiducial/'
#   modelname = 'A'
#   writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)


#   pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/fiducial/'
#   modelname = 'B'
#   writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=True)


#   pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/zeroBHkick/'
#   modelname = 'G'
#   writeToRatesFile(modelname=modelname, pa thCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=True)

# INITIALIZE_FormationChannels = False
# INITIALIZE_NSBH= False #False#True
# INITIALIZE=False #False #True 
# INITIALIZE_GW190814 = False
# INITIALIZE_EM =False

# # INITIALIZE_NSBH= True #False#True
# # INITIALIZE=True #False #True 
# # INITIALIZE_GENERAL = True #False#True #False


# # ['.0.0.0', '.1.1.1', '.2.1.1', '.3.1.1', '.1.1.2', '.2.1.2', '.3.1.2', '.1.1.3', '.2.1.3', '.3.1.3', '.1.2.1', '.2.2.1', '.3.2.1', '.1.2.2', '.2.2.2', '.3.2.2', '.1.2.3', '.2.2.3', '.3.2.3', '.1.3.1', '.2.3.1', '.3.3.1', '.1.3.2', '.2.3.2', '.3.3.2', '.1.3.3', '.2.3.3', '.3.3.3']


# if INITIALIZE_FormationChannels==True:

#   # namesEMlist=[]

#   DCOname ='NSNS' 
#   iii=0
    

#   # CREATE PANDAS FILE 
#   nModels=26
#   BPSnameslist = list(string.ascii_uppercase)[0:nModels]

#   NAMES = []
#   stringgg =  'AllDCOsimulation_formation_channels'

#   for ind_l, BPSmodelName in enumerate(BPSnameslist):
#       # str_z0 = str(L + ' intrinsic (z=0) [Gpc^{-3} yr^{-1}]')
#       # str_obs = str(L + ' observed (design LVK) [yr^{-1}]')

#       namez0 = BPSmodelName + 'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#       nameObs = BPSmodelName + 'All observed (design LVK) [yr^{-1}]'

#       namez0_I = BPSmodelName + 'channel I intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#       nameObs_I = BPSmodelName + 'channel I observed (design LVK) [yr^{-1}]'
#       namez0_II = BPSmodelName + 'channel II intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#       nameObs_II = BPSmodelName + 'channel II observed (design LVK) [yr^{-1}]'
#       namez0_III = BPSmodelName + 'channel III intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#       nameObs_III = BPSmodelName + 'channel III observed (design LVK) [yr^{-1}]'
#       namez0_IV = BPSmodelName + 'channel IV intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#       nameObs_IV = BPSmodelName + 'channel IV observed (design LVK) [yr^{-1}]'
#       namez0_V = BPSmodelName + 'channel V intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#       nameObs_V = BPSmodelName + 'channel V observed (design LVK) [yr^{-1}]'

#       NAMES.append(namez0)
#       NAMES.append(nameObs)

#       NAMES.append(namez0_I)
#       NAMES.append(nameObs_I)
#       NAMES.append(namez0_II)
#       NAMES.append(nameObs_II)
#       NAMES.append(namez0_III)
#       NAMES.append(nameObs_III)
#       NAMES.append(namez0_IV)
#       NAMES.append(nameObs_IV)
#       NAMES.append(namez0_V)
#       NAMES.append(nameObs_V)


        
        


#   datas=[]

#   for i in range(len(BPSnameslist)):
#       for ii in range(6):
#           datas.append(np.zeros_like(MSSFRnameslist))
#           datas.append(np.zeros_like(MSSFRnameslist))
        
        
#   df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
#   df.columns =   df.columns.map(str)
#   df.index.names = ['xyz']
#   df.columns.names = ['m']

#   # print(df) 

#   df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/Double-Compact-Object-Mergers/dataFiles/summary_data_Fig_2/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv')





# if INITIALIZE_GW190814==True:

#   for dcotype in ['NSNS', 'BHBH', 'BHNS']:

#       namesEMlist=[]

#       DCOname=dcotype
#       iii=0
        

#       # CREATE PANDAS FILE 
#       nModels=26
#       BPSnameslist = list(string.ascii_uppercase)[0:nModels]

#       NAMES = []
#       stringgg = 'GW190814rate'

#       for ind_l, L in enumerate(BPSnameslist):
#           str_z0 = str(L + ' intrinsic (z=0) [Gpc^{-3} yr^{-1}]')
#           str_obs = str(L + ' observed (design LVK) [yr^{-1}]')
#           NAMES.append(str_z0)
#           NAMES.append(str_obs)
            
            


#       datas=[]

#       for i in range(len(BPSnameslist)):
#           datas.append(np.zeros_like(MSSFRnameslist))
#           datas.append(np.zeros_like(MSSFRnameslist))
            
            
#       df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
#       df.columns =   df.columns.map(str)
#       df.index.names = ['xyz']
#       df.columns.names = ['m']

#       # print(df) 

#       df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/Double-Compact-Object-Mergers/dataFiles/summary_data_Fig_2/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv')





# if INITIALIZE_NSBH==True:


#   namesEMlist=[]

#   DCOname ='BHNS'
#   iii=0
    

#   # CREATE PANDAS FILE 
#   nModels=26
#   BPSnameslist = list(string.ascii_uppercase)[0:nModels]

#   NAMES = []
#   stringgg = 'NSBH'

#   for ind_l, L in enumerate(BPSnameslist):
#       str_z0 = str(L + ' intrinsic (z=0) [Gpc^{-3} yr^{-1}]')
#       str_obs = str(L + ' observed (design LVK) [yr^{-1}]')
#       NAMES.append(str_z0)
#       NAMES.append(str_obs)
        
        


#   datas=[]

#   for i in range(len(BPSnameslist)):
#       datas.append(np.zeros_like(MSSFRnameslist))
#       datas.append(np.zeros_like(MSSFRnameslist))
        
        
#   df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
#   df.columns =   df.columns.map(str)
#   df.index.names = ['.x.y.z']
#   df.columns.names = ['m']

#   # print(df) 

#   df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/Double-Compact-Object-Mergers/dataFiles/summary_data_Fig_2/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv')









# # #### INITIALIZE::: 
# if INITIALIZE_EM==True:

#   namesEMlist=[]

#   DCOname ='BHNS'
#   iii=0
#   for ind_chi, chi in enumerate([0.0, .5, 'Qin']):
#       # print(chi)
#       iii+=1
#       BH_chi   = chi 
#       for ind_Rns, NSradii in enumerate([11.5,13.0]):
#           iii+=1
#           Rns = NSradii
#           # if ind_mssfr ==0:
#           #   # print(chi)
#           stringg = 'Rns_'+ str(NSradii) + 'km_' + 'spinBH_' + str(chi) 
#           namesEMlist.append(stringg)


#           # CREATE PANDAS FILE 
#           nModels=26
#           BPSnameslist = list(string.ascii_uppercase)[0:nModels]

#           NAMES = []

#           for ind_l, L in enumerate(BPSnameslist):
#               str_z0 = str(L + ' intrinsic (z=0) [Gpc^{-3} yr^{-1}]')
#               str_obs = str(L + ' observed (design LVK) [yr^{-1}]')
#               NAMES.append(str_z0)
#               NAMES.append(str_obs)
                


#           datas=[]

#           for i in range(len(BPSnameslist)):
#               datas.append(np.zeros_like(MSSFRnameslist))
#               datas.append(np.zeros_like(MSSFRnameslist))
                
#           print(MSSFRnameslist)
#           df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
#           df.columns =   df.columns.map(str)
#           df.index.names = ['xyz']
#           df.columns.names = ['m']

#           # print(df) 

#           df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/Double-Compact-Object-Mergers/dataFiles/summary_data_Fig_2/rates_MSSFR_Models_'+DCOname+ '_' + stringg + '.csv')


# # print(namesEMlist)


# # for DCOtype in ['BHNS']:

# #     for Rns in enumerate()

# #     pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/zeroBHkick/'
# #     modelname = 'G'
# #     print('modelname')
# #     writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=True)

# #     print('at DCOtype =', DCOtype)
# #     pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/fiducial/'
# #     modelname = 'A'
# #     print('modelname')
# #     writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)


# #     pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/fiducial/'
# #     modelname = 'B'
# #     print('modelname')
# #     writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=True)




# # for DCOtype in ['BHNS']:

# #     for Rns in enumerate()

# #     pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/zeroBHkick/'
# #     modelname = 'G'
# #     print('modelname')
# #     writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=True)

# #     print('at DCOtype =', DCOtype)
# #     pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/fiducial/'
# #     modelname = 'A'
# #     print('modelname')
# #     writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)


# #     pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/fiducial/'
# #     modelname = 'B'
# #     print('modelname')
# #     writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=True)





# # pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/alpha0_5/'
# # modelname, DCOtype = 'M', 'BNS'
# # writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)

# # DCOtype='BHNS'
# # writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)

# # DCOtype='BBH'
# # writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)



# # pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/alpha2_0/'
# # modelname, DCOtype = 'N', 'BNS'
# # writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)

# # DCOtype='BHNS'
# # writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)

# # DCOtype='BBH'
# # writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)



