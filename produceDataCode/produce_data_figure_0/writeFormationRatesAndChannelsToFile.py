
# from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import time
import sys
import copy
#Quick fudge to make import from ../Scripts work

sys.path.append('../Scripts')


# import gc


# import ClassCosmicIntegrator  as CI #Given settings and redshifts returns rates (2D arrays) Loads the data
# import coencodeVarious        as CV
from PostProcessingScripts import * 
from formation_channels import * 
import ClassCOMPAS     as CC ###


import pandas as pd
from astropy import units as u
from astropy import constants as const



##### I THINK THIS IS AN OLD FILE #### 



dictDCOtypeDCOlabel = {'BBH':'BHBH', 'BNS':'NSNS', 'BHNS':'BHNS'}





# # seedsChannelSCCE, percentageChannelSCCE = returnSeedsPercentageOther()

def createEmptyCSVplaceholder(DCOtype='BBH', nBPSmodels=15):


   


    DCOname=dictDCOtypeDCOlabel[DCOtype]
    BPSnameslist = list(string.ascii_uppercase)[0:nBPSmodels]   
    channel_names = ['total', 'I_classic', 'II_only_stable_MT', 'III_single_core_CE', 'IV_double_core_CE', 'V_other']



    NAMES = []
    # stringgg = 'GW190814rate'

    for ind_m, m_ in enumerate(BPSnameslist):
        for ind_c, c_ in enumerate(channel_names):
            str_ = m_ + ' ' + c_ + '  [Msun^{-1}]'

            NAMES.append(str_)

            
            


    datas=[]
    nMetallicities = 53
    Zlist=['0_0001','0_00011', '0_00012', '0_00014', '0_00016', '0_00017',\
    '0_00019', '0_00022', '0_00024', '0_00027', '0_0003', '0_00034',\
    '0_00037', '0_00042', '0_00047', '0_00052', '0_00058', '0_00065',\
    '0_00073', '0_00081', '0_0009', '0_00101', '0_00113', '0_00126',\
    '0_0014', '0_00157', '0_00175', '0_00195', '0_00218', '0_00243',\
    '0_00272', '0_00303', '0_00339', '0_00378', '0_00422', '0_00471',\
    '0_00526', '0_00587', '0_00655', '0_00732', '0_00817', '0_00912',\
    '0_01018', '0_01137', '0_01269', '0_01416', '0_01581', '0_01765', '0_01971', '0_022', '0_0244', '0_02705', '0_03']

    for i in range(len(NAMES)):
        datas.append(np.zeros(nMetallicities))
        # datas.append(np.zeros(nMetallicities))
        
            
    df = pd.DataFrame(data=datas, index=NAMES, columns=Zlist).T
    df.columns =   df.columns.map(str)
    df.index.names = ['Z_i']
    df.columns.names = ['model']

        # print(df) 

    df.to_csv('formationRatesTotalAndPerChannel_'+DCOname+ '_' +  '.csv')
    return 


INITIALIZE=False 

if INITIALIZE == True:
    createEmptyCSVplaceholder(DCOtype='BNS', nBPSmodels=15)
    createEmptyCSVplaceholder(DCOtype='BNS', nBPSmodels=15)
    createEmptyCSVplaceholder(DCOtype='BHNS', nBPSmodels=15)



def writeFormationRatesAndChannelsToFile(DCOtype='BBH', \
    pathCOMPASOutput='/Volumes/Andromeda/DATA/AllDCO_bugfix/',\
     alphabetDirDict=[], nBPSmodels=15):
    
    
    BPSnameslist = list(string.ascii_uppercase)[0:nBPSmodels]   
    channel_names = ['total', 'I_classic', 'II_only_stable_MT', 'III_single_core_CE', 'IV_double_core_CE', 'V_other']
    # temp = range(nModels+3)
    DCOname=dictDCOtypeDCOlabel[DCOtype]
    
 

    print('now at DCO type  ', DCOtype)
        
    for ind_m, bps_model in enumerate(BPSnameslist):    

        print()
        print('now at model ', alphabetDirDict[bps_model])
            
        # set always optimistic CE false, unless we are doing the optimistic variation
        OPTIMISTIC=False
        if bps_model=='H':
            OPTIMISTIC=True
            print('doing optimistic version of fiducial')
            
        # path to datafile 
        path = pathCOMPASOutput+alphabetDirDict[bps_model] + '/'

            
        #But I want only within Hubble time 
        Data            = CC.COMPASData(path=path, lazyData=True, Mlower=5., \
                         Mupper=150, binaryFraction=1)
        Data.setCOMPASDCOmask(types=DCOtype,  withinHubbleTime=True, optimistic=OPTIMISTIC)
        Data.setCOMPASData()
        
        metallicities = Data.metallicitySystems
        seeds    = Data.seeds[Data.Hubble==True]
        weights = Data.weight
            
            
        
       

        

        seedsPercentageClassic, seedsPercentageOnlyStableMT = returnSeedsPercentageClassicAndOnlyStableMT(pathCOMPASOutput=path,\
                                        types=DCOtype,  withinHubbleTime=True, optimistic=OPTIMISTIC, \
                                        binaryFraction=1)
        seedsClassic, percentageClassic = seedsPercentageClassic
        seedsOnlyStableMT, percentageOnlyStableMT = seedsPercentageOnlyStableMT



        seedsDoubleCE, percentageDoubleCE = returnSeedsPercentageDoubleCoreCEE(pathCOMPASOutput=path,\
                                        types=DCOtype,  withinHubbleTime=True, optimistic=OPTIMISTIC, \
                                        binaryFraction=1)


        seedsSingleCE, percentageSingleCE = returnSeedsPercentageSingleCoreCEE(pathCOMPASOutput=path,\
                                        types=DCOtype,  withinHubbleTime=True, optimistic=OPTIMISTIC, \
                                        binaryFraction=1)



        seedschannels = [seedsClassic, seedsOnlyStableMT, seedsSingleCE, seedsDoubleCE]


        seedsOther, percentageOther = returnSeedsPercentageOther(pathCOMPASOutput=path,\
                                        types=DCOtype,  withinHubbleTime=True, optimistic=OPTIMISTIC, \
                                        binaryFraction=1, channelsSeedsList=seedschannels)


        seedschannels = [seedsClassic, seedsOnlyStableMT, seedsSingleCE, seedsDoubleCE, seedsOther]




        dictChannelsBHNS = { 'classic':seedsClassic, \
                            'immediate CE':seedsSingleCE,\
                                 'stable B no CEE':seedsOnlyStableMT, \
                             r'double-core CE':seedsDoubleCE,  \
                                'other':seedsOther\
                               }



        formationRateTotal           = np.zeros(len(Data.metallicityGrid))  
        formationRateClassic         = np.zeros(len(Data.metallicityGrid)) 
        formationRateOnlyStableMT    = np.zeros(len(Data.metallicityGrid)) 
        formationRateSingleCE        = np.zeros(len(Data.metallicityGrid)) 
        formationRateDoubleCE        = np.zeros(len(Data.metallicityGrid)) 
        formationRateOther           = np.zeros(len(Data.metallicityGrid)) 


        for nrZ, Z in enumerate(Data.metallicityGrid):
            maskZ = (metallicities == Z)
            formationRateTotal[nrZ] = np.sum(weights[maskZ]) # //floor weights


            # mask different channels
            InClassic       = np.in1d(seeds, seedsClassic)
            InOnlyStableMT  = np.in1d(seeds, seedsOnlyStableMT)
            InSingleCE      = np.in1d(seeds, seedsSingleCE)
            InDoubleCE      = np.in1d(seeds, seedsDoubleCE)
            InOther         = np.in1d(seeds, seedsOther)

            maskClassic         = (metallicities == Z) & (InClassic==1)
            maskOnlyStableMT    = (metallicities == Z) & (InOnlyStableMT==1)
            maskSingleCE        = (metallicities == Z) & (InSingleCE==1)
            maskDoubleCE        = (metallicities == Z) & (InDoubleCE==1)
            maskOther           = (metallicities == Z) & (InOther==1)

            formationRateClassic[nrZ]         = np.sum(weights[maskClassic])
            formationRateOnlyStableMT[nrZ]    = np.sum(weights[maskOnlyStableMT])
            formationRateSingleCE[nrZ]        = np.sum(weights[maskSingleCE]) 
            formationRateDoubleCE[nrZ]        = np.sum(weights[maskDoubleCE])
            formationRateOther[nrZ]           = np.sum(weights[maskOther])

        formationRateTotal = np.divide(formationRateTotal, Data.totalMassEvolvedPerZ) + 0 #lowerY        
        formationRateClassic = np.divide(formationRateClassic, Data.totalMassEvolvedPerZ)
        formationRateOnlyStableMT = np.divide(formationRateOnlyStableMT, Data.totalMassEvolvedPerZ)
        formationRateSingleCE = np.divide(formationRateSingleCE, Data.totalMassEvolvedPerZ)
        formationRateDoubleCE = np.divide(formationRateDoubleCE, Data.totalMassEvolvedPerZ)
        formationRateOther = np.divide(formationRateOther, Data.totalMassEvolvedPerZ)


        df = pd.read_csv('formationRatesTotalAndPerChannel_'+DCOname+ '_' +  '.csv', index_col=0)
        # namez0 = bps_model +' total  [Msun^{-1}]'
        for ind_c, c_ in enumerate(channel_names):
            str_ = bps_model + ' ' + c_ + '  [Msun^{-1}]'

            # total rates 
            if c_=='total':
                df[str_] = formationRateTotal 
            elif c_=='I_classic':
                df[str_] = formationRateClassic
            elif c_=='II_only_stable_MT':
                df[str_] = formationRateOnlyStableMT
            elif c_=='III_single_core_CE':
                df[str_] = formationRateSingleCE
            elif c_=='IV_double_core_CE':
                df[str_] = formationRateDoubleCE
            elif c_=='V_other':
                df[str_] = formationRateOther


        df.to_csv('formationRatesTotalAndPerChannel_'+DCOname+ '_' +  '.csv')


    print('finished')

    return


import string
nModels=15
BPSnameslist = list(string.ascii_uppercase)[0:nModels]
modelDirList = ['fiducial', 'massTransferEfficiencyFixed_0_25', 'massTransferEfficiencyFixed_0_5', 'massTransferEfficiencyFixed_0_75', \
               'unstableCaseBB', 'alpha0_5', 'alpha2_0', 'fiducial', 'rapid', 'maxNSmass2_0', 'maxNSmass3_0', 'noPISN',  'ccSNkick_100km_s', 'ccSNkick_30km_s', 'noBHkick' ]

alphabetDirDict =  {BPSnameslist[i]: modelDirList[i] for i in range(len(BPSnameslist))}


writeFormationRatesAndChannelsToFile(DCOtype='BHNS', \
    pathCOMPASOutput='/Volumes/Andromeda/DATA/AllDCO_bugfix/',\
     alphabetDirDict=alphabetDirDict, nBPSmodels=nModels)


writeFormationRatesAndChannelsToFile(DCOtype='BNS', \
    pathCOMPASOutput='/Volumes/Andromeda/DATA/AllDCO_bugfix/',\
     alphabetDirDict=alphabetDirDict, nBPSmodels=nModels)

writeFormationRatesAndChannelsToFile(DCOtype='BBH', \
    pathCOMPASOutput='/Volumes/Andromeda/DATA/AllDCO_bugfix/',\
     alphabetDirDict=alphabetDirDict, nBPSmodels=nModels)

# plotFormationRatePerZ(pathCOMPASOutput='/Volumes/Andromeda/DATA/AllDCO_bugfix/', alphabetDirDict=alphabetDirDict)    
    
    
