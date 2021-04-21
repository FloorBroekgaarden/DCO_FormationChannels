import numpy as np
import sys
import h5py as h5

#Quick fudge to make import from ../Scripts work
sys.path.append('../../Scripts')


# for reading datafiles 
import pandas as pd

# import script that has formation channel classification functions:
from PostProcessingScripts import * 
from formation_channels import * 
import ClassCOMPAS     as CC ###


def createEmptyCSVplaceholder(DCOtype='BBH'):


   


    DCOname=DCOname_dict[DCOtype]
    channel_names = ['total', 'I_classic', 'II_only_stable_MT', 'III_single_core_CE', 'IV_double_core_CE', 'V_other']



    NAMES = []
    # stringgg = 'GW190814rate'
    NAMES.append('Z_i')
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



    Zinitial = [0.0001, 0.00011, 0.00012, 0.00014, 0.00016, 0.00017,\
       0.00019, 0.00022, 0.00024, 0.00027, 0.0003, 0.00034, \
       0.00037, 0.00042, 0.00047, 0.00052, 0.00058, 0.00065,\
       0.00073, 0.00081, 0.0009, 0.00101, 0.00113, 0.00126,\
       0.0014, 0.00157, 0.00175, 0.00195, 0.00218, 0.00243, \
       0.00272, 0.00303, 0.00339, 0.00378, 0.00422, 0.00471, \
       0.00526, 0.00587, 0.00655, 0.00732, 0.00817, 0.00912, \
       0.01018, 0.01137, 0.01269, 0.01416, 0.01581, 0.01765, 0.01971, 0.022, 0.0244, 0.02705, 0.03]

    datas.append(Zinitial)

    for i in range(len(NAMES)-1):
        datas.append(np.zeros(nMetallicities))
        # datas.append(np.zeros(nMetallicities))
        
            
    df = pd.DataFrame(data=datas, index=NAMES, columns=Zlist).T
    df.columns =   df.columns.map(str)
    df.index.names = ['Z_i_label']
    df.columns.names = ['model']

        # print(df) 

    df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_0/formationRatesTotalAndPerChannel_'+DCOname+    '.csv')
    return 











def writeFormationRatesAndChannelsToFile(DCOtype='BBH', \
    pathCOMPASOutput='/Volumes/Andromeda/DATA/AllDCO_bugfix/'):
    
    
  
    channel_names = ['total', 'I_classic', 'II_only_stable_MT', 'III_single_core_CE', 'IV_double_core_CE', 'V_other']
    # temp = range(nModels+3)
    DCOname=DCOname_dict[DCOtype]
    
 

    print('now at DCO type  ', DCOtype)

    writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_0/formationRatesTotalAndPerChannel_' + DCOname + '_.csv'     

    Zinitial = [0.0001, 0.00011, 0.00012, 0.00014, 0.00016, 0.00017,\
       0.00019, 0.00022, 0.00024, 0.00027, 0.0003, 0.00034, \
       0.00037, 0.00042, 0.00047, 0.00052, 0.00058, 0.00065,\
       0.00073, 0.00081, 0.0009, 0.00101, 0.00113, 0.00126,\
       0.0014, 0.00157, 0.00175, 0.00195, 0.00218, 0.00243, \
       0.00272, 0.00303, 0.00339, 0.00378, 0.00422, 0.00471, \
       0.00526, 0.00587, 0.00655, 0.00732, 0.00817, 0.00912, \
       0.01018, 0.01137, 0.01269, 0.01416, 0.01581, 0.01765, 0.01971, 0.022, 0.0244, 0.02705, 0.03]
        
    for ind_m, bps_model in enumerate(BPSnameslist):    

        print()
        print('now at model ', bps_model)

        OPTIMISTIC=False
        if bps_model=='H':
            OPTIMISTIC=True
            print('doing optimistic version of fiducial')



        path_dir = '/Volumes/Andromeda/DATA/AllDCO_bugfix/'
        path_ = path_dir
        path_ = path_ + alphabetDirDict[bps_model] +'/'
        path  = path_ + 'COMPASCompactOutput_'+ DCOtype +'_' + bps_model + '.h5'


        # read in data 
        fdata = h5.File(path) 


        #But I want only within Hubble time 
        Data            = CC.COMPASData(path=path, lazyData=True, Mlower=5., \
                         Mupper=150, binaryFraction=1)
        Data.setCOMPASDCOmask(types=DCOtype,  withinHubbleTime=True, optimistic=OPTIMISTIC)
        Data.setCOMPASData()
        totalMassEvolvedPerZ = Data.totalMassEvolvedPerZ
        del Data

        seeds = fdata['doubleCompactObjects']['seed'][...].squeeze()
        channels = identify_formation_channels(seeds=seeds, file=fdata)
        
        metallicities = fdata['doubleCompactObjects']['Metallicity1'][...].squeeze()

        # get intrinsic weights
        weights = fdata['doubleCompactObjects']['weight'][...].squeeze()




        formationRateTotal           = np.zeros(len(Zinitial))  
        formationRateClassic         = np.zeros(len(Zinitial)) 
        formationRateOnlyStableMT    = np.zeros(len(Zinitial)) 
        formationRateSingleCE        = np.zeros(len(Zinitial)) 
        formationRateDoubleCE        = np.zeros(len(Zinitial)) 
        formationRateOther           = np.zeros(len(Zinitial)) 

        # to keep track of where we are in totalMassEvolvedPerZ, since the length of this array only takes into account the Z where there are systems
        ZinTotalMassEvolvedPerZ_ind =0
        for nrZ, Z in enumerate(Zinitial):
            maskZ = (metallicities == Z)
            

            # if we have systems at this metallicity do the following (calculate the formation yield)
            if (np.sum(maskZ))>0:
                 
                # mask different channels
                maskClassic         = (metallicities == Z) & (channels==1)
                maskOnlyStableMT    = (metallicities == Z) & (channels==2)
                maskSingleCE        = (metallicities == Z) & (channels==3)
                maskDoubleCE        = (metallicities == Z) & (channels==4)
                maskOther           = (metallicities == Z) & (channels==0)

                formationRateTotal[nrZ]           = np.sum(weights[maskZ]) / totalMassEvolvedPerZ[ZinTotalMassEvolvedPerZ_ind] # TOTAL rates 
                formationRateClassic[nrZ]         = np.sum(weights[maskClassic]) / totalMassEvolvedPerZ[ZinTotalMassEvolvedPerZ_ind]
                formationRateOnlyStableMT[nrZ]    = np.sum(weights[maskOnlyStableMT]) / totalMassEvolvedPerZ[ZinTotalMassEvolvedPerZ_ind]
                formationRateSingleCE[nrZ]        = np.sum(weights[maskSingleCE])  / totalMassEvolvedPerZ[ZinTotalMassEvolvedPerZ_ind]
                formationRateDoubleCE[nrZ]        = np.sum(weights[maskDoubleCE]) / totalMassEvolvedPerZ[ZinTotalMassEvolvedPerZ_ind]
                formationRateOther[nrZ]           = np.sum(weights[maskOther]) / totalMassEvolvedPerZ[ZinTotalMassEvolvedPerZ_ind]
                ZinTotalMassEvolvedPerZ_ind +=1

            # otherwise add 0, 
            else:
                print('no systems at metallicity ', Z)
                formationRateTotal[nrZ]           = 0
                formationRateClassic[nrZ]         = 0
                formationRateOnlyStableMT[nrZ]    = 0
                formationRateSingleCE[nrZ]        = 0
                formationRateDoubleCE[nrZ]        = 0
                formationRateOther[nrZ]           = 0              


        # formationRateTotal        = np.divide(formationRateTotal, totalMassEvolvedPerZ) + 0 #lowerY        
        # formationRateClassic      = np.divide(formationRateClassic, totalMassEvolvedPerZ)
        # formationRateOnlyStableMT = np.divide(formationRateOnlyStableMT,totalMassEvolvedPerZ)
        # formationRateSingleCE     = np.divide(formationRateSingleCE, totalMassEvolvedPerZ)
        # formationRateDoubleCE     = np.divide(formationRateDoubleCE, totalMassEvolvedPerZ)
        # formationRateOther        = np.divide(formationRateOther, totalMassEvolvedPerZ)


        df = pd.read_csv('/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_0/formationRatesTotalAndPerChannel_'+DCOname+    '.csv', index_col=0)
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


        df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_0/formationRatesTotalAndPerChannel_'+DCOname+    '.csv')


    print('finished')

    return



#######################################################
#
#  The code below first initializes BNS, BBH and BHNS empty data files through "INITIALIZE = TRUE"
#  Then it goes for BNS, BHNS and BBH through the hdf5 files, per metallicity it calculates which formation channel each DCO type formed through, 
#  calculates its formation rate per Solar mass evolved, and writes that into a data file. 
#  You will end up with a datafile for each DCO type, that for each BPS model (A, B, C, ... ) lists for each metallicity the formation channel contribution (total rate per mass evolved)
#
###############################################


INITIALIZE=True

if INITIALIZE == True:
    createEmptyCSVplaceholder(DCOtype='BNS')
    createEmptyCSVplaceholder(DCOtype='BBH')
    createEmptyCSVplaceholder(DCOtype='BHNS')



writeFormationRatesAndChannelsToFile(DCOtype='BNS', \
    pathCOMPASOutput='/Volumes/Andromeda/DATA/AllDCO_bugfix/')


writeFormationRatesAndChannelsToFile(DCOtype='BHNS', \
    pathCOMPASOutput='/Volumes/Andromeda/DATA/AllDCO_bugfix/')


writeFormationRatesAndChannelsToFile(DCOtype='BBH', \
    pathCOMPASOutput='/Volumes/Andromeda/DATA/AllDCO_bugfix/')



