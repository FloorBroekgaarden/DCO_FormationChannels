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








def writeToRatesFile_FormationChannels(BPSmodelName='Z', DCOtype='BHNS'):



    DCOname = DCOname_dict[DCOtype]

    # path for files 
    path_dir = '/Volumes/Andromeda2/DATA/AllDCO_bugfix/'
    path_ = path_dir
    path_ = path_ + alphabetDirDict[BPSmodelName] +'/'
    path  = path_ + 'COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'




    # read in data 
    fdata = h5.File(path)



    # set optimistic true if that is the variation (H) 
    # OPTIMISTIC=False
    # if (bps_model=='F') or (bps_model=='K'):
    #     OPTIMISTIC=True
    #     print('doing optimistic version of %s'%alphabetDirDict[bps_model])


    seeds = fdata['doubleCompactObjects']['seed'][...].squeeze()
    channels = identify_formation_channels(seeds=seeds, file=fdata)

    
    
    

    stringgg =  'formation_channels'
    writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_1/Formation_yields_'  + stringgg + '_'  + DCOname + 'fast.csv'     
    headerDict_intrinsic = {5:'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  0:'channel V intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  1:'channel I intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 2:'channel II intrinsic (z=0) [Gpc^{-3} yr^{-1}]',3:'channel III intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 4:'channel IV intrinsic (z=0) [Gpc^{-3} yr^{-1}]'}
    
    headerDict_observed  = {5:'All observed (design LVK) [yr^{-1}]',  0:'channel V observed (design LVK) [yr^{-1}]', 1:'channel I observed (design LVK) [yr^{-1}]', 2:'channel II observed (design LVK) [yr^{-1}]', 3:'channel III observed (design LVK) [yr^{-1}]', 4:'channel IV observed (design LVK) [yr^{-1}]'}    


    
    # get intrinsic weights
    fparam_intrinsic = 'weights_intrinsic'
    # get detected weights
    fparam_detected = 'weights_detected'


    ####################################################
    ######### ITERATE  OVER  MSSFR  MODELS #############
    ####################################################





    df = pd.read_csv(writePath, index_col=0)
    

    for nrC, Channel in enumerate(range(6)):          

    #           #Get the seeds that relate to sorted indices
        mask_C  = (channels==Channel)
        
    
        intrinsicRates = np.zeros(len(MSSFRnameslist))
        detectedRates = np.zeros(len(MSSFRnameslist))       


        for ind_mssfr, mssfr in enumerate(MSSFRnameslist):


            weightheader = 'w_' + mssfr
            w_int = fdata[fparam_intrinsic][weightheader][...].squeeze()
            w_det = fdata[fparam_detected][weightheader][...].squeeze()



            if nrC==5: 
                # TOTAL RATE
                intrinsicRates[ind_mssfr] = np.sum(w_int)
                detectedRates[ind_mssfr]  = np.sum(w_det)  
            else:
                # CHANNEL RATE 
                intrinsicRates[ind_mssfr] = np.sum(w_int[mask_C])
                detectedRates[ind_mssfr]  = np.sum(w_det[mask_C])          

            

            
        
        namez0  = BPSmodelName + '_' + headerDict_intrinsic[Channel]
        nameObs = BPSmodelName + '_' + headerDict_observed[Channel]
        df[namez0] = intrinsicRates
        df[nameObs] = detectedRates



    df.to_csv(writePath)


    fdata.close() 

    return



def initalize_formationChannels(DCOname):

    stringgg =  'formation_channels'
    writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_1/Formation_yields_'  + stringgg + '_'  + DCOname + 'fast.csv' 

    iii=0


    # CREATE PANDAS FILE 
    nModels=26
    BPSnameslist = list(string.ascii_uppercase)[0:nModels]

    NAMES = []
    stringgg =  'AllDCOsimulation_formation_channels'
    



    for ind_l, BPSmodelName in enumerate(BPSnameslist):
        # str_z0 = str(L + ' intrinsic (z=0) [Gpc^{-3} yr^{-1}]')
        # str_obs = str(L + ' observed (design LVK) [yr^{-1}]')

        namez0 = BPSmodelName + '_' +'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
        nameObs = BPSmodelName + '_' +'All observed (design LVK) [yr^{-1}]'


        # namez0_oth = BPSmodelName + '_' +'Other intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
        # nameObs_oth = BPSmodelName +'_' + 'Other observed (design LVK) [yr^{-1}]'

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

        # NAMES.append(namez0_oth)
        # NAMES.append(nameObs_oth)

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
        for ii in range(6):
            datas.append(np.zeros_like(np.zeros(len(MSSFRnameslist))))
            datas.append(np.zeros_like(np.zeros(len(MSSFRnameslist))))


    df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
    df.columns =   df.columns.map(str)
    df.index.names = ['xyz']
    df.columns.names = ['m']

    # print(df) 

    df.to_csv(writePath)








def writeToRatesFile_FormationChannelsRatio(BPSmodelName='Z', DCOtype='BHNS'):



    DCOname = DCOname_dict[DCOtype]

    # path for files 
    path_dir = '/Volumes/Andromeda2/DATA/AllDCO_bugfix/'
    path_ = path_dir
    path_ = path_ + alphabetDirDict[BPSmodelName] +'/'
    path  = path_ + 'COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'




    # read in data 
    fdata = h5.File(path)



    # set optimistic true if that is the variation (H) 
    # OPTIMISTIC=False
    # if (bps_model=='F') or (bps_model=='K'):
    #     OPTIMISTIC=True
    #     print('doing optimistic version of %s'%alphabetDirDict[bps_model])


    seeds = fdata['doubleCompactObjects']['seed'][...].squeeze()
    channels = identify_formation_channels(seeds=seeds, file=fdata)

    
    
    

    stringgg =  'formation_channels_ratio'
    writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_1/Formation_yields_'  + stringgg + '_'  + DCOname + 'fast.csv'     
    headerDict_intrinsic = {5:'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  0:'channel V intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  1:'channel I intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 2:'channel II intrinsic (z=0) [Gpc^{-3} yr^{-1}]',3:'channel III intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 4:'channel IV intrinsic (z=0) [Gpc^{-3} yr^{-1}]'}
    
    headerDict_observed  = {5:'All observed (design LVK) [yr^{-1}]',  0:'channel V observed (design LVK) [yr^{-1}]', 1:'channel I observed (design LVK) [yr^{-1}]', 2:'channel II observed (design LVK) [yr^{-1}]', 3:'channel III observed (design LVK) [yr^{-1}]', 4:'channel IV observed (design LVK) [yr^{-1}]'}    


    
    # get intrinsic weights
    fparam_intrinsic = 'weights_intrinsic'
    # get detected weights
    fparam_detected = 'weights_detected'


    ####################################################
    ######### ITERATE  OVER  MSSFR  MODELS #############
    ####################################################





    df = pd.read_csv(writePath, index_col=0)
    

    for nrC, Channel in enumerate(range(6)):          

    #           #Get the seeds that relate to sorted indices
        mask_C  = (channels==Channel)
        
    
        intrinsicRates = np.zeros(len(MSSFRnameslist))
        detectedRates = np.zeros(len(MSSFRnameslist))       


        for ind_mssfr, mssfr in enumerate(MSSFRnameslist):


            weightheader = 'w_' + mssfr
            w_int = fdata[fparam_intrinsic][weightheader][...].squeeze()
            w_det = fdata[fparam_detected][weightheader][...].squeeze()



            if nrC==5: 
                # TOTAL RATE
                intrinsicRates[ind_mssfr] = np.sum(w_int)
                detectedRates[ind_mssfr]  = np.sum(w_det)  
            else:
                # CHANNEL RATE 
                intrinsicRates[ind_mssfr] = np.sum(w_int[mask_C] / np.sum(w_int))
                detectedRates[ind_mssfr]  = np.sum(w_det[mask_C] / np.sum(w_det))          

            

            
        
        namez0  = BPSmodelName + '_' + headerDict_intrinsic[Channel]
        nameObs = BPSmodelName + '_' + headerDict_observed[Channel]
        df[namez0] = intrinsicRates
        df[nameObs] = detectedRates



    df.to_csv(writePath)


    fdata.close() 

    return



def initalize_formationChannelsRatio(DCOname):

    stringgg =  'formation_channels_ratio'
    writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_1/Formation_yields_'  + stringgg + '_'  + DCOname + 'fast.csv' 

    iii=0


    # CREATE PANDAS FILE 
    nModels=26
    BPSnameslist = list(string.ascii_uppercase)[0:nModels]

    NAMES = []
    # stringgg =  'AllDCOsimulation_formation_channels'
    



    for ind_l, BPSmodelName in enumerate(BPSnameslist):
        # str_z0 = str(L + ' intrinsic (z=0) [Gpc^{-3} yr^{-1}]')
        # str_obs = str(L + ' observed (design LVK) [yr^{-1}]')

        namez0 = BPSmodelName + '_' +'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
        nameObs = BPSmodelName + '_' +'All observed (design LVK) [yr^{-1}]'


        # namez0_oth = BPSmodelName + '_' +'Other intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
        # nameObs_oth = BPSmodelName +'_' + 'Other observed (design LVK) [yr^{-1}]'

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

        # NAMES.append(namez0_oth)
        # NAMES.append(nameObs_oth)

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
        for ii in range(6):
            datas.append(np.zeros_like(np.zeros(len(MSSFRnameslist))))
            datas.append(np.zeros_like(np.zeros(len(MSSFRnameslist))))


    df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
    df.columns =   df.columns.map(str)
    df.index.names = ['xyz']
    df.columns.names = ['m']

    # print(df) 

    df.to_csv(writePath)








#### RUN different simulation summaries : 









INITIALIZE_FormationChannelsRatio = True
INITIALIZE_FormationChannels = False




if INITIALIZE_FormationChannelsRatio==True:
    for DCOtype in ['BHNS','BBH', 'BNS']:
        DCOname = DCOname_dict[DCOtype]
        initalize_formationChannelsRatio(DCOname)


    print('done')


if INITIALIZE_FormationChannels==True:
	for DCOtype in ['BHNS','BBH', 'BNS']:
		DCOname = DCOname_dict[DCOtype]
		initalize_formationChannels(DCOname)


	print('done')



runFormationChannels=False
runFormationChannelsRatio=True



if runFormationChannelsRatio ==True:
    for BPS in  BPSnameslist[:]:
        print(BPS)
        for DCOtype in ['BHNS','BBH', 'BNS']:
            print('at DCOtype =', DCOtype)
            writeToRatesFile_FormationChannelsRatio(BPSmodelName=BPS, DCOtype=DCOtype)
            print('done with ', BPS)






if runFormationChannels ==True:
    for BPS in  BPSnameslist[:]:
        print(BPS)
        for DCOtype in ['BHNS','BBH', 'BNS']:
            print('at DCOtype =', DCOtype)
            writeToRatesFile_FormationChannels(BPSmodelName=BPS, DCOtype=DCOtype)
            print('done with ', BPS)



