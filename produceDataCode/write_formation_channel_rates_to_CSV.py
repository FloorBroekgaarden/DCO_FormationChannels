import numpy as np
import sys
import h5py as h5

#Quick fudge to make import from ../Scripts work
sys.path.append('../../common_code/')


from PostProcessingScripts import * 
from formation_channels import * 

import pandas as pd
from astropy import units as u
from astropy import constants as const





def initalize_formationChannels(DCOname):

    stringgg =  'formation_channels'
    writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_1/Formation_yields_'  + stringgg + '_'  + DCOname + 'fast.csv' 

    iii=0


    # CREATE PANDAS FILE 
    nModels=26
    BPSnameslist = list(string.ascii_uppercase)[0:nModels]

    NAMES = []
    stringgg =  'AllDCOsimulation_formation_channels'



    headerDict_intrinsic = { 5:'Channel VI intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  6:'Channel VII intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 7:'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  0:'channel V intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  1:'channel I intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 2:'channel II intrinsic (z=0) [Gpc^{-3} yr^{-1}]',3:'channel III intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 4:'channel IV intrinsic (z=0) [Gpc^{-3} yr^{-1}]'}
    headerDict_observed  = { 5:'channel VI observed (design LVK) [yr^{-1}]',     6:'channel VII observed (design LVK) [yr^{-1}]',    7:'All observed (design LVK) [yr^{-1}]',     0:'channel V observed (design LVK) [yr^{-1}]', 1:'channel I observed (design LVK) [yr^{-1}]', 2:'channel II observed (design LVK) [yr^{-1}]', 3:'channel III observed (design LVK) [yr^{-1}]', 4:'channel IV observed (design LVK) [yr^{-1}]'}    
    enumerate_list = range(8)
    
    for ind_l, BPSmodelName in enumerate(BPSnameslist):
        for ind_c, Channel in enumerate(enumerate_list):
            namez0 = BPSmodelName + '_' + headerDict_intrinsic[Channel]
            nameObs = BPSmodelName + '_' + headerDict_observed[Channel]            
            NAMES.append(namez0)
            NAMES.append(nameObs)

    datas=[]

    for i in range(len(BPSnameslist)):
        for ii in range(8):
            datas.append(np.zeros_like(MSSFRnameslist))
            datas.append(np.zeros_like(MSSFRnameslist))
        
        
    df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
    df.columns =   df.columns.map(str)
    df.index.names = ['xyz']
    df.columns.names = ['m']


    df.to_csv(writePath)





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

    headerDict_intrinsic = { 5:'Channel VI intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  6:'Channel VII intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 7:'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  0:'channel V intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  1:'channel I intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 2:'channel II intrinsic (z=0) [Gpc^{-3} yr^{-1}]',3:'channel III intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 4:'channel IV intrinsic (z=0) [Gpc^{-3} yr^{-1}]'}
    headerDict_observed  = { 5:'channel VI observed (design LVK) [yr^{-1}]',     6:'channel VII observed (design LVK) [yr^{-1}]',    7:'All observed (design LVK) [yr^{-1}]',     0:'channel V observed (design LVK) [yr^{-1}]', 1:'channel I observed (design LVK) [yr^{-1}]', 2:'channel II observed (design LVK) [yr^{-1}]', 3:'channel III observed (design LVK) [yr^{-1}]', 4:'channel IV observed (design LVK) [yr^{-1}]'}    
    enumerate_list = range(8)

    
    # get intrinsic weights
    fparam_intrinsic = 'weights_intrinsic'
    # get detected weights
    fparam_detected = 'weights_detected'


    ####################################################
    ######### ITERATE  OVER  MSSFR  MODELS #############
    ####################################################





    df = pd.read_csv(writePath, index_col=0)
    

    for nrC, Channel in enumerate(enumerate_list):          

    #           #Get the seeds that relate to sorted indices
        mask_C  = (channels==Channel)
        
    
        intrinsicRates = np.zeros(len(MSSFRnameslist))
        detectedRates = np.zeros(len(MSSFRnameslist))       


        for ind_mssfr, mssfr in enumerate(MSSFRnameslist):


            weightheader = 'w_' + mssfr
            w_int = fdata[fparam_intrinsic][weightheader][...].squeeze()
            w_det = fdata[fparam_detected][weightheader][...].squeeze()



            if nrC==enumerate_list[-1]: 
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











 ##### RATIO #######



# def writeToRatesFile_FormationChannelsRatio(BPSmodelName='Z', DCOtype='BHNS'):



#     DCOname = DCOname_dict[DCOtype]

#     # path for files 
#     path_dir = '/Volumes/Andromeda2/DATA/AllDCO_bugfix/'
#     path_ = path_dir
#     path_ = path_ + alphabetDirDict[BPSmodelName] +'/'
#     path  = path_ + 'COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'




#     # read in data 
#     fdata = h5.File(path)



#     # set optimistic true if that is the variation (H) 
#     # OPTIMISTIC=False
#     # if (bps_model=='F') or (bps_model=='K'):
#     #     OPTIMISTIC=True
#     #     print('doing optimistic version of %s'%alphabetDirDict[bps_model])


#     seeds = fdata['doubleCompactObjects']['seed'][...].squeeze()
#     channels = identify_formation_channels(seeds=seeds, file=fdata)

    
    
    

#     stringgg =  'formation_channels_ratio'
#     writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_1/Formation_yields_'  + stringgg + '_'  + DCOname + 'fast.csv'     

#     headerDict_intrinsic = { 5:'Channel VI intrinsic (z=0) ratio',  6:'Channel VII intrinsic (z=0) ratio', 7:'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  0:'channel V intrinsic (z=0) ratio',  1:'channel I intrinsic (z=0) ratio', 2:'channel II intrinsic (z=0) ratio',3:'channel III intrinsic (z=0) ratio', 4:'channel IV intrinsic (z=0) ratio'}
#     headerDict_observed  = { 5:'channel VI observed (design LVK) ratio',     6:'channel VII observed (design LVK) ratio',    7:'All observed (design LVK) [yr^{-1}]',     0:'channel V observed (design LVK) ratio', 1:'channel I observed (design LVK) ratio', 2:'channel II observed (design LVK) ratio', 3:'channel III observed (design LVK) ratio', 4:'channel IV observed (design LVK) ratio'}    
#     enumerate_list = range(8)

    
#     # get intrinsic weights
#     fparam_intrinsic = 'weights_intrinsic'
#     # get detected weights
#     fparam_detected = 'weights_detected'


#     ####################################################
#     ######### ITERATE  OVER  MSSFR  MODELS #############
#     ####################################################





#     df = pd.read_csv(writePath, index_col=0)
    

#     for nrC, Channel in enumerate(enumerate_list):     
#         print('now at Channel ', Channel)     

#     #           #Get the seeds that relate to sorted indices
#         mask_C  = (channels==Channel)
        
    
#         intrinsicRates = np.zeros(len(MSSFRnameslist))
#         detectedRates = np.zeros(len(MSSFRnameslist))       


#         for ind_mssfr, mssfr in enumerate(MSSFRnameslist):


#             weightheader = 'w_' + mssfr
#             w_int = fdata[fparam_intrinsic][weightheader][...].squeeze()
#             w_det = fdata[fparam_detected][weightheader][...].squeeze()



#             if nrC==enumerate_list[-1]: 
#                 # TOTAL RATE
#                 intrinsicRates[ind_mssfr] = np.sum(w_int)
#                 detectedRates[ind_mssfr]  = np.sum(w_det)  
#             else:
#                 # CHANNEL RATE 
#                 intrinsicRates[ind_mssfr] = np.sum(w_int[mask_C])/np.sum(w_int)
#                 detectedRates[ind_mssfr]  = np.sum(w_det[mask_C])/np.sum(w_det)          

            

            
        
#         namez0  = BPSmodelName + '_' + headerDict_intrinsic[Channel]
#         nameObs = BPSmodelName + '_' + headerDict_observed[Channel]
#         df[namez0] = intrinsicRates
#         df[nameObs] = detectedRates



#     df.to_csv(writePath)


#     fdata.close() 

#     return




def writeToRatesFile_FormationChannelsRatio(BPSmodelName='Z', DCOtype='BHNS'):



    DCOname = DCOname_dict[DCOtype]

    # path for files 
    path_dir = '/Volumes/Andromeda2/DATA/AllDCO_bugfix/'
    path_ = path_dir
    path_ = path_ + alphabetDirDict[BPSmodelName] +'/'
    path  = path_ + 'COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'



    headerDict_intrinsic = { 5:'Channel VI intrinsic (z=0) ratio',  6:'Channel VII intrinsic (z=0) ratio', 7:'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  0:'channel V intrinsic (z=0) ratio',  1:'channel I intrinsic (z=0) ratio', 2:'channel II intrinsic (z=0) ratio',3:'channel III intrinsic (z=0) ratio', 4:'channel IV intrinsic (z=0) ratio'}
    enumerate_list = range(8)
    headerDict_Z  = { 5:'channelVI',     6:'channelVII',    7:'All',     0:'channelV', 1:'channelI', 2:'channelII', 3:'channelIII', 4:'channelIV'} 
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
    

    
    # get intrinsic weights
    fparam_intrinsic = 'weights_intrinsic'
    # get detected weights
    fparam_detected = 'weights_detected'


    ####################################################
    ######### ITERATE  OVER  MSSFR  MODELS #############
    ####################################################





    
    

    for nrC, Channel in enumerate(enumerate_list):   


        
        writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_1/Formation_yields_'  + stringgg + '_'  + DCOname + '_'  + headerDict_Z[Channel] + '.csv' 
        print(writePath)
        df = pd.read_csv(writePath, index_col=0)

        print('now at Channel ', Channel)     

    #           #Get the seeds that relate to sorted indices
        mask_C  = (channels==Channel)
        
    
        intrinsicRates = np.zeros(len(MSSFRnameslist))
        intrinsicRates_all = np.zeros(len(MSSFRnameslist))       


        for ind_mssfr, mssfr in enumerate(MSSFRnameslist):


            weightheader = 'w_' + mssfr
            w_int = fdata[fparam_intrinsic][weightheader][...].squeeze()
            # w_det = fdata[fparam_detected][weightheader][...].squeeze()



            if nrC==enumerate_list[-1]: 
                # TOTAL RATE
                intrinsicRates[ind_mssfr] = np.sum(w_int)
                intrinsicRates_all[ind_mssfr]  = np.sum(w_int)  
            else:
                # CHANNEL RATE 
                intrinsicRates[ind_mssfr] = np.sum(w_int[mask_C])/np.sum(w_int)
                intrinsicRates_all[ind_mssfr]  = np.sum(w_int)         

            

            
        
        namez0  = BPSmodelName + '_' + headerDict_intrinsic[Channel]
        nameAll = BPSmodelName + '_' + headerDict_intrinsic[7]   
        df[namez0] = intrinsicRates
        df[nameAll] = intrinsicRates_all



        df.to_csv(writePath)


    fdata.close() 

    return




def initalize_formationChannelsRatio(DCOname):






    headerDict_intrinsic = { 5:'Channel VI intrinsic (z=0) ratio',  6:'Channel VII intrinsic (z=0) ratio', 7:'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  0:'channel V intrinsic (z=0) ratio',  1:'channel I intrinsic (z=0) ratio', 2:'channel II intrinsic (z=0) ratio',3:'channel III intrinsic (z=0) ratio', 4:'channel IV intrinsic (z=0) ratio'}
    enumerate_list = range(8)
    headerDict_Z  = { 5:'channelVI',     6:'channelVII',    7:'All',     0:'channelV', 1:'channelI', 2:'channelII', 3:'channelIII', 4:'channelIV'} 

    for ind_c, Channel in enumerate(enumerate_list):


        stringgg =  'formation_channels_ratio'
        writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_1/Formation_yields_'  + stringgg + '_'  + DCOname + '_'  + headerDict_Z[ind_c] + '.csv' 

        iii=0


        # CREATE PANDAS FILE 
        nModels=26
        BPSnameslist = list(string.ascii_uppercase)[0:nModels]

        NAMES = []
        stringgg =  'AllDCOsimulation_formation_channels'



        # headerDict_intrinsic = { 5:'Channel VI intrinsic (z=0) ratio',  6:'Channel VII intrinsic (z=0) ratio', 7:'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  0:'channel V intrinsic (z=0) ratio',  1:'channel I intrinsic (z=0) ratio', 2:'channel II intrinsic (z=0) ratio',3:'channel III intrinsic (z=0) ratio', 4:'channel IV intrinsic (z=0) ratio'}
        # headerDict_observed  = { 5:'channel VI observed (design LVK) ratio',     6:'channel VII observed (design LVK) ratio',    7:'All observed (design LVK) [yr^{-1}]',     0:'channel V observed (design LVK) ratio', 1:'channel I observed (design LVK) ratio', 2:'channel II observed (design LVK) ratio', 3:'channel III observed (design LVK) ratio', 4:'channel IV observed (design LVK) ratio'}    
        # enumerate_list = range(8)
        
        for ind_l, BPSmodelName in enumerate(BPSnameslist):
            # for ind_c, Channel in enumerate(enumerate_list):
            namez0 = BPSmodelName + '_' + headerDict_intrinsic[Channel]
            nameAll = BPSmodelName + '_' + headerDict_intrinsic[7]           
            NAMES.append(namez0)
            NAMES.append(nameAll)

        datas=[]

        for i in range(len(BPSnameslist)):
        # for ii in range(8):
            datas.append(np.zeros_like(MSSFRnameslist))
            datas.append(np.zeros_like(MSSFRnameslist))

            
        df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
        df.columns =   df.columns.map(str)
        df.index.names = ['xyz']
        df.columns.names = ['m']


        df.to_csv(writePath)









#### RUN different simulation summaries : 









INITIALIZE_FormationChannelsRatio = True #False #True
INITIALIZE_FormationChannels = False #True




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
runFormationChannelsRatio=True# False #True





if runFormationChannelsRatio ==True:
    for BPS in  BPSnameslist[:]:
        print(BPS)
        for DCOtype in ['BBH', 'BHNS','BNS']:
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



