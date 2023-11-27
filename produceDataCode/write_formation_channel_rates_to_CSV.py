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

    headerDict_intrinsic = { 5:'All intrinsic (z=0) [Gpc^-3 yr^-1]',  0:'channel V intrinsic (z=0)',  1:'channel I intrinsic (z=0)', 2:'channel II intrinsic (z=0)',3:'channel III intrinsic (z=0)', 4:'channel IV intrinsic (z=0)'}
    stringgg =  'Formation_Channels_Local_Rates'

    for MSSFR in MSSFRnameslist:
        writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_1/'  + stringgg + '_'  + DCOname + '_xyz_' + MSSFR + '.csv' 

        NAMES = []
        for Channel in range(len(headerDict_intrinsic)):
            namez0 = headerDict_intrinsic[Channel]      
            NAMES.append(namez0)

        datas=[]
        for ii in range(len(headerDict_intrinsic)):
            datas.append(np.zeros_like(BPSnameslist))    
        
        df = pd.DataFrame(data=datas, index=NAMES, columns=BPSnameslist).T
        df.columns =   df.columns.map(str)
        df.index.names = ['formation channel']
        df.columns.names = ['model']
        df.to_csv(writePath)

    return 



def writeToRatesFile_FormationChannels(pathData='/Volumes/SimonsFoundation/DataDCO/', DCOtype='BHNS'):
    """ DCOType = ['BBH', 'BHNS', BNS]"""

    headerDict_intrinsic = { 5:'All intrinsic (z=0) [Gpc^-3 yr^-1]',  0:'channel V intrinsic (z=0)',  1:'channel I intrinsic (z=0)', 2:'channel II intrinsic (z=0)',3:'channel III intrinsic (z=0)', 4:'channel IV intrinsic (z=0)'}
    DCOname = DCOname_dict[DCOtype]
    which_z_ind=0 # gives lowest redshift. 
    stringgg =  'Formation_Channels_Local_Rates'
    fparam_key = "formationchannel_z_rates"


    for ind_L, MSSFRname in enumerate(MSSFRnameslist[:]):
        print('at MSSFR: ', MSSFRname)
        #iterate over the formation channels 
        for ind_c, whichChannel in enumerate(['classic', 'stable B no CEE',  'immediate CE',  r'double-core CE', 'other']): #'vii', , 'vi'
            merger_ratio_z = np.zeros(len(BPSnameslist))
            index_channel = dictFormationChannelIndex[whichChannel]
            #iterate over stellar evolution models 
            for ind_m, BPSmodelName in enumerate(BPSnameslist):
                full_data_path  = pathData + alphabetDirDict[BPSmodelName] +'/COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'
                fparam_key = "formationchannel_z_rates"
                header = "fraction_" + whichChannel + "_" + MSSFRname
                fdata = h5.File(full_data_path,'r')
                fc_fraction = fdata[fparam_key][header][...].squeeze()[which_z_ind]  # [which_z_ind = 0] gives at lowest redshift
                
                if whichChannel=='classic':
                    header = "fraction_" + 'vi' + "_" + MSSFRname
                    fc_fraction += fdata[fparam_key][header][...].squeeze()[which_z_ind]  # [which_z_ind = 0] gives at lowest redshift  
                elif whichChannel=='stable B no CEE':
                    header = "fraction_" + 'vii' + "_" + MSSFRname
                    fc_fraction += fdata[fparam_key][header][...].squeeze()[which_z_ind]  # [which_z_ind = 0] gives at lowest redshift  
                merger_ratio_z[ind_m] = fc_fraction
                fdata.close()
                
            writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_1/'  + stringgg + '_'  + DCOname + '_xyz_' + MSSFRname + '.csv' 
            df = pd.read_csv(writePath, index_col=0)                
            column_name = headerDict_intrinsic[index_channel] 
            df[column_name] = merger_ratio_z
            df.to_csv(writePath)

        # add total rate: 
        totals_z = np.zeros(len(BPSnameslist))
        for ind_m, BPSmodelName in enumerate(BPSnameslist):
            full_data_path  = pathData + alphabetDirDict[BPSmodelName] +'/COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'
            header = "rate_" + "total" + "_" + MSSFRname
            fdata = h5.File(full_data_path,'r')
            totals_z[ind_m] =  fdata[fparam_key][header][...].squeeze()[which_z_ind] # [which_z_ind = 0] gives at lowest redshift
            fdata.close()

        writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_1/'  + stringgg + '_'  + DCOname + '_xyz_' + MSSFRname + '.csv' 
        df = pd.read_csv(writePath, index_col=0)                
        column_name = 'All intrinsic (z=0) [Gpc^-3 yr^-1]' 
        df[column_name] = totals_z
        df.to_csv(writePath)


    return



INITIALIZE_FormationChannels = True #True

if INITIALIZE_FormationChannels==True:
    for DCOtype in ['BHNS','BBH', 'BNS']:
        DCOname = DCOname_dict[DCOtype]
        initalize_formationChannels(DCOname)
    print('done')

runFormationChannels=True# #False

if runFormationChannels ==True:
    for DCOtype in ['BHNS', 'BNS', 'BBH']:
        print('at DCOtype =', DCOtype)
        writeToRatesFile_FormationChannels(DCOtype=DCOtype)
        print('done with ', DCOtype)






# def initalize_formationChannels(DCOname):


#     stringgg =  'Formation_Channels_Local_Rates'
#     writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_1/'  + stringgg + '_'  + DCOname + '.csv' 

#     iii=0


#     # CREATE PANDAS FILE 
#     nModels=26
#     BPSnameslist = list(string.ascii_uppercase)[0:nModels]

#     NAMES = []

#     # headerDict_intrinsic = { 5:'Channel VI intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  6:'Channel VII intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 7:'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  0:'channel V intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  1:'channel I intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 2:'channel II intrinsic (z=0) [Gpc^{-3} yr^{-1}]',3:'channel III intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 4:'channel IV intrinsic (z=0) [Gpc^{-3} yr^{-1}]'}
#     # headerDict_observed  = { 5:'channel VI observed (design LVK) [yr^{-1}]',     6:'channel VII observed (design LVK) [yr^{-1}]',    7:'All observed (design LVK) [yr^{-1}]',     0:'channel V observed (design LVK) [yr^{-1}]', 1:'channel I observed (design LVK) [yr^{-1}]', 2:'channel II observed (design LVK) [yr^{-1}]', 3:'channel III observed (design LVK) [yr^{-1}]', 4:'channel IV observed (design LVK) [yr^{-1}]'}    
#     # enumerate_list = range(8)

#     headerDict_intrinsic = { 5:'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  0:'channel V intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  1:'channel I intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 2:'channel II intrinsic (z=0) [Gpc^{-3} yr^{-1}]',3:'channel III intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 4:'channel IV intrinsic (z=0) [Gpc^{-3} yr^{-1}]'}
#     # headerDict_observed  = { 5:'All observed (design LVK) [yr^{-1}]',     0:'channel V observed (design LVK) [yr^{-1}]', 1:'channel I observed (design LVK) [yr^{-1}]', 2:'channel II observed (design LVK) [yr^{-1}]', 3:'channel III observed (design LVK) [yr^{-1}]', 4:'channel IV observed (design LVK) [yr^{-1}]'}    
    
#     enumerate_list = range(len(headerDict_intrinsic))
    
#     for ind_c, Channel in enumerate(enumerate_list):
#         for ind_l, BPSmodelName in enumerate(BPSnameslist):
#             namez0 = BPSmodelName + '_' + headerDict_intrinsic[Channel]      
#             NAMES.append(namez0)

#     datas=[]
#     for ii in range(len(headerDict_intrinsic)):
#         for i in range(len(BPSnameslist)):
#             datas.append(np.zeros_like(MSSFRnameslist))    
        
#     df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
#     df.columns =   df.columns.map(str)
#     df.index.names = ['xyz']
#     df.columns.names = ['m']


#     df.to_csv(writePath)





# def writeToRatesFile_FormationChannels(BPSmodelName='Z', DCOtype='BHNS'):



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

#     immediateRLOFAfterCEE = fdata["commonEnvelopes"]["immediateRLOFAfterCEE"][...].squeeze()
#     immediateRLOFAfterCEEmask = (immediateRLOFAfterCEE==1)
#     CErandomSeed = fdata["commonEnvelopes"]["randomSeed"][...].squeeze()
    
#     seeds_withImmediateRLOFAfterCEE = np.in1d(fdata['doubleCompactObjects']['seed'][...].squeeze(), CErandomSeed[immediateRLOFAfterCEEmask])
#     print('systems with withImmediateRLOFAfterCEE',np.sum(seeds_withImmediateRLOFAfterCEE), ' which is a fraction', np.sum(seeds_withImmediateRLOFAfterCEE)/len(seeds_withImmediateRLOFAfterCEE))

#     seeds = fdata['doubleCompactObjects']['seed'][...].squeeze()
#     channels = identify_formation_channels(seeds=seeds, file=fdata)

    
    
    

#     stringgg =  'formation_channels'
#     writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_1/Formation_yields_'  + stringgg + '_'  + DCOname + 'fast.csv'     

#     headerDict_intrinsic = { 5:'Channel VI intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  6:'Channel VII intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 7:'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  0:'channel V intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  1:'channel I intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 2:'channel II intrinsic (z=0) [Gpc^{-3} yr^{-1}]',3:'channel III intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 4:'channel IV intrinsic (z=0) [Gpc^{-3} yr^{-1}]'}
#     headerDict_observed  = { 5:'channel VI observed (design LVK) [yr^{-1}]',     6:'channel VII observed (design LVK) [yr^{-1}]',    7:'All observed (design LVK) [yr^{-1}]',     0:'channel V observed (design LVK) [yr^{-1}]', 1:'channel I observed (design LVK) [yr^{-1}]', 2:'channel II observed (design LVK) [yr^{-1}]', 3:'channel III observed (design LVK) [yr^{-1}]', 4:'channel IV observed (design LVK) [yr^{-1}]'}    
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

#     #           #Get the seeds that relate to sorted indices
#         mask_C  = (channels==Channel)[~seeds_withImmediateRLOFAfterCEE]# everything except these seeds  
        
    
#         intrinsicRates = np.zeros(len(MSSFRnameslist))
#         detectedRates = np.zeros(len(MSSFRnameslist))       


#         for ind_mssfr, mssfr in enumerate(MSSFRnameslist):


#             weightheader = 'w_' + mssfr
#             w_int = fdata[fparam_intrinsic][weightheader][...].squeeze()[~seeds_withImmediateRLOFAfterCEE]# everything except these seeds
#             w_det = fdata[fparam_detected][weightheader][...].squeeze()[~seeds_withImmediateRLOFAfterCEE]# everything except these seeds



#             if nrC==enumerate_list[-1]: 
#                 # TOTAL RATE
#                 intrinsicRates[ind_mssfr] = np.sum(w_int)
#                 detectedRates[ind_mssfr]  = np.sum(w_det)  
#             else:
#                 # CHANNEL RATE 
#                 intrinsicRates[ind_mssfr] = np.sum(w_int[mask_C])
#                 detectedRates[ind_mssfr]  = np.sum(w_det[mask_C])          

            

            
        
#         namez0  = BPSmodelName + '_' + headerDict_intrinsic[Channel]
#         nameObs = BPSmodelName + '_' + headerDict_observed[Channel]
#         df[namez0] = intrinsicRates
#         df[nameObs] = detectedRates



#     df.to_csv(writePath)


#     fdata.close() 

#     return









# def writeToRatesFile_FormationChannelsRatio(BPSmodelName='Z', DCOtype='BHNS'):



#     DCOname = DCOname_dict[DCOtype]

#     # path for files 
#     path_dir = '/Volumes/Andromeda2/DATA/AllDCO_bugfix/'
#     path_ = path_dir
#     path_ = path_ + alphabetDirDict[BPSmodelName] +'/'
#     path  = path_ + 'COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'



#     headerDict_intrinsic = { 5:'Channel VI intrinsic (z=0) ratio',  6:'Channel VII intrinsic (z=0) ratio', 7:'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  0:'channel V intrinsic (z=0) ratio',  1:'channel I intrinsic (z=0) ratio', 2:'channel II intrinsic (z=0) ratio',3:'channel III intrinsic (z=0) ratio', 4:'channel IV intrinsic (z=0) ratio'}
#     enumerate_list = range(8)
#     headerDict_Z  = { 5:'channelVI',     6:'channelVII',    7:'All',     0:'channelV', 1:'channelI', 2:'channelII', 3:'channelIII', 4:'channelIV'} 
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
    

    
#     # get intrinsic weights
#     fparam_intrinsic = 'weights_intrinsic'
#     # get detected weights
#     fparam_detected = 'weights_detected'


#     ####################################################
#     ######### ITERATE  OVER  MSSFR  MODELS #############
#     ####################################################





    
    

#     for nrC, Channel in enumerate(enumerate_list):   


        
#         writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_1/Formation_yields_'  + stringgg + '_'  + DCOname + '_'  + headerDict_Z[Channel] + '.csv' 
#         print(writePath)
#         df = pd.read_csv(writePath, index_col=0)

#         print('now at Channel ', Channel)     

#     #           #Get the seeds that relate to sorted indices
#         mask_C  = (channels==Channel)
        
    
#         intrinsicRates = np.zeros(len(MSSFRnameslist))
#         intrinsicRates_all = np.zeros(len(MSSFRnameslist))       


#         for ind_mssfr, mssfr in enumerate(MSSFRnameslist):


#             weightheader = 'w_' + mssfr
#             w_int = fdata[fparam_intrinsic][weightheader][...].squeeze()
#             # w_det = fdata[fparam_detected][weightheader][...].squeeze()



#             if nrC==enumerate_list[-1]: 
#                 # TOTAL RATE
#                 intrinsicRates[ind_mssfr] = np.sum(w_int)
#                 intrinsicRates_all[ind_mssfr]  = np.sum(w_int)  
#             else:
#                 # CHANNEL RATE 
#                 intrinsicRates[ind_mssfr] = np.sum(w_int[mask_C])/np.sum(w_int)
#                 intrinsicRates_all[ind_mssfr]  = np.sum(w_int)         

            

            
        
#         namez0  = BPSmodelName + '_' + headerDict_intrinsic[Channel]
#         nameAll = BPSmodelName + '_' + headerDict_intrinsic[7]   
#         df[namez0] = intrinsicRates
#         df[nameAll] = intrinsicRates_all



#         df.to_csv(writePath)


#     fdata.close() 

#     return




# def initalize_formationChannelsRatio(DCOname):






#     headerDict_intrinsic = { 5:'Channel VI intrinsic (z=0) ratio',  6:'Channel VII intrinsic (z=0) ratio', 7:'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  0:'channel V intrinsic (z=0) ratio',  1:'channel I intrinsic (z=0) ratio', 2:'channel II intrinsic (z=0) ratio',3:'channel III intrinsic (z=0) ratio', 4:'channel IV intrinsic (z=0) ratio'}
#     enumerate_list = range(8)
#     headerDict_Z  = { 5:'channelVI',     6:'channelVII',    7:'All',     0:'channelV', 1:'channelI', 2:'channelII', 3:'channelIII', 4:'channelIV'} 

#     for ind_c, Channel in enumerate(enumerate_list):


#         stringgg =  'formation_channels_ratio'
#         writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_1/Formation_yields_'  + stringgg + '_'  + DCOname + '_'  + headerDict_Z[ind_c] + '.csv' 

#         iii=0


#         # CREATE PANDAS FILE 
#         nModels=26
#         BPSnameslist = list(string.ascii_uppercase)[0:nModels]

#         NAMES = []
#         stringgg =  'AllDCOsimulation_formation_channels'



#         # headerDict_intrinsic = { 5:'Channel VI intrinsic (z=0) ratio',  6:'Channel VII intrinsic (z=0) ratio', 7:'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]',  0:'channel V intrinsic (z=0) ratio',  1:'channel I intrinsic (z=0) ratio', 2:'channel II intrinsic (z=0) ratio',3:'channel III intrinsic (z=0) ratio', 4:'channel IV intrinsic (z=0) ratio'}
#         # headerDict_observed  = { 5:'channel VI observed (design LVK) ratio',     6:'channel VII observed (design LVK) ratio',    7:'All observed (design LVK) [yr^{-1}]',     0:'channel V observed (design LVK) ratio', 1:'channel I observed (design LVK) ratio', 2:'channel II observed (design LVK) ratio', 3:'channel III observed (design LVK) ratio', 4:'channel IV observed (design LVK) ratio'}    
#         # enumerate_list = range(8)
        
#         for ind_l, BPSmodelName in enumerate(BPSnameslist):
#             # for ind_c, Channel in enumerate(enumerate_list):
#             namez0 = BPSmodelName + '_' + headerDict_intrinsic[Channel]
#             nameAll = BPSmodelName + '_' + headerDict_intrinsic[7]           
#             NAMES.append(namez0)
#             NAMES.append(nameAll)

#         datas=[]

#         for i in range(len(BPSnameslist)):
#         # for ii in range(8):
#             datas.append(np.zeros_like(MSSFRnameslist))
#             datas.append(np.zeros_like(MSSFRnameslist))

            
#         df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
#         df.columns =   df.columns.map(str)
#         df.index.names = ['xyz']
#         df.columns.names = ['m']


#         df.to_csv(writePath)









#### RUN different simulation summaries : 









# INITIALIZE_FormationChannelsRatio = False #True
# INITIALIZE_FormationChannels = True #True




# if INITIALIZE_FormationChannelsRatio==True:
#     for DCOtype in ['BHNS','BBH', 'BNS']:
#         DCOname = DCOname_dict[DCOtype]
#         initalize_formationChannelsRatio(DCOname)


#     print('done')


# if INITIALIZE_FormationChannels==True:
# 	for DCOtype in ['BHNS','BBH', 'BNS']:
# 		DCOname = DCOname_dict[DCOtype]
# 		initalize_formationChannels(DCOname)


# 	print('done')



# runFormationChannels=False# #False
# runFormationChannelsRatio= False #True





# if runFormationChannelsRatio ==True:
#     for BPS in  BPSnameslist[:]:
#         print(BPS)
#         for DCOtype in ['BBH', 'BHNS','BNS']:
#             print('at DCOtype =', DCOtype)
#             writeToRatesFile_FormationChannelsRatio(BPSmodelName=BPS, DCOtype=DCOtype)
#             print('done with ', BPS)

# if runFormationChannels ==True:
#     for BPS in  BPSnameslist[:]:
#         print(BPS)
#         for DCOtype in ['BHNS','BBH', 'BNS']:
#             print('at DCOtype =', DCOtype)
#             writeToRatesFile_FormationChannels(BPSmodelName=BPS, DCOtype=DCOtype)
#             print('done with ', BPS)



