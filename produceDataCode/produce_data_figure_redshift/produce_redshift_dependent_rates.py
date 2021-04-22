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


MSSFRnameslistWantedOrder = []
MSSFRnameslistWantedOrder.append('000') # add phenomenological 

for ind_GSMF, GSMF in enumerate([0,1,2]):
    ind_x = ind_GSMF + 1
    for ind_MZ, MZ in enumerate([0,1,2]):
        ind_y = ind_MZ +1
        for ind_SFR, SFR in enumerate([0,1,2]):
            ind_z = ind_SFR+1
            
            
            
            
        

            MSSFRnameslistWantedOrder.append('%s%s%s'%(ind_x, ind_y, ind_z))




def writeToRatesFile_FormationChannels(BPSmodelName='Z', DCOtype='BHNS'):



    DCOname = DCOname_dict[DCOtype]

    # path for files 
    path_dir = '/Volumes/Andromeda/DATA/AllDCO_bugfix/'
    path_ = path_dir
    path_ = path_ + alphabetDirDict[BPSmodelName] +'/'
    path  = path_ + 'COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'




    # read in data 
    fdata = h5.File(path)



    # set optimistic true if that is the variation (H) 
    OPTIMISTIC=False
    if BPSmodelName=='H':
        OPTIMISTIC=True 


    seeds = fdata['doubleCompactObjects']['seed'][...].squeeze()
    channels = identify_formation_channels(seeds=seeds, file=fdata)

    
    
    

    stringgg =  'redshift'
    writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_redshift/Formation_yields_'  + stringgg + '_'  + DCOname + 'fast.csv'     
    headerDict_intrinsic = {0:'channel I intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 1:'channel II intrinsic (z=0) [Gpc^{-3} yr^{-1}]',2:'channel III intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 3:'channel IV intrinsic (z=0) [Gpc^{-3} yr^{-1}]', 4:'channel V intrinsic (z=0) [Gpc^{-3} yr^{-1}]'}
    
    headerDict_observed  = {1:'channel I observed (design LVK) [yr^{-1}]', 2:'channel II observed (design LVK) [yr^{-1}]', 3:'channel III observed (design LVK) [yr^{-1}]', 4:'channel IV observed (design LVK) [yr^{-1}]', 0:'channel V observed (design LVK) [yr^{-1}]'}    

    
    # get intrinsic weights
    fparam_intrinsic = 'weights_intrinsic'
    # get detected weights
    fparam_detected = 'weights_detected'


    ####################################################
    ######### ITERATE  OVER  MSSFR  MODELS #############
    ####################################################





    df = pd.read_csv(writePath, index_col=0)
    

    for nrC, Channel in enumerate(list(np.unique(channels))):              
    #           #Get the seeds that relate to sorted indices
        mask_C  = (channels==Channel)
        
        intrinsicRates = np.zeros(len(MSSFRnameslist))
        detectedRates = np.zeros(len(MSSFRnameslist))
        intrinsicRates_all = np.zeros(len(MSSFRnameslist))
        detectedRates_all  = np.zeros(len(MSSFRnameslist))        


        for ind_mssfr, mssfr in enumerate(MSSFRnameslist):

            weightheader = 'w_' + mssfr
            w_int = fdata[fparam_intrinsic][weightheader][...].squeeze()
            w_det = fdata[fparam_detected][weightheader][...].squeeze()

            # CHANNEL RATE 
            intrinsicRates[ind_mssfr] = np.sum(w_int[mask_C])
            detectedRates[ind_mssfr]  = np.sum(w_det[mask_C])

            if nrC==0: 
                # TOTAL RATE
                intrinsicRates_all[ind_mssfr] = np.sum(w_int)
                detectedRates_all[ind_mssfr]  = np.sum(w_det)            

            
        if nrC==0: 
            namez0_all  = BPSmodelName + '_' + 'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
            nameObs_all = BPSmodelName + '_' + 'All observed (design LVK) [yr^{-1}]'
            df[namez0_all] = intrinsicRates_all
            df[nameObs_all] = detectedRates_all
            
        
        namez0  = BPSmodelName + '_' + headerDict_intrinsic[Channel]
        nameObs = BPSmodelName + '_' + headerDict_observed[Channel]

        df[namez0] = intrinsicRates
        df[nameObs] = detectedRates





    df.to_csv(writePath)


    fdata.close() 

    return



def initalize_formationChannels(DCOname):


    iii=0


    # CREATE PANDAS FILE 
    nModels=17
    BPSnameslist = list(string.ascii_uppercase)[0:nModels]


    minz = 0.
    if DCOtype=='BHNS':
        redshifts = np.arange(0.001, 0.499, step=0.002)
    elif DCOtype=='BNS':
        redshifts = np.arange(0.0012, 0.2488, step=0.0038-0.0012)
    elif DCOtype=='BBH': 
        redshifts =  np.arange(0.0125, 2.4875, step=0.0375-0.0125)




    for ind_l, BPSmodelName in enumerate(BPSnameslist):






        NAMES = []
        stringgg =  'redshift'
        writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_redshift/Formation_yields_'  + stringgg + '_'  + DCOname + '_' +  BPSmodelName+ '.csv' 

        NAMES.append('redshift')

        datas=[]
        datas.append(redshifts)

        for ind_mssfr, mssfr_name in enumerate(MSSFRnameslistWantedOrder):

            namez0     =  "R_%s"%(mssfr_name)+ '_' +'All intrinsic [Gpc^{-3} yr^{-1}]'
            namez0_I   =  "R_%s"%(mssfr_name) + '_' +'channel I intrinsic [Gpc^{-3} yr^{-1}]'
            namez0_II  =  "R_%s"%(mssfr_name) + '_' +'channel II intrinsic [Gpc^{-3} yr^{-1}]'
            namez0_III =  "R_%s"%(mssfr_name) +'_' + 'channel III intrinsic [Gpc^{-3} yr^{-1}]'
            namez0_IV  =  "R_%s"%(mssfr_name) + '_' +'channel IV intrinsic [Gpc^{-3} yr^{-1}]'
            namez0_V   =  "R_%s"%(mssfr_name) + '_' +'channel V intrinsic [Gpc^{-3} yr^{-1}]'

            NAMES.append(namez0)
            NAMES.append(namez0_I)
            NAMES.append(namez0_II)
            NAMES.append(namez0_III)
            NAMES.append(namez0_IV)
            NAMES.append(namez0_V)



            for ii in range(6):
                datas.append(np.zeros(len(redshifts)))




        df = pd.DataFrame(data=datas, index=NAMES).T
        df.index.names = ['redshift']



        df.to_csv(writePath)







#### RUN different simulation summaries : 










INITIALIZE_FormationChannels = True



if INITIALIZE_FormationChannels==True:
    for DCOtype in ['BHNS','BBH', 'BNS']:
        DCOname = DCOname_dict[DCOtype]
        initalize_formationChannels(DCOname)


    print('done')



# runFormationChannels=True 
# if runFormationChannels ==True:
#     for BPS in  ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q' ]:
#         print(BPS)
#         for DCOtype in ['BHNS','BBH', 'BNS']:
#             print('at DCOtype =', DCOtype)
#             writeToRatesFile_FormationChannels(BPSmodelName=BPS, DCOtype=DCOtype)
#             print('done with ', BPS)



