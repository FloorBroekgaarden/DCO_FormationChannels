
import numpy as np
import sys
import h5py as h5

sys.path.append('../../common_code/')


from PostProcessingScripts import * 
from formation_channels import * 

import ClassCOMPAS     as CC ###


import pandas as pd
from astropy import units as u
from astropy import constants as const






def createEmptyCSVplaceholder(DCOtype='BBH'):


    DCOname = DCOname_dict[DCOtype]

    stringgg =  'formation_channels_per_Z_hdf5weights'
    # stringgg =  'formation_channels_per_Z'
    writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_1/'  + stringgg + '_'  + DCOname + '.csv' 


    BPSnameslist = list(string.ascii_uppercase)[0:26]
    headerDict_Z  = { 5:'channel VI',     6:'channel VII',    7:'All',     0:'channel V', 1:'channel I', 2:'channel II', 3:'channel III', 4:'channel IV'}    
    enumerate_list = range(8)


    # make list with header names 
    NAMES = []
    for ind_m, m_ in enumerate(BPSnameslist):
        for ind_c, c_ in enumerate(enumerate_list):
            str_ = m_ + ' ' + headerDict_Z[c_]+ '  [Msun^{-1}]'
            NAMES.append(str_)

            
    # metallicities of simulation 
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

    # make dataframe with zeroes as placeholders 
    for i in range(len(NAMES)):
        datas.append(np.zeros(nMetallicities))
 
    
    df = pd.DataFrame(data=datas, index=NAMES, columns=Zlist).T
    df.columns =   df.columns.map(str)
    df.index.names = ['Z_i']
    df.columns.names = ['model']


    df.to_csv(writePath)
    
    return 










def writeFormationRatesAndChannelsToFile(DCOtype='BBH', BPSmodelName='Z', runquick=True):
    """
    simple function that opens the data file for given DCO type and BPS model, and then calculates for each metallicity the total yield
    for each formation channel, as well as the total yield. This uses the weights from the hdf5 file (stroopwafel weights), as well as 
    the calculated formation channels from the DoubleCompactObjects file (though you can self-consistently calculate the formation channels
    by uncommenting the channels = identify_formation_channels() instead)
    It returns a csv file with the saved calculated yields per formation channel and metallicity. 



    """

    DCOname = DCOname_dict[DCOtype]


    # path for files 

    path_ = '/Volumes/SimonsFoundation/DataDCO/'+ alphabetDirDict[BPSmodelName] +'/'
    path  = path_ + 'COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'
    stringgg =  'formation_channels_per_Z_hdf5weights'
    writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_1/'  + stringgg + '_'  + DCOname + '.csv' 


    # read in data 
    fdata = h5.File(path, 'r')



    seeds = fdata['doubleCompactObjects']['seed'][...].squeeze()
    channels = fdata['doubleCompactObjects']["formaton channel"][...].squeeze()  # this is faster, but below will always work, even if fc are not ready yet 
    # channels = identify_formation_channels(seeds=seeds, file=fdata) # update hdf5 weights 
    metallicities = fdata['doubleCompactObjects']['Metallicity1'][...].squeeze()
    weights = fdata['doubleCompactObjects']['weight'][...].squeeze()

    headerDict_Z  = { 5:'channel VI',     6:'channel VII',    7:'All',     0:'channel V', 1:'channel I', 2:'channel II', 3:'channel III', 4:'channel IV'}    
    enumerate_list = range(8)

    # temp = range(nModels+3)




    print('now at DCO type  ', DCOtype)

    Zlist_labelnames=['0_0001','0_00011', '0_00012', '0_00014', '0_00016', '0_00017',\
    '0_00019', '0_00022', '0_00024', '0_00027', '0_0003', '0_00034',\
    '0_00037', '0_00042', '0_00047', '0_00052', '0_00058', '0_00065',\
    '0_00073', '0_00081', '0_0009', '0_00101', '0_00113', '0_00126',\
    '0_0014', '0_00157', '0_00175', '0_00195', '0_00218', '0_00243',\
    '0_00272', '0_00303', '0_00339', '0_00378', '0_00422', '0_00471',\
    '0_00526', '0_00587', '0_00655', '0_00732', '0_00817', '0_00912',\
    '0_01018', '0_01137', '0_01269', '0_01416', '0_01581', '0_01765', '0_01971', '0_022', '0_0244', '0_02705', '0_03']

    Zlist=[0.0001, 0.00011, 0.00012, 0.00014, 0.00016, 0.00017,\
           0.00019, 0.00022, 0.00024, 0.00027, 0.0003, 0.00034, \
           0.00037, 0.00042, 0.00047, 0.00052, 0.00058, 0.00065,\
           0.00073, 0.00081, 0.0009, 0.00101, 0.00113, 0.00126,\
           0.0014, 0.00157, 0.00175, 0.00195, 0.00218, 0.00243, \
           0.00272, 0.00303, 0.00339, 0.00378, 0.00422, 0.00471, \
           0.00526, 0.00587, 0.00655, 0.00732, 0.00817, 0.00912, \
           0.01018, 0.01137, 0.01269, 0.01416, 0.01581, 0.01765, 0.01971, 0.022, 0.0244, 0.02705, 0.03]



    print()
    # print('now at model ', alphabetDirDict[bps_model])

    # set always optimistic CE false, unless we are doing the optimistic variation
    OPTIMISTIC=False
    if (BPSmodelName=='F') | (BPSmodelName=='K'):
        OPTIMISTIC=True
        print('doing optimistic version of fiducial')

    # path to datafile 
    # path = pathCOMPASOutput+alphabetDirDict[bps_model] + '/' + 'COMPASCompactOutput_'+DCOtype +'_'+bps_model+'.h5'
    




    if runquick==True: 
        print('running quick Data calculation')
        meanMassEvolved = 77708655 # hack to save time, real calculation can be done by using the elif option below 
    
        for nrC, Channel in enumerate(enumerate_list):  

            print('at channel ', Channel)      

            formationRate = np.zeros(len(Zlist))  

            mask_C  = (channels==Channel)

            ind_inMetallicityGrid = 0
            for nrZ, Z in enumerate(Zlist):

                maskZ = (metallicities == Z)
                mask_systems = (mask_C==1) & (maskZ==1)

                # -1 channel is the total rate, so we skip the channel mask
                if nrC==enumerate_list[-1]:
                    formationRate[nrZ]         = np.sum(weights[maskZ])/ (meanMassEvolved) 
                # channel specific rate
                else:
                    formationRate[nrZ]         = (np.sum(weights[mask_systems])) / (meanMassEvolved)

                ind_inMetallicityGrid +=1


            df = pd.read_csv(writePath, index_col=0)
            str_ = BPSmodelName + ' ' + headerDict_Z[nrC]+ '  [Msun^{-1}]'
            df[str_] = formationRate
            df.to_csv(writePath)



    else:
        path = path_ + '/' + 'COMPASOutput.h5'
        print('doing difficult Data calculation')
        #But I want only within Hubble time 
        Data            = CC.COMPASData(path=path, lazyData=True, Mlower=5., \
                         Mupper=150, binaryFraction=1)
        Data.setCOMPASDCOmask(types=DCOtype,  withinHubbleTime=True, optimistic=OPTIMISTIC)
        Data.setCOMPASData()
        metallicityGrid = Data.metallicityGrid
        
        totalMassEvolvedPerZ = Data.totalMassEvolvedPerZ
        print('totalMassEvolvedPerZ', totalMassEvolvedPerZ)
        meanMassEvolved = np.mean(totalMassEvolvedPerZ)
        print('meanMassEvolved', meanMassEvolved)
        del Data
        del totalMassEvolvedPerZ 

        for nrC, Channel in enumerate(enumerate_list):  

            print('at channel ', Channel)      

            formationRate = np.zeros(len(Zlist))  

            mask_C  = (channels==Channel)

            ind_inMetallicityGrid = 0
            for nrZ, Z in enumerate(Zlist):


                if Z in metallicityGrid:


                    maskZ = (metallicities == Z)
                    mask_systems = (mask_C==1) & (maskZ==1)

                    # -1 channel is the total rate, so we skip the channel mask
                    if nrC==enumerate_list[-1]:
                        formationRate[nrZ]         = np.sum(weights[maskZ])/ (meanMassEvolved) 
                    # channel specific rate
                    else:
                        formationRate[nrZ]         = (np.sum(weights[mask_systems])) / (meanMassEvolved)

                    ind_inMetallicityGrid +=1
                # no DCO systems at this metallicty
                else:
                    formationRate[nrZ]         = 0

            df = pd.read_csv(writePath, index_col=0)
            str_ = BPSmodelName + ' ' + headerDict_Z[nrC]+ '  [Msun^{-1}]'
            df[str_] = formationRate
            df.to_csv(writePath)


    print('finished')

    return









INITIALIZE_FormationChannels_per_Z = False #True




if INITIALIZE_FormationChannels_per_Z==True:
    for DCOtype in ['BHNS','BBH', 'BNS']:
        # DCOname = DCOname_dict[DCOtype]
        createEmptyCSVplaceholder(DCOtype=DCOtype)


    print('done')



runFormationChannels_per_Z=True



if runFormationChannels_per_Z==True:
    for BPS in  BPSnameslist[:]:
        print(BPS)
        for DCOtype in [ 'BNS', 'BHNS', 'BBH']: #'BHNS',, 'BBH'
            print('at DCOtype =', DCOtype)
            writeFormationRatesAndChannelsToFile(BPSmodelName=BPS, DCOtype=DCOtype, runquick=True)
            print('done with ', BPS)













