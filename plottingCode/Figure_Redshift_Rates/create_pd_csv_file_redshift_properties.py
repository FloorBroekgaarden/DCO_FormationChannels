import sys

sys.path.append('../../../common_code')
from PostProcessingScripts import * 
from formation_channels import * 
import astropy.stats

##########
# define colors for formation channels 
# channelColorDict = {'classic':'#118AB2', 'stable B no CEE':'orange',  'immediate CE': '#EF476F'  , r'double-core CE':'#073B4C', 'other':'gray', 'vi':'cyan', 'vii':'#FFD166'}
List_formationchannelOptions = ['All',  'classic',  'stable B no CEE',  'immediate CE',  r'double-core CE', 'vi', 'vii', 'other']
ind_formationchannelOptions = [7,  1, 2, 3, 4, 5, 6, 0]
dictFormationChannelIndex =  {List_formationchannelOptions[i]: ind_formationchannelOptions[i] for i in range(len(List_formationchannelOptions))}

# channelColorDict_lighter = {'classic':adjust_lightness(color='#118AB2', amount=1.6),'stable B no CEE':adjust_lightness(color='orange', amount=1.4), 'immediate CE':adjust_lightness(color='#EF476F', amount=1.2),\
                            # r'double-core CE':adjust_lightness(color='#073B4C', amount=1.8), 'other':adjust_lightness(color='gray', amount=1.5),  'vi':adjust_lightness(color='cyan', amount=1.5), 'vii':adjust_lightness(color='#FFD166', amount=1.2)}
# channelList = ['classic', 'stable B no CEE', 'vii',  'immediate CE',  r'double-core CE', 'other'] #, 'vi']
#######

import pandas as pd
from pathlib import Path


def obtain_redshiftsruns(pathData = '/Volumes/SimonsFoundation/DataDCO/'):
    BPSmodelName='B'
    DCOtype='BNS'
    path_ = '/Volumes/SimonsFoundation/DataDCO/' + alphabetDirDict[BPSmodelName] +'/'
    path  = path_ + 'COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'
    fdata = h5.File(path, 'r')
    redshifts = fdata['redshifts']['redshift'][...].squeeze()
    fdata.close()
    return redshifts 





# adjustedChannelList = ['classic', 'stable B no CEE',  'immediate CE',  r'double-core CE', 'other']
pathData='/Volumes/SimonsFoundation/DataDCO/' # path to datafiles 
redshifts_runs = obtain_redshiftsruns(pathData = pathData) # obtain redshifts that were run 
# print('available redshifts:', redshifts_runs) 
# print(MSSFRnameslist)





def create_pd_redshift_from_xparam_for_fc(DCOtype='BHNS', BPSmodelName='A', pathData='/Volumes/SimonsFoundation/DataDCO/', quantile_values=[0.5, 0.25, 0.75],\
                                               weights_type='merger'):
    """
    whichplot='rate', 'ratio'
    
    """
    redshifts_runs = obtain_redshiftsruns(pathData = pathData) # obtain redshifts that were run    
    adjustedChannelList, DCOname = dict_channel_list[DCOtype], DCOname_dict[DCOtype]

    full_data_path = pathData + alphabetDirDict[BPSmodelName] +'/COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'     # path for files 
    fdata = h5.File(full_data_path,'r')     # read in data 
    
    massCO_ZAMSM1 = fdata['doubleCompactObjects']['M1'][...].squeeze()
    massCO_ZAMSM2 = fdata['doubleCompactObjects']['M2'][...].squeeze()
    # M1 will be the most massive, M2 the least massive compact object. 
    massCO_LVKM1, massCO_LVKM2 = obtainM1BHandM2BHassymetric(m1=fdata['doubleCompactObjects']['M1'][...].squeeze(), m2=fdata['doubleCompactObjects']['M2'][...].squeeze()) 
    MassRatioCO_LVK = massCO_LVKM2/massCO_LVKM1
    
    channels = fdata['doubleCompactObjects']['formaton channel'][...].squeeze()
    
    xparams_to_run = ['t_delay']#['mass_1_LVK', 'mass_2_LVK',  'mass_tot', 't_delay', 'mass_ratio_LVK', 'qZAMS', 'separationInitial', 'M1ZAMS', 'M2ZAMS']
    for xparam in xparams_to_run: 
        print('running xparam ', xparam)
        
        pd_file_path = './formation_median/'+ xparam + '_' + DCOtype + '_' + BPSmodelName + '_' + xparam + '_w_' + weights_type + '.csv'

        # make this a dictionary instead !!
        if xparam=='chirp_mass_LVK':
            param_x = chirpmass(massCO_LVKM1, massCO_LVKM2)
        elif xparam=='mass_tot':
            param_x = massCO_LVKM1 + massCO_LVKM2
        elif xparam=='mass_ratio_LVK':
            param_x = MassRatioCO_LVK
        elif xparam=='mass_1_LVK':
            param_x = massCO_LVKM1
        elif xparam=='mass_2_LVK':
            param_x = massCO_LVKM2
        elif xparam=='log10_t_delay':
            param_x = (fdata['doubleCompactObjects']['tform'][...].squeeze() +  fdata['doubleCompactObjects']['tc'][...].squeeze() ) / 1000 # divide by 1000 to make it in [Gyr]
            param_x = np.log10(param_x)
        elif xparam=='t_delay':
            param_x = (fdata['doubleCompactObjects']['tform'][...].squeeze() +  fdata['doubleCompactObjects']['tc'][...].squeeze() ) / 1000 # divide by 1000 to make it in [Gyr]
        elif xparam=='M1ZAMS':
            param_x = fdata['doubleCompactObjects']['M1ZAMS'][...].squeeze()
        elif xparam=='M2ZAMS':
            param_x = fdata['doubleCompactObjects']['M2ZAMS'][...].squeeze()
        elif xparam=='qZAMS':
            param_x = fdata['doubleCompactObjects']['M2ZAMS'][...].squeeze() / fdata['doubleCompactObjects']['M1ZAMS'][...].squeeze()
        elif xparam=='separationInitial':
            param_x = fdata['doubleCompactObjects']['separationInitial'][...].squeeze()
        del massCO_ZAMSM1
        del massCO_ZAMSM2
        del massCO_LVKM2
        del massCO_LVKM1
        del MassRatioCO_LVK
    
        nZ = len(redshifts_runs)
        fc_data = np.zeros((nZ, 1)) # create empty dataset for redshifts 
        fc_data[:,0] = np.round(redshifts_runs,4)  # fill with redshift data rounded to 4 digits 
        df = pd.DataFrame(fc_data, columns=["redshift"])

    
        for nrC, Channel in enumerate(adjustedChannelList): 
            # obtain fc_mask for the requested channel name 
            print('at Channel: ', Channel)
            ind_wanted = dictFormationChannelIndex[Channel]
            mask_MRR = (channels==ind_wanted)
            if Channel=='stable B no CEE': # add subchannel vii (case A MT to this channel)
                mask_MRR=  (channels==ind_wanted) | (channels==6)
                print('adding case A MT systems (channel 6)')
            elif Channel=='classic': # add subchannel vi (case A MT to this channel)
                mask_MRR=  (channels==ind_wanted) | (channels==5)
                print('adding case A MT systems (channel 5)')
            for ind_mssfr, mssfr in enumerate(MSSFRnameslist[1:]):

                median_at_redshifts = np.zeros_like(redshifts_runs)

                data_to_add = np.zeros((nZ, len(quantile_values)))
                data_to_add[:] = np.nan # start with nan values so that we do not plot datapoint if it doesnt exist 

                column_names = [Channel + 'xyz_' + mssfr+ ' q_' + str(quantile_values[i]) for i in range(len(quantile_values))] # creates array with header names 

                for z_ind, redshift in enumerate(redshifts_runs[:]):
                    redshift = np.round(redshift,4)

                    if weights_type=='merger':
                        fparam_key = 'weights_intrinsicPerRedshift'
                        weightheader = 'w_' + mssfr + '_z_' +  str(redshift)
                    elif weights_type=='formation':
                        fparam_key = "weights_intrinsicFormationPerRedshift"
                        weightheader = 'wform_' + mssfr + '_z_' +  str(redshift)
                    weights_ = fdata[fparam_key][weightheader][...].squeeze()    
                    
                    
                    # we do not care about the very small rates, so do not report medians for formation channels at redshifts that dont pass the minimum channel contribution
                    if np.sum(weights_)>1E-4:  # and we dont want to calculate medians when the total rate is just super low 
                        fractional_contribution_fc = np.sum(weights_[mask_MRR])/np.sum(weights_) # contribution of the channel at this redshift 
                        if fractional_contribution_fc>=minimum_fractional_contribution_fc: 
                            data_to_add[z_ind] = weighted_quantile(values=param_x[mask_MRR], quantiles=quantile_values, sample_weight=weights_[mask_MRR])
#                     else:
#                         print('skipping channel', Channel, ' because fractional contribution is ', fractional_contribution_fc)

                        
                df_to_add = pd.DataFrame(data_to_add, columns=column_names)
                df = pd.concat([df, df_to_add], axis=1)        
        

        df.to_csv(pd_file_path)
        print('saved data, ', DCOtype, xparam, )
    fdata.close()

    return 




def create_pd_redshift_from_xparam_for_total(DCOtype='BHNS', BPSmodelName='A', pathData='/Volumes/SimonsFoundation/DataDCO/', quantile_values=[0.5, 0.25, 0.75],  weights_type='merger'):
    """
    whichplot='rate', 'ratio'
    
    """
    
    print('doing total median')

    redshifts_runs = obtain_redshiftsruns(pathData = pathData) # obtain redshifts that were run    
    adjustedChannelList, DCOname = dict_channel_list[DCOtype], DCOname_dict[DCOtype]

    full_data_path = pathData + alphabetDirDict[BPSmodelName] +'/COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'     # path for files 
    fdata = h5.File(full_data_path,'r')     # read in data 
    
    massCO_ZAMSM1 = fdata['doubleCompactObjects']['M1'][...].squeeze()
    massCO_ZAMSM2 = fdata['doubleCompactObjects']['M2'][...].squeeze()
    # M1 will be the most massive, M2 the least massive compact object. 
    massCO_LVKM1, massCO_LVKM2 = obtainM1BHandM2BHassymetric(m1=fdata['doubleCompactObjects']['M1'][...].squeeze(), m2=fdata['doubleCompactObjects']['M2'][...].squeeze()) 
    MassRatioCO_LVK = massCO_LVKM2/massCO_LVKM1
    
    channels = fdata['doubleCompactObjects']['formaton channel'][...].squeeze()
    weights_ =  fdata['doubleCompactObjects']['weight'][...].squeeze()
    
    xparams_to_run = ['t_delay'] #, 'mass_2_LVK',  'mass_tot', 't_delay', 'mass_ratio_LVK', 'qZAMS', 'separationInitial', 'M1ZAMS', 'M2ZAMS']
    for xparam in xparams_to_run: 

        pd_file_path = './formation_median/'+ xparam + '_' + DCOtype + '_' + BPSmodelName + '_' + xparam + '_w_' + weights_type + 'total.csv'

        # make this a dictionary instead !!
        if xparam=='chirp_mass_LVK':
            param_x = chirpmass(massCO_LVKM1, massCO_LVKM2)
        elif xparam=='mass_tot':
            param_x = massCO_LVKM1 + massCO_LVKM2
        elif xparam=='mass_ratio_LVK':
            param_x = MassRatioCO_LVK
        elif xparam=='mass_1_LVK':
            param_x = massCO_LVKM1
        elif xparam=='mass_2_LVK':
            param_x = massCO_LVKM2
        elif xparam=='log10_t_delay':
            param_x = (fdata['doubleCompactObjects']['tform'][...].squeeze() +  fdata['doubleCompactObjects']['tc'][...].squeeze() ) / 1000 # divide by 1000 to make it in [Gyr]
            param_x = np.log10(param_x)
        elif xparam=='t_delay':
            param_x = (fdata['doubleCompactObjects']['tform'][...].squeeze() +  fdata['doubleCompactObjects']['tc'][...].squeeze() ) / 1000 # divide by 1000 to make it in [Gyr]
        elif xparam=='M1ZAMS':
            param_x = fdata['doubleCompactObjects']['M1ZAMS'][...].squeeze()
        elif xparam=='M2ZAMS':
            param_x = fdata['doubleCompactObjects']['M2ZAMS'][...].squeeze()
        elif xparam=='qZAMS':
            param_x = fdata['doubleCompactObjects']['M2ZAMS'][...].squeeze() / fdata['doubleCompactObjects']['M1ZAMS'][...].squeeze()
        elif xparam=='separationInitial':
            param_x = fdata['doubleCompactObjects']['separationInitial'][...].squeeze()


        nZ = len(redshifts_runs)
        fc_data = np.zeros((nZ, 1)) # create empty dataset for redshifts 
        fc_data[:,0] = np.round(redshifts_runs,4)  # fill with redshift data rounded to 4 digits 
        df = pd.DataFrame(fc_data, columns=["redshift"])

        

#         for ind_mssfr, mssfr in enumerate([MSSFRnameslist[1]]):
        for ind_mssfr, mssfr in enumerate(MSSFRnameslist[1:3]):

            median_at_redshifts = np.zeros_like(redshifts_runs)

            data_to_add = np.zeros((nZ, len(quantile_values)))
            data_to_add[:] = np.nan # start with nan values so that we do not plot datapoint if it doesnt exist 
            column_names = ['all_' + 'xyz_' + mssfr+ ' q_' + str(quantile_values[i]) for i in range(len(quantile_values))] # creates array with header names 

            for z_ind, redshift in enumerate(redshifts_runs[0:]):
                redshift = np.round(redshift,4)

                if weights_type=='merger':
                    fparam_key = 'weights_intrinsicPerRedshift'
                    weightheader = 'w_' + mssfr + '_z_' +  str(redshift)
                elif weights_type=='formation':
                    fparam_key = "weights_intrinsicFormationPerRedshift"
                    weightheader = 'wform_' + mssfr + '_z_' +  str(redshift)
                weights_ = fdata[fparam_key][weightheader][...].squeeze()    

                # we do not care about the very small rates, keep those to zero (nan)
                if np.sum(weights_)>1E-4: 
                    data_to_add[z_ind] = weighted_quantile(values=param_x, quantiles=quantile_values, sample_weight=weights_)

            df_to_add = pd.DataFrame(data_to_add, columns=column_names)
            df = pd.concat([df, df_to_add], axis=1)        


        df.to_csv(pd_file_path)
        print('saved data, ', DCOtype, xparam, )
    fdata.close()

    return 







################ CHANGE THE THINGS BELOW ################
quantile_values=[0.5, 0.25, 0.75, 0.1, 0.9]
DCOTypeList = ['BBH'] #, 'BHNS'] #['BHNS', 'BNS', 'BBH'] #, https://arxiv.org/pdf/2402.00935.pdf
pathData='/Volumes/SimonsFoundation/DataDCO/'
single_model=True
# only_channels_with_min_contribution = 0.01 # percent
dict_channel_list = {'BBH':['classic', 'stable B no CEE', 'other', 'immediate CE', r'double-core CE'],\
                     'BHNS':['classic', 'stable B no CEE', 'other', 'immediate CE', r'double-core CE'],\
                     'BNS':['classic', 'stable B no CEE', 'other', 'immediate CE', r'double-core CE'] } 
weights_type= 'formation'# 'merger'# 'formation' #
create_df_file = True
minimum_fractional_contribution_fc = 0.05
###########################################################


        
# TOTAL MEDIAN 

def run_create_median(RunCreatePDfilesQuantilesDCO_form_channels=False, RunCreatePDfilesQuantilesDCO_totals=False):

    
    if RunCreatePDfilesQuantilesDCO_form_channels==True: 
        for DCOtype in DCOTypeList:
            print()
            print('------------------------------')
            print('at ',  DCOtype)
            for ind_bps, BPSmodelName in enumerate(BPSnameslist[8:10]):
#             for ind_bps, BPSmodelName in enumerate([BPSnameslist[0]]):
                print('at BPS model ', BPSmodelName)

                df = create_pd_redshift_from_xparam_for_fc(quantile_values=quantile_values, BPSmodelName=BPSmodelName, DCOtype=DCOtype, weights_type=weights_type) #, pd_file_path=pd_file_path)

    
    
    if RunCreatePDfilesQuantilesDCO_totals==True:
        for DCOtype in DCOTypeList:
            print()
            print('------------------------------')
            print('at ',  DCOtype)
            for ind_bps, BPSmodelName in enumerate([BPSnameslist[0]]):
#             for ind_bps, BPSmodelName in enumerate(BPSnameslist):
                print('at BPS model ', BPSmodelName)

                df = create_pd_redshift_from_xparam_for_total(quantile_values=quantile_values, BPSmodelName=BPSmodelName, DCOtype=DCOtype, weights_type=weights_type) #, pd_file_path=pd_file_path)






    return print('DONE')


# run create pandas efficiently 
run_create_median(RunCreatePDfilesQuantilesDCO_form_channels=True)
# run totals:
# run_create_median(RunCreatePDfilesQuantilesDCO_totals=True) 
    
    
    
