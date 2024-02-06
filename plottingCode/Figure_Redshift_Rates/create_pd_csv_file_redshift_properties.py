import sys

sys.path.append('../../../common_code')
from PostProcessingScripts import * 
from formation_channels import * 
import astropy.stats

# from IPython.display import display, HTML
# display(HTML("<style>.container { width:100% !important; }</style>"))

##########
# define colors for formation channels 
channelColorDict = {'classic':'#118AB2', 'stable B no CEE':'orange',  'immediate CE': '#EF476F'  , r'double-core CE':'#073B4C', 'other':'gray', 'vi':'cyan', 'vii':'#FFD166'}
List_formationchannelOptions = ['All',  'classic',  'stable B no CEE',  'immediate CE',  r'double-core CE', 'vi', 'vii', 'other']
ind_formationchannelOptions = [7,  1, 2, 3, 4, 5, 6, 0]
dictFormationChannelIndex =  {List_formationchannelOptions[i]: ind_formationchannelOptions[i] for i in range(len(List_formationchannelOptions))}

channelColorDict_lighter = {'classic':adjust_lightness(color='#118AB2', amount=1.6),'stable B no CEE':adjust_lightness(color='orange', amount=1.4), 'immediate CE':adjust_lightness(color='#EF476F', amount=1.2),\
                            r'double-core CE':adjust_lightness(color='#073B4C', amount=1.8), 'other':adjust_lightness(color='gray', amount=1.5),  'vi':adjust_lightness(color='cyan', amount=1.5), 'vii':adjust_lightness(color='#FFD166', amount=1.2)}
channelList = ['classic', 'stable B no CEE', 'vii',  'immediate CE',  r'double-core CE', 'other'] #, 'vi']
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





adjustedChannelList = ['classic', 'stable B no CEE', 'vii', 'immediate CE',  r'double-core CE', 'other']
pathData='/Volumes/SimonsFoundation/DataDCO/' # path to datafiles 
redshifts_runs = obtain_redshiftsruns(pathData = pathData) # obtain redshifts that were run 
print('available redshifts:', redshifts_runs) 

print(MSSFRnameslist)





def create_pd_redshift_from_xparam(DCOtype='BHNS', BPSmodelName='A', pathData='/Volumes/SimonsFoundation/DataDCO/', quantile_values=[0.5, 0.25, 0.75], xparam='log10_t_delay', weights_type='merger', pd_file_path='path_to_pd_file_to_create'):
    """
    whichplot='rate', 'ratio'
    
    """

    pathData='/Volumes/SimonsFoundation/DataDCO/' # path to datafiles 
    redshifts_runs = obtain_redshiftsruns(pathData = pathData) # obtain redshifts that were run    
    adjustedChannelList, DCOname = dict_channel_list[DCOtype], DCOname_dict[DCOtype]

    full_data_path = pathData + alphabetDirDict[BPSmodelName] +'/COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'     # path for files 
    fdata = h5.File(full_data_path,'r')     # read in data 
    
    massCO_ZAMSM1 = fdata['doubleCompactObjects']['M1'][...].squeeze()
    massCO_ZAMSM2 = fdata['doubleCompactObjects']['M2'][...].squeeze()
    # M1 will be the most massive, M2 the least massive compact object. 
    massCO_LVKM1, massCO_LVKM2 = obtainM1BHandM2BHassymetric(m1=fdata['doubleCompactObjects']['M1'][...].squeeze(), m2=fdata['doubleCompactObjects']['M2'][...].squeeze()) 
    MassRatioCO_LVK = massCO_LVKM2/massCO_LVKM1
    
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


    
    channels = fdata['doubleCompactObjects']['formaton channel'][...].squeeze()
    weights_ =  fdata['doubleCompactObjects']['weight'][...].squeeze()
    #     param_x = np.log10(param_x)
    
    
    nZ = len(redshifts_runs)
    fc_data = np.zeros((nZ, 1)) # create empty dataset for redshifts 
    fc_data[:,0] = np.round(redshifts_runs,4)  # fill with redshift data rounded to 4 digits 
    df = pd.DataFrame(fc_data, columns=["redshift"])

    # HERE 

            
            
    
    for nrC, Channel in enumerate(adjustedChannelList): 
        # obtain fc_mask for the requested channel name 
        ind_wanted = dictFormationChannelIndex[Channel]
        mask_MRR = (channels==ind_wanted)
        print('now at channel: ', Channel)
        
        
        for ind_mssfr, mssfr in enumerate([MSSFRnameslist[1]]):

            median_at_redshifts = np.zeros_like(redshifts_runs)

            

            data_to_add = np.zeros((nZ, len(quantile_values)))
            data_to_add[:] = np.nan # start with nan values so that we do not plot datapoint if it doesnt exist 
            column_names = [Channel + 'xyz_' + mssfr+ ' q_' + str(quantile_values[i]) for i in range(len(quantile_values))] # creates array with header names 
            
            
            for z_ind, redshift in enumerate(redshifts_runs[0:]):
                redshift = np.round(redshift,4)
      
    
                if weights_type=='merger':
                    fparam_key = 'weights_intrinsicPerRedshift'
                    weightheader = 'w_' + mssfr + '_z_' +  str(redshift)
                    weights_ = fdata[fparam_key][weightheader][...].squeeze()
                elif weights_type=='formation':
                    fparam_key = "weights_intrinsicFormationPerRedshift"
                    weightheader = 'wform_' + mssfr + '_z_' +  str(redshift)
                    weights_ = fdata[fparam_key][weightheader][...].squeeze()    
    
                # we do not care about the very small rates, keep those to zero (nan)
                if np.sum(weights_[mask_MRR])>1E-6: 
                    data_to_add[z_ind] = weighted_quantile(values=param_x[mask_MRR], quantiles=quantile_values, sample_weight=weights_[mask_MRR])

            df_to_add = pd.DataFrame(data_to_add, columns=column_names)
            df = pd.concat([df, df_to_add], axis=1)

    # always close the dataset
    fdata.close()
    
    
    df.to_csv(pd_file_path)

    return df




def plot_xparam_formation_channels_redshift_for_quantiles(axe='None', DCOtype='BHNS', BPS_models_to_run_list=['A'], 
                                                  pathData='/Volumes/SimonsFoundation/DataDCO/',\
                                                mask_specific_mssfr=None,\
                                                whichQuantity='median', value_for_fraction=False, \
                                              add_model_label=True, quantile_values=[0.5, 0.25, 0.75], xparam='log10_t_delay' , weights_type='merger', single_model=True, create_df_file=True):
    
    
    
    adjustedChannelList, DCOname = dict_channel_list[DCOtype], DCOname_dict[DCOtype]
    
    
    pd_file_path = './formation_median/'+ xparam_wanted + '_' + DCOtype + '_' + DCOtype + '_' + xparam_wanted + '_w_' + weights_type + '.csv'
    for ind_bps, BPSmodelName in enumerate(BPS_models_to_run_list):
        if create_df_file==True: df = create_pd_redshift_from_xparam(quantile_values=quantile_values, BPSmodelName=BPSmodelName, DCOtype=DCOtype, xparam=xparam, weights_type=weights_type, pd_file_path=pd_file_path)
        else: df = pd.read_csv(pd_file_path)
        redshifts = df["redshift"]
        
        

        # plot the channel 
        for nrC, Channel in enumerate(adjustedChannelList): 
            
            for ind_mssfr, mssfr in enumerate([MSSFRnameslist[1]]):
            
                c_FC = channelColorDict[Channel]
                colors_lighter_FC =  channelColorDict_lighter[Channel]

                column_names = [Channel + 'xyz_' + mssfr+ ' q_' + str(quantile_values[i]) for i in range(len(quantile_values))] # creates array with header names 
                qvalues = [df[column_names[i]] for i in range(len(quantile_values))]

                axe.scatter((redshifts), qvalues[0].values, color=c_FC, marker=dictMarkerShape[BPSmodelName], s=40) #/norm_classic_tdelay
                axe.plot(   (redshifts), qvalues[0].values, color=c_FC, lw=4) #, ls=linestyles_mssfrind[ind_mssfr_zind])

#                 axe.fill_between((redshifts), y1=qvalues[1], y2=qvalues[0], color=colors_lighter_FC, alpha=0.5)
#                 axe.fill_between((redshifts), y1=qvalues[0], y2=qvalues[2], color=colors_lighter_FC, alpha=0.5)

#                 axe.fill_between((redshifts), y1=qvalues[3], y2=qvalues[1], color=colors_lighter_FC, alpha=0.2)
#                 axe.fill_between((redshifts), y1=qvalues[2], y2=qvalues[4], color=colors_lighter_FC, alpha=0.2)



    xlabel = r'\textbf{redshift} $z$'
    
    if xparam=='chirp_mass_LVK':
        ylabel = r'$\mathcal{M}_{\rm{c}} \ [M_{\odot}]$'
    elif xparam=='mass_tot':
        ylabel = r'$\rm{M}_{\rm{tot}} \ [M_{\odot}]$'
    elif xparam=='mass_ratio_LVK':
        ylabel = r'$q$'
    elif xparam=='mass_1_LVK':
        ylabel = r'$m_1 [M_{\odot}]$'
    elif xparam=='mass_2_LVK':
        ylabel = r'$m_2 [M_{\odot}]$'
    elif xparam=='log10_t_delay':
        ylabel = r'$\log_{10} t_{\rm{delay}} \ [\rm{Gyr}]$'
    elif xparam=='t_delay':
        ylabel = r'$t_{\rm{delay}} \ [\rm{Gyr}]$'
    elif xparam=='M1ZAMS':
        ylabel = r'$\rm{M}_{\rm{ZAMS, 1}} \ [M_{\odot}]$'
    elif xparam=='M2ZAMS':
        ylabel = r'$\rm{M}_{\rm{ZAMS, 2}} \ [M_{\odot}]$'
    elif xparam=='qZAMS':
        ylabel = r'${q}_{\rm{ZAMS}} $'
    elif xparam=='separationInitial':
        ylabel = r'${a}_{\rm{ZAMS}} [\rm{AU}] $'      
    

    # axes properties 
    axe = layoutAxes(axe, nameX=xlabel, nameY=ylabel, setMinor=True)
    axe.set_xlim(0,9.7)  # redshift range 
    if xparam in ['t_delay', 'separationInitial']:
        axe.set_yscale('log')
    
    ## add label in a legend for plot 
    if single_model==True: annotate_label = r'\textbf{model %s:}'%(BPSmodelName) +'\n' + alphabetPhysicalNameDict[BPSmodelName]
    else: annotate_label = r'\textbf{%s}'%(DCOtype)        
    bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.95)
    axe.annotate(annotate_label, xy=(0.042, .95), xycoords='axes fraction', fontsize = fs-6, weight = 'bold', ha='left', va="top",bbox=bbox_props, zorder=100)             
            
            
            
    return axe








################ CHANGE THE THINGS BELOW ################
quantile_values=[0.5, 0.25, 0.75, 0.1, 0.9 ]
DCOTypeList = ['BBH'] #, https://arxiv.org/pdf/2402.00935.pdf
whichQuantity = 'median'
pathData='/Volumes/SimonsFoundation/DataDCO/'
single_model=True
boot_strap_uncertainties=False # True 
only_channels_with_min_contribution = 0.01 # percent
dict_channel_list = {'BBH':['classic', 'stable B no CEE'],\
                     'BHNS':['classic', 'stable B no CEE', 'immediate CE'],\
                     'BNS':['classic', r'double-core CE', 'other'] } 
weights_type='merger'
create_df_file = True
###########################################################




# for xparam_wanted in [ 'mass_tot', 't_delay', 'mass_ratio_LVK', 'mass_1_LVK', 'mass_2_LVK','chirp_mass_LVK','qZAMS', 'separationInitial']:
for xparam_wanted in [ 'mass_tot', 't_delay', 'mass_ratio_LVK', 'qZAMS', 'separationInitial', 'M1ZAMS', 'M2ZAMS']: #, 'M1ZAMS', 'M2ZAMS']:
    print('at xparam ', xparam_wanted)
    
    if single_model==True:enumerate_list = BPSnameslist[3:]
    else: enumerate_list = [0]
    for ind_m, BPSmodelName in  enumerate(enumerate_list):
        if single_model==True:
            BPS_models_to_run_list=[BPSmodelName]
            print('at BPS model ', BPSmodelName)
            save_fig_string = BPSmodelName
        else: 
            BPS_models_to_run_list = BPSnameslist
            save_fig_string = 'all'
            
        for DCOtype in DCOTypeList: #'BNS', 'BHNS', 
            print('at DCOtype =', DCOtype)


            ncols, nrows= 1,1
            f, ax= plt.subplots(ncols=ncols,nrows=nrows,figsize=(10,6), gridspec_kw={"width_ratios":1.5*np.ones(ncols), "height_ratios":1*np.ones(nrows)})

            ax = plot_xparam_formation_channels_redshift_for_quantiles(axe=ax, DCOtype=DCOtype, BPS_models_to_run_list=BPS_models_to_run_list,\
                                                                  pathData=pathData,  single_model=single_model, quantile_values=quantile_values,\
                                                                               xparam=xparam_wanted, weights_type=weights_type, create_df_file=create_df_file) 


            
            ##  SAVE FIG  ###
            plt.tight_layout()
            plt.subplots_adjust(wspace=0., hspace=0.18)  
            plt.savefig('./formation_median/'+ xparam_wanted + '/zQuantile_' +  DCOtype + '_' + save_fig_string + '_' + xparam_wanted + '_w_' + weights_type + '.png', transparent=False, dpi=300)
            plt.show()
            plt.close()
            print()











