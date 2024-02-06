#Needed in general
import sys

from KDEpy import FFTKDE 
sys.path.append('../../../common_code')
from PostProcessingScripts import * 
from formation_channels import * 
import astropy.stats


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

    





def plot_FC_distribution(axe='None', xparam='chiEff', BPSmodelName='A', mode='pdf',\
                          spin_threshold='None', bw=0.01, xlim=[0,1], ylim=[0,1],\
                          plotYlog='False', ylim_threshold=0.02,DCOtype='BBH', histtype='kde', Metallicity_lower=False, Metallicity_upper=False, \
                          pathData='/Volumes/SimonsFoundation/DataDCO/', channelList=['classic']):#, mssfr='112'):
    
    fs_l = 28 # label fontsize
    fs_major=34
            






    # path for files 
    path_ = pathData + alphabetDirDict[BPSmodelName] +'/'
    path  = path_ + 'COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'
    
    fdata = h5.File(path, 'r')
    massCO_ZAMSM1 = fdata['doubleCompactObjects']['M1'][...].squeeze()
    massCO_ZAMSM2 = fdata['doubleCompactObjects']['M2'][...].squeeze()
    # M1 will be the most massive, M2 the least massive compact object. 
    massCO_LVKM1, massCO_LVKM2 = obtainM1BHandM2BHassymetric(m1=fdata['doubleCompactObjects']['M1'][...].squeeze(), m2=fdata['doubleCompactObjects']['M2'][...].squeeze()) 
    MassRatioCO_LVK = massCO_LVKM2/massCO_LVKM1
    
    seedsDCO  = fdata['doubleCompactObjects']['seed'][...].squeeze()  # get the seeds in the DCO file 
    seedsSN   = fdata['supernovae']['randomSeed'][...].squeeze()    # get the seeds in the SN file 
    indices   = np.sort(np.unique(seedsSN[1::2], return_index=True)[1])
    maskSNdco = np.in1d(seedsSN,  seedsDCO) # mask in the SNe files the SNe that correspond to our DCO
    whichSN   = fdata['supernovae']['whichStar'][...].squeeze()[maskSNdco]  # this is 1 if the initially primary star goes SN and 2 if the secondary goes supernova     
    whichSN2  = whichSN[1::2][indices]

    # either SN2 = primary (1) and M1 is > M2, or SN2 = secondary & M1 < M2 
    # this takes into account (first term) rejuvenation 
    mask_MRR = ((whichSN2==1) & (massCO_ZAMSM1>massCO_ZAMSM2) ) | ((whichSN2==2) & (massCO_ZAMSM1<massCO_ZAMSM2)) 

    # obtain formation channels 
    seeds = fdata['doubleCompactObjects']['seed'][...].squeeze()
    channels = identify_formation_channels(seeds=seeds, file=fdata)
    Metallicity =  fdata['doubleCompactObjects']['Metallicity1'][...].squeeze()
    
    
    del massCO_ZAMSM1
    del massCO_ZAMSM2
    del whichSN2
    del whichSN
    del maskSNdco
    del indices
    del seedsSN
    del seedsDCO


    
    #if (mode in ['spin_PDF', 'spin_fraction', 'spinOne_fraction', 'spinTwo_fraction', 'm1spin_or_m2spin_fraction', 'spin_CDF' ]) | (xparam in ['chi_of_spinning_BH', 'chi_effective', 'log10_t_delay']):
    
    if xparam in ['spin_1_LVK', 'spin_2_LVK', 'chi_of_spinning_BH', 'chi_effective']:
        spin = COspin(data_path=path, state='he_depletion')  # set class 
        spin.setCOMPASData() # reads in the COMPAS DCO parameters 
        spinZAMSM1, spinZAMSM2  = spin.BaveraSpin() #ZAMS M1 SPIN 

        spinLVKM1, spinLVKM2 = np.zeros_like(spinZAMSM1), np.zeros_like(spinZAMSM1)
        spinLVKM1[mask_MRR] = spinZAMSM2[mask_MRR]  # MRR so M1 comes from M2ZAMS, we assign it spin from M2ZAMS
        spinLVKM1[~mask_MRR] = spinZAMSM1[~mask_MRR]  # no MRR so M1 comes from M1ZAMS, we assign it spin from M1ZAMS
        spinLVKM2[mask_MRR] = spinZAMSM1[mask_MRR]   # MRR so M2 comes from M1ZAMS, we assign it spin from M1ZAMS
        spinLVKM2[~mask_MRR] = spinZAMSM2[~mask_MRR]   # no MRR so M2 comes from M2ZAMS, we assign it spin from M2ZAMS     

        # spin_threshold = 0.05 # definition of "spinning BH"
        mask_LVKM1_spinning = (spinLVKM1 > spin_threshold ) 
        mask_LVKM2_spinning = (spinLVKM2 > spin_threshold ) # definition of "spinning BH"
        mask_anySpin = (spinLVKM1 > spin_threshold ) | (spinLVKM2 > spin_threshold )    

    

    if xparam=='chirp_mass_LVK':
        param_x = chirpmass(massCO_LVKM1, massCO_LVKM2)
        nameX = r'$\mathcal{M}_{\rm{c}} \ [M_{\odot}]$'
        nameY = r'\textbf{PDF}'
        xx = np.linspace(1,100,1000)
        print(param_x)
        print(np.shape(param_x))
        
    elif xparam=='mass_ratio_LVK':
        param_x = MassRatioCO_LVK
        nameX = r'$q$'
        nameY = r'\textbf{PDF}'
        xx = np.linspace(-0.2,1.2,1000)

        

    elif xparam=='mass_1_LVK':
        param_x = massCO_LVKM1
        nameX = r'$m_1 [M_{\odot}]$'
        nameY = r'\textbf{PDF}'
        if DCOtype=='BBH':
            xx = np.linspace(-1,150,1000) # needs to be a little bit larger
        elif DCOtype=='BHNS':
            xx = np.linspace(-1,50,1000)
        elif DCOtype=='BNS':
            xx = np.linspace(-2,5,500)


    elif xparam=='mass_2_LVK':
        param_x = massCO_LVKM2
        nameX = r'$m_2 [M_{\odot}]$'
        nameY = r'\textbf{PDF}'
        if DCOtype=='BBH':
            xx = np.linspace(-1,150,1000) # needs to be a little bit larger
        elif DCOtype=='BHNS':
            xx = np.linspace(-2,5,1000)
        elif DCOtype=='BNS':
            xx = np.linspace(-2,5,500)

    elif xparam=='spin_1_LVK':
        param_x = spinLVKM1
        nameX = r'$\chi_1$'
        nameY = r'\textbf{PDF}'
        xx = np.linspace(-2,2,1000)
        # xx = np.linspace(0,1,100)

    elif xparam=='spin_2_LVK':
        param_x = spinLVKM1
        nameX = r'$\chi_2$'
        nameY = r'\textbf{PDF}'
        xx = np.linspace(-2,2,1000)


    elif xparam=='chi_of_spinning_BH':
        param_x = spinLVKM1 + spinLVKM2

        nameX = r'$\chi_{\rm{i}}$'
        nameY = r'\textbf{PDF}'
        xx = np.linspace(-0.2,2,1000)  

        
    elif xparam=='chi_effective':
        param_x = ((spinLVKM1*massCO_LVKM1) + (spinLVKM2*massCO_LVKM2)) / (massCO_LVKM1+massCO_LVKM2)
        nameX = r'$\chi_{\rm{eff}}$'
        nameY = r'\textbf{PDF}'
        xx = np.linspace(-0.2,2,1000)          
        
    elif xparam=='log10_t_delay':
        param_x = (fdata['doubleCompactObjects']['tform'][...].squeeze() +  fdata['doubleCompactObjects']['tc'][...].squeeze() ) / 1000 # divide by 1000 to make it in [Gyr]
        param_x = np.log10(param_x)
        nameX = r'$t_{\rm{delay}} \ [\rm{Gyr}]$'
        nameY = r'\textbf{PDF}'  
        xx = np.linspace(-4,2,500)  
        print('obtained params')

    
    if (mode=='MRR_PDF') | (mode=='spin_PDF'):
        nameY = r'\textbf{PDF}'
    elif (mode=='MRR_fraction'):
        nameY = r'$\rm{f}_{\rm{channel}}$'
    elif (mode=='notMRR_fraction'):
        nameY = r'$\rm{f}_{\rm{non MRR}}$'
    elif (mode=='spin_fraction'):
        nameY =   r'$\rm{f}_{\chi_{1} > %s}$'%spin_threshold
    elif (mode=='m1spin_or_m2spin_fraction'):
        nameY =   r'$\rm{f}_{\chi_{\rm{i}} > %s}$'%spin_threshold 
    elif (mode=='MRR_CDF') | (mode=='spin_CDF'):
        nameY = r'\textbf{CDF}'
    elif (mode=='spinOne_fraction'):
        nameY =   r'$\rm{f}_{\chi_{1} > %s}$'%spin_threshold
    elif (mode=='spinTwo_fraction'):
        nameY =   r'$\rm{f}_{\chi_{2} > %s}$'%spin_threshold
    else:
        raise ValueError("the provided `mode` is not recognized, the value given is %s"%mode)
        
    
      

        
    adjustedChannelList = ['classic', 'stable B no CEE', 'vii', 'immediate CE',  r'double-core CE', 'other']
    
    # adjustedChannelList = ['classic', 'stable B no CEE', 'immediate CE',  'other']#, 'stable B no CEE', 'immediate CE']#
    for nrC, Channel in enumerate(channelList): 
        

        yy_max = np.zeros_like(xx)
        yy_min = np.ones_like(xx) *5

        yy_max_hist = np.zeros_like(xx[1:])
        yy_min_hist = np.ones_like(xx[1:]) *5
        # yy_max_tot = np.zeros_like(xx)
        # yy_min_tot = np.ones_like(xx) *5

        print('now at channel: ', Channel)
        ind_wanted = dictFormationChannelIndex[Channel]
        # set color         
        c_FC = channelColorDict[Channel]
        colors_lighter_FC =  channelColorDict_lighter[Channel]

        colors_metallicity_bins =  sns.color_palette("husl", 7)

        metallicity_alphas_list = np.linspace(0.3,1, 6) 
        # metallicity_linestyles_list = ['solid', 'long dash with offset', 'dashed',  'dashdot', 'loosely dashed' ':']
        metallicity_linestyles_list = [ 'dotted', 'dotted','dashed', 'dashed', 'solid', 'solid']        


        for ind_Zbin, lowerZbin in enumerate(Metallicity_lower):
            upperZbin = Metallicity_upper[ind_Zbin]    
            # print('running for Zbin: [',lowerZbin, ',', upperZbin, ']' )        
            mask_MRR = (channels==ind_wanted) & (Metallicity>=lowerZbin) & (Metallicity<=upperZbin)
        
            if np.sum(mask_MRR)>10:
                ls_ = '-'
                ### read in Metallicity weights: ###

                weights_  = fdata['doubleCompactObjects']['weight'][...].squeeze() 

                w = weights_
                bandwidth=1

                # we want at least 10 datapoint for KDE
                # if len(w[mask_MRR])>10:
            
                # print('running kde ')
                # plot for formation channel 

                # edges = astropy.stats.bayesian_blocks(param_x[mask_MRR], fitness='events', p0=0.01)
                estimator = FFTKDE(kernel='biweight', bw=bw)    
                # print('bandwidth=', bw)

                yy_MRR = estimator.fit(param_x[mask_MRR], weights=w[mask_MRR]).evaluate(xx) 
                # yy_total = estimator.fit(param_x, weights=w).evaluate(xx) 
                rel_weight_MRR    = np.sum(w[mask_MRR]) # / (np.sum(w))
                yy_MRR *= rel_weight_MRR
                yy_min = np.minimum(yy_min, yy_MRR)
                yy_max = np.maximum(yy_max, yy_MRR)
                # axe.plot(xx[mask_kde_inside], yy_MRR[mask_kde_inside],    color=c_FC, lw=3, zorder=16, alpha=1, ls=ls_)
                axe.plot(xx, yy_MRR, color=colors_metallicity_bins[ind_Zbin], lw=4, zorder=16, alpha=1 ) #, ls=metallicity_linestyles_list[ind_Zbin])


                        
            else:
                print('%s data points for KDE, this is below the threshold 10, not drawing KDE'%len(mask_MRR))


    axe = layoutAxes(axe, nameX=nameX, nameY=nameY, setMinor=True, labelpad=0, fontsize=fs_l, labelSizeMajor=fs_major) 


    axe.set_xlim(xlim[0], xlim[1])
    axe.set_ylim(ylim[0], ylim[1])
    axe.grid(True)
    
    
    if plotYlog==True:       
        axe.set_yscale('log')



    return axe 




def set_limits(xparam, DCOtype):
    if xparam=='mass_1_LVK':
        if DCOtype=='BBH':
            bw=0.7
            ylim_threshold = 0.0004  
            xlim, ylim = [0, 48], [0.0003, 0.15]     
        elif DCOtype=='BHNS':
            bw=0.40
            ylim_threshold = 0.0004 
            xlim, ylim = [0, 25], [0.0003, .35]
        elif DCOtype=='BNS':
            bw=0.035
            ylim_threshold = 0.0004 
            xlim, ylim = [1, 3], [0.0003, 8]
    elif xparam=='mass_2_LVK':
        if DCOtype=='BBH':
            bw=0.7
            ylim_threshold = 0.0004  
            xlim, ylim = [0, 48], [0.001, 0.15]     
        elif DCOtype=='BHNS':
            bw=0.035
            ylim_threshold = 0.0004 
            xlim, ylim = [1, 3], [0.001, 8]
        elif DCOtype=='BNS':
            bw=0.035
            ylim_threshold = 0.0004 
            xlim, ylim = [1, 3], [0.001, 8]
    elif xparam=='mass_ratio_LVK':
            bw=0.035
            ylim_threshold = 0.0004 
            xlim, ylim = [0, 1], [0.001, 10]  

    elif xparam=='spin_1_LVK':
            bw=0.025
            ylim_threshold = 0.0004 
            xlim, ylim = [0, 1], [0.001, 10]  
    elif xparam=='spin_2_LVK':
            bw=0.025
            ylim_threshold = 0.0004 
            xlim, ylim = [0, 1], [0.001, 10]  

    elif xparam=='log10_t_delay':
            ylim_threshold = 0.0004 
            xlim, ylim = [-4, 2], [1, 100000]  

            if DCOtype=='BHNS':
                bw = 0.05
            elif DCOtype=='BNS':
                bw = 0.028
            elif DCOtype=='BBH':
                bw = 0.05


    return xlim, ylim,  bw, ylim_threshold


def plot_param_fc(xparam = 'mass_1_LVK', plotYlog = True, DCOtypeList = ['BNS', 'BBH', 'BHNS'], histtype='kde', Metallicity_lower=0.0, Metallicity_upper=0.03, channelList=['classic']):


    for DCOtype in DCOtypeList: 

        print('running for DCOtype = %s, for %s parameter '%(DCOtype, xparam))
        for ind_BPS, BPSmodel in enumerate(BPSnameslist[0:]):

            fig, ax = plt.subplots(1,1, figsize=(14,9))#,\
            fs_major=34
            
            print('running model %s, for %s'%(BPSmodel, DCOtype))
            

            s_text = r'$\textbf{model}$ $\textbf{%s}$: '%str(BPSmodel) + alphabetPhysicalNameDictWithEnter[BPSmodel]
            if xparam!='mass_ratio_LVK':
                ax.text(0.98, 0.97, s=s_text , rotation = 0, fontsize=24, color = 'k', va='top', ha = 'right',transform=ax.transAxes)#, weight = 'bold')
            else: 
                ax.text(0.02, 0.97, s=s_text , rotation = 0, fontsize=24, color = 'k', va='top', ha = 'left',transform=ax.transAxes)#, weight = 'bold')

            xlim, ylim,  bw, ylim_threshold = set_limits(xparam, DCOtype) 



            ax = plot_FC_distribution(axe=ax, xparam=xparam, BPSmodelName=BPSmodel, mode='MRR_PDF',\
                                          spin_threshold=0.05, bw=bw, xlim=xlim, ylim=ylim, plotYlog=plotYlog, ylim_threshold=ylim_threshold, DCOtype=DCOtype, histtype=histtype, \
                                          Metallicity_lower=Metallicity_lower, Metallicity_upper=Metallicity_upper, channelList=channelList)
            # xx = np.linspace(1000,20000,100)
            # yy = np.linspace(1000,20000,100)
            # lw=12
            # fs = 30
            # ax.plot(xx, yy,    color=colors_lighter[1], lw=lw, zorder=16, alpha=1, label=r'\textbf{MRR}')#, ls=ls_)
            # ax[2,1].plot(xx, yy, color=colors_lighter[0], lw=lw, zorder=15, alpha=1, label=r'\textbf{non-MRR}')#, ls=ls_)
            # ax[2,1].legend(fontsize=fs, frameon=False, loc='upper left')

            plt.tight_layout()  
            # plt.subplots_adjust(wspace=0.2, hspace=0.1)

            plt.savefig('./singlemodel_perZ/'+ xparam+'/super_FC_split_panel_%s_%s_%s.pdf'%(DCOtype, BPSmodel, channelList[0]), transparent=False, bbox_inches="tight",  format='pdf')
            plt.savefig('./singlemodel_perZ/'+ xparam+'/super_FC_split_panel_%s_%s_%s.png'%(DCOtype, BPSmodel, channelList[0]), transparent=False, bbox_inches="tight", dpi=600, format='png')

            print('done')
            # plt.show()  
            # plt.close()

    return 


##### PLOT LVKM1 
# plot_param_fc(xparam = 'mass_1_LVK', plotYlog = True) 

#### PLOT LVKM2 
# plot_param_fc(xparam = 'mass_2_LVK', plotYlog = True)

##### plot mass ratios 
# plot_param_fc(xparam = 'mass_ratio_LVK', plotYlog = True)

####### plot spins 
# plot_param_fc(xparam = 'spin_1_LVK', plotYlog = True, DCOtypeList=['BHNS', 'BBH'])
# plot_param_fc(xparam = 'spin_2_LVK', plotYlog = True, DCOtypeList=['BBH'])

# Metallicity_lower=0
# Metallicity_upper=0.0014


nZbins = 6
Metallicity_lower_bins, Metallicity_upper_bins = np.zeros(nZbins), np.zeros(nZbins)
metallicity_list = np.concatenate((np.asarray(metallicities_list),np.asarray([0.03])))

channelList=['classic']
for n_rough_Z_bins in range(6): 
    lowerZ, upperZ = 0 + 9*n_rough_Z_bins, 8 + 9*n_rough_Z_bins
    Metallicity_lower_bins[n_rough_Z_bins], Metallicity_upper_bins[n_rough_Z_bins]  = metallicity_list[lowerZ], metallicity_list[upperZ]

plot_param_fc(xparam = 'log10_t_delay', plotYlog = True, DCOtypeList=['BBH'], Metallicity_lower=Metallicity_lower_bins, Metallicity_upper=Metallicity_upper_bins,\
    channelList=channelList)

# plot_param_fc(xparam = 'log10_t_delay', plotYlog = True, DCOtypeList=['BNS', 'BHNS', 'BBH'], Metallicity_lower=Metallicity_lower, Metallicity_upper=Metallicity_upper)











######## TRY TO MIRROR KDE #############

# bw = 0.1
# random_data = np.random.uniform(0, 1, 10000)
# x_vals = np.linspace(0, 1, 100)

# plt.hist(random_data, bins=25, density=True)


# estimator = FFTKDE(kernel='biweight', bw=bw) 
# regular_kde = estimator.fit(random_data, weights=np.ones_like(random_data)).evaluate(x_vals) 
# # regular_kde = scipy.stats.gaussian_kde(random_data)
# plt.plot(x_vals, regular_kde, lw=3)


# kde_vals = estimator.fit(random_data, weights=np.ones_like(random_data)).evaluate(x_vals)
# unmirrored_x_vals = np.copy(x_vals)
# lower_bound=0.0
# upper_bound=1.0

# if lower_bound is not None:
#     x_vals = 2.0 * lower_bound - x_vals
#     x_vals = np.concatenate((x_vals, unmirrored_x_vals))
#     print(np.sort(x_vals))
#     kde_vals += estimator.fit(random_data, weights=np.ones_like(random_data)).evaluate(np.sort(x_vals)) 
#     x_vals = unmirrored_x_vals

# if upper_bound is not None:
#         x_vals = 2.0 * upper_bound - x_vals
#         kde_vals += estimator.fit(random_data, weights=np.ones_like(random_data)).evaluate(np.sort(x_vals)) 
#         x_vals = unmirrored_x_vals

# plt.plot(x_vals, kde_vals, lw=3)

# plt.show() 




# def mirror_KDE_values(x_vals, kde_val, lower_bound=None, upper_bound=None):

#     unmirrored_x_vals = np.copy(x_vals)

#     # evaluate the kde at the original x values
#     # kde_vals = super().evaluate(x_vals)

#     # if either bound is present then mirror the data and
#     # add the evaluated kde for the mirrored data to the original
#     if lower_bound is not None:
#         x_vals = 2.0 * lower_bound - x_vals
#         kde_vals += super().evaluate(x_vals)
#         x_vals = unmirrored_x_vals

#     if upper_bound is not None:
#         x_vals = 2.0 * upper_bound - x_vals
#         kde_vals += super().evaluate(x_vals)
#         x_vals = unmirrored_x_vals

#     return x_vals, kde_vals

















    
    
    
