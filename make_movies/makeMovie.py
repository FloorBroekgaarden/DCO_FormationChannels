import sys

sys.path.append('/usr/local/lib/python2.7/site-packages')

import cv2
import os
import string
import os
import moviepy.video.io.ImageSequenceClip
from PIL import Image
import glob

Image.MAX_IMAGE_PIXELS = 933120000



def makeMovie_merger_rates(DCOtype='BBH', fps=.4, duration=300):
	'''
	whichRate = 'intrinsic' or 'observed'
	fps=0.4, frames per second
	duration = duration of the movie 
	'''




	GSMFs = ['Panter et al. (2004) Single', 'Furlong et al. (2015) Single', 'Furlong et al. (2015) Double']
	MZs   = [ 'Langer et al. (2006)'      , 'Langer et al. +offset (2006)', 'Ma et al. (2015)']
	SFRs  = ['Madau et al. (2014)'         ,'Strolger et al. (2004)',     'Madau et al. (2017)']


	MSSFRnameslist = []
	MSSFRnameslist.append('000') # add phenomenological 



	for ind_SFR, SFR in enumerate(SFRs):
		ind_x = ind_SFR + 1
		for ind_GSMF, GSMF in enumerate(GSMFs):
			ind_y = ind_GSMF + 1
			for ind_MZ, MZ in enumerate(MZs):
				ind_z = ind_MZ + 1
				
				MSSFRnameslist.append('%s%s%s'%(ind_x, ind_y, ind_z))



                
  
	image_folder = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/plottingCode/paper-Fig_X_total-redshift-rates/zrate_bigpanel_per_mssfr/'

	images = []
	# print('get here')
	for ind_m, SFRD_model in enumerate(MSSFRnameslist):
		images.append(image_folder +   'RatesZ_' + DCOtype + '_' + SFRD_model + '.png')



	image_files = images
	# clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
	# clip.write_videofile(image_folder+'movie_'+ 'redshift_rates_' + DCOtype + '.mp4')

	# make also gif:
 
	# Create the frames
	frames = []
	# imgs = glob.glob("*.png")
	for i in images:
	    new_frame = Image.open(i)
	    frames.append(new_frame)
	 
	# Save into a GIF file that loops forever
	frames[0].save(image_folder+'gif_'+ 'redshiftrates_' + DCOtype +  '.gif', format='GIF',
	               append_images=frames[1:],
	               save_all=True,
	               duration=duration, loop=0)


	print('done')
	return 



#####



import sys

sys.path.append('/usr/local/lib/python2.7/site-packages')
# import cv2
import os
import string
import os
import moviepy.video.io.ImageSequenceClip
from PIL import Image
import glob

Image.MAX_IMAGE_PIXELS = 933120000

def obtain_movie_image_dict(path_to_folder_figures='path_to_figures_folder', plot_param_name='qZAMS_separationInitial_logx', loop_modus='loop_by_model', adjustedChannelList=['classic'], BPSnameslist=['A']):
    
    image_folder = path_to_folder_figures
    images = []
    if loop_modus=='loop_by_model':
        for ind_bps, BPSmodelName in enumerate(BPSnameslist[0:]): # loop over population synthesis models
            for nrC, _ in enumerate(adjustedChannelList):
                images.append(image_folder +   'initial2D_' + BPSmodelName  + '_'+ str(nrC) + '.png')
                
    elif loop_modus=='loop_by_formation_channel':
        for nrC, _ in enumerate(adjustedChannelList):
            for ind_bps, BPSmodelName in enumerate(BPSnameslist[0:]): # loop over population synthesis models
                images.append(image_folder +   'initial2D_' + BPSmodelName  + '_'+ str(nrC) + '.png')
    
    movie_name = plot_param_name + '_' + loop_modus 
    
    return images, movie_name 
    

def makeMovie(fps=.4, duration=300, image_files=[], path_to_movie_name='path_to_movie'):
    '''
    fps=0.4, frames per second
    duration = duration of the movie 
    '''

    clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
    clip.write_videofile(path_to_movie_name + '.mp4')

    # make also gif:
    # Create the frames
    # frames = []
    # # imgs = glob.glob("*.png")
    # for i in images:
    #     new_frame = Image.open(i)
    #     frames.append(new_frame)

    # # Save into a GIF file that loops forever
    # frames[0].save(path_to_movie_name+  '.gif', format='GIF',append_images=frames[1:],save_all=True,duration=duration, loop=0)


    print('done')
    
    return 



sys.path.append('../../common_code')
from PostProcessingScripts import * 


# plot_param_name = 'qZAMS_separationInitial_logx'
# xparam_key, xparam_function = None,'qZAMS'
# yparam_key, yparam_function = 'separationInitial', None
 


plot_param_name = 'qLVK_tdelay_logy'
xparam_key, xparam_function = None, 'qLVK'
yparam_key, yparam_function = None, 'tdelay'

fps = 4
duration=60

## GENERAL 
pathData='/Volumes/SimonsFoundation/DataDCO/'
adjustedChannelList= ['classic', 'stable B no CEE', 'immediate CE',  r'double-core CE', 'other', 'vii']
base_path = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/plottingCode/Figure_plot_2D_properties_fc_formation_efficiency'
path_to_folder_figures = base_path + '/figures/' + plot_param_name +  '/' 
###
for loop_modus in [ 'loop_by_formation_channel', 'loop_by_model']:
	image_files, movie_name = obtain_movie_image_dict(path_to_folder_figures, plot_param_name, loop_modus, adjustedChannelList, BPSnameslist)
	path_to_movie_name = path_to_folder_figures + movie_name
	makeMovie(fps, duration, image_files, path_to_movie_name)








makeMovie_redshiftRates=False



# Run rhis using python 3!! 

if makeMovie_redshiftRates==True:
	fps = 0.4
	duration=60
	for DCOtype in ['NSNS' ]: # , 'BHNS']:
		print('at DCOtype ', DCOtype)
		makeMovie_merger_rates(DCOtype=DCOtype, fps=fps, duration=duration)




