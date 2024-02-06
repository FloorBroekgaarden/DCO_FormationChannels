import sys

sys.path.append('/usr/local/lib/python2.7/site-packages')

import cv2
import os
import string
import os
import moviepy.video.io.ImageSequenceClip
from PIL import Image
import glob
import numpy as np 





def makeMovie_masses(DCOtype='BHBH', fps=.4, duration=300):
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



                
	BPSnameslist = ['A','B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',  'P', 'Q', 'R', 'S', 'T'] #list(string.ascii_uppercase)[0:nModels]
	image_folder = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/plottingCode/paper-Figure_X_FC_Properties/singlemodel/'

	image_files = []

	for ind_m, BPSmodel in enumerate(BPSnameslist):
		image_files.append(image_folder +   'super_FC_split_panel_'+DCOtype + '_' + BPSmodel +'.png' )


	clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
	clip.write_videofile(image_folder+'movie_'+ 'models_' + DCOtype + '.mp4')


	# make also gif:
 
	# Create the frames
	frames = []
	# imgs = glob.glob("*.png")
	for i in image_files:
	    new_frame = Image.open(i)
	    frames.append(new_frame)
	 
	# Save into a GIF file that loops forever
	frames[0].save(image_folder+'gif_'+ 'models_' + DCOtype +  '.gif', format='GIF',
	               append_images=frames[1:],
	               save_all=True,
	               duration=duration, loop=0)


	print('done')
	return 




def makeMovie_delaytimes(DCOtype='BHBH', fps=.4, duration=300):
	'''
	whichRate = 'intrinsic' or 'observed'
	fps=0.4, frames per second
	duration = duration of the movie 
	'''

                
	BPSnameslist = ['A','B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',  'P', 'Q', 'R', 'S', 'T'] #list(string.ascii_uppercase)[0:nModels]
	image_folder = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/plottingCode/paper-Figure_X_FC_Properties/singlemodel/'

	image_files = []

	for ind_m, BPSmodel in enumerate(BPSnameslist):
		for redshift in [0.19230769230769232, 2.1153846153846154, 3.269230769230769, 6.346153846153847, 7.5]:
			redshift = np.round(redshift,4)
			image_files.append(image_folder +   'log10_t_delay/fcplot_'+DCOtype + '_' + BPSmodel + '_z' + str(redshift) +  '_wformation.png' )

	clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
	clip.write_videofile(image_folder+'movie_'+ 'models_' + DCOtype + '_log10delay_wform.mp4')


	# make also gif:
 
	# Create the frames
	frames = []
	# imgs = glob.glob("*.png")
	for i in image_files:
	    new_frame = Image.open(i)
	    frames.append(new_frame)
	 
	# Save into a GIF file that loops forever
	frames[0].save(image_folder+'gif_'+ 'models_' + DCOtype +  '_log10delay_wform.gif', format='GIF',
	               append_images=frames[1:],
	               save_all=True,
	               duration=duration, loop=0)


	print('done')
	return 




makeMovie_mass=False
makeMovie_tdelay = True 
# makeMovie_intrinsicM1=False
# makeMovie_Mtotq_FormChannels=False
# makeMovieChiBH=False









# Run rhis using python 3!! 

if makeMovie_mass==True:
	makeMovie_masses(DCOtype='BBH', fps=.4, duration=300)
	makeMovie_masses(DCOtype='BHNS', fps=.4, duration=300)
	makeMovie_masses(DCOtype='BNS', fps=.4, duration=300)

if makeMovie_tdelay==True:
	# makeMovie_delaytimes(DCOtype='BBH', fps=.4, duration=300)
	# makeMovie_delaytimes(DCOtype='BHNS', fps=.4, duration=300)
	makeMovie_delaytimes(DCOtype='BNS', fps=.4, duration=300)

# if makeMovie_intrinsicRates==True:
# 	makeMovie_rates(whichRate='intrinsic')


# if makeMovie_intrinsicM1==True:
# 	makeMovie_MBH1(whichRate='intrinsic')

# if makeMovie_Mtotq_FormChannels==True: 
# 	makeMovie_Mtotq_FormChannel(whichRate='intrinsic')

