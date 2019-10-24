#!/usr/bin/env python3


#from glue.ligolw import ligolw
#from glue.ligolw import lsctables
#from glue.ligolw import utils as ligolw_utils

#import lal

#import lalsimulation

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

import numpy as np
import cmath
import random
import pickle
import scipy
import copy

import pycbc.noise
import pycbc.psd
import pylab
from pycbc.types import TimeSeries, zeros
from pycbc.types import complex_same_precision_as, FrequencySeries
from pycbc.filter import matched_filter
from pycbc.waveform import get_td_waveform
from pycbc.psd import interpolate, inverse_spectrum_truncation

#from lalsimulation import SimNoise
#import lal

from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt, iirdesign, zpk2tf, freqz
import scipy.io.wavfile as wavf
from optparse import OptionParser

from tqdm import tqdm
import multiprocessing as mp

'''
From injection_data :

tab_waveform[i][0] = the ith waveform (array) of type : h_plus*F_plus+h_cross*F_cross
                   = tab_h (in read_injections.py) 
tab_waveform[i][1] = the ith time (array) (useful for plot tab_waverform[0][i])
                   = tab_time (in read_injections.py) 
tab_waveform[i][2] = the waveform h_plus (array) of the ith GW
                   = tab_h_plus (in read_injections.py)
tab_waveform[i][3] = the waveform h_cross (array) of the ith GW
                   = tab_h_cross (in read_injections.py)
tab_waveform[i][4] = the angle theta of F_plus and F_cross of the ith GW
                   = tab_theta (in read_injections.py)
tab_waveform[i][5] = the angle phi of F_plus and F_cross of the ith GW
                   = tab_phi (in read_injections.py)
tab_waveform[i][6] = the mass 1 of the ith GW
                   = tab_m1 (in read_injections.py)
tab_waveform[i][7] = the mass 2 of the ith GW
                   = tab_m2 (in read_injections.py)
tab_waveform[i][8] = the componante x du spin1 of the ith GW
                   = tab_S1x (in read_injections.py)
tab_waveform[i][9] = the componante y du spin1 of the ith GW
                   = tab_S1y (in read_injections.py)
tab_waveform[i][10] = the componante z du spin1 of the ith GW
                    = tab_S1z (in read_injections.py)
tab_waveform[i][11] = the componante x du spin2 of the ith GW
                    = tab_S2x (in read_injections.py)
tab_waveform[i][12] = the componante y du spin2 of the ith GW
                    = tab_S2y (in read_injections.py)
tab_waveform[i][13] = the componante z du spin2 of the ith GW
                    = tab_S2z (in read_injections.py)
tab_waveform[i][14] = the distance of the source of the ith GW
                    = tab_distance (in read_injections.py)
tab_waveform[i][15] = inclination of the ith GW
                    = tab_inclination (in read_injections.py)
tab_waveform[i][16] = the parametre phiRef of the ith GW
                    = tab_phiRef (in read_injections.py)
tab_waveform[i][17] = the parametre deltaT (distance between each point) of the ith GW
                    = tab_deltaT (in read_injections.py)
tab_waveform[i][18] = minimal frequence of the ith GW
                    = tab_f_min (in read_injections.py)

After function RandomSecond
tab_waveform[i][19] = zero_GW (1 if there are GW 0 if not)


This description is the array as use in data_generation.py. So tab_waveform[i][j] = parametre j of the ith GW

In the folder waveform_wavenet/noisy_waveform/numpy/ files are different is a big array but define like this tab_data = [tab_wf, chirp_mass_tab]

tab_wf[i] = array represnent the ith GW
chirp_mass_tab[i] =represent the chirp masse of the ith GW

'''

#######################################################################
#                        Les arguments                                #
#######################################################################
def parse_command_line():
	parser = OptionParser(description = __doc__)
	parser.set_defaults(nbr_run=1)
	#parser.set_defaults(nbr_file=10)
	parser.set_defaults(nbr_initial_file=0)
	parser.set_defaults(mass_min_1=20)
	parser.set_defaults(mass_max_1=30)
	parser.set_defaults(mass_min_2=20)
	parser.set_defaults(mass_max_2=30)
	parser.set_defaults(jump_mass=1)
	parser.set_defaults(wanted_ratio=5)
	
	parser.add_option("--nbr_run", metavar = "1", type="int", nargs=1,  help="nbr of run, i.e. 3 run use 3 times all the GW but in different noise")
	#parser.add_option("--nbr_file", metavar = "1", type="int", nargs=1,  help="nbr of GW that you want generate. If you want more than 500 GW, it is suggest to run different time this code.")
	parser.add_option("--nbr_initial_file", metavar = "1", type="int", nargs=1,  help="nbr of GW initial. If you want more GW but conserve what you have")
	parser.add_option("--mass_min_1", metavar = "1", type="float", nargs=1,  help="mass minimal 1")
	parser.add_option("--mass_max_1", metavar = "1", type="float", nargs=1,  help="mass maximal 1")
	parser.add_option("--mass_min_2", metavar = "1", type="float", nargs=1,  help="mass minimal 1")
	parser.add_option("--mass_max_2", metavar = "1", type="float", nargs=1,  help="mass maximal 2")
	parser.add_option("--jump_mass", metavar = "1", type="float", nargs=1,  help="increment de la masse a chaque iteration")
	parser.add_option("--wanted_ratio", metavar = "1", type="float", nargs=1,  help="This ratio will define the range of SNR of the creating files")	


	args,filenames=parser.parse_args()
	
	return args,filenames


##########################################################
#         definition des F_plus et F_cross               #
##########################################################
#by definition : h = F_cross*h_cross + F_plus*h_plus
def F_cross(theta,phi):
	f_cross=(1/2)*(1+cmath.cos(theta)**2)*cmath.cos(2*phi)
	return f_cross.real

def F_plus(theta, phi):
	f_plus=cmath.cos(theta)*cmath.sin(2*phi)
	return f_plus.real


######################################################################
#                     generate GW data                               #
######################################################################
# This finction will generate wf, it will generate a big array with the GW, h_plus, h_cross,... (see description)
# This function generate GW with random parameter between the two born
# It didn t take in account the spin
# mass is on solar mass, distance is on solar mass
# remind : theta [0, 2 pi] et phi

def Generate_GW_data(distance_min, distance_max, m1_min, m1_max, m2_min, m2_max, f_lower, freq, theta_min, theta_max, phi_min, phi_max):


	tab_waveform = []
	waveform = []
	
	# /!\ les lignes en commentaires sont bonnes mais donne masse1 et masse2
	# aleatoire, distance aleatoire (le tout compris entre les bornes en arguments
	# Vu que l on desire garder un ceraine SNR on ne va pas faire ca 
	
	# /100 to have something with 2 nbr after the decimal point
	m1 = random.randint(m1_min*100, m1_max*100)/100
	m2 = random.randint(m2_min*100, m2_max*100)/100 
	dist = random.randint(distance_min, distance_max)
	
	# /!\ on decide que m1=m2 et reste aleatoire
	# ensuite on appelle une fonction qui va genere une masse pseudo aleatoire
	# etant donne la masse pour conserver un certaine SNR
	
	#m1 = random.randint(m1_min*100, m1_max*100)/100
	#m2 = m1
	#dist = Dist_goodSNR(m1) 
	
	theta = random.uniform(theta_min, theta_max)
	phi = random.uniform(phi_min, phi_max)
	# spin is not use for the GW, those parameters are  for the future
	s1x = 0
	s1y = 0
	s1z = 0
	s2x = 0
	s2y = 0
	s2z = 0
	# inclinaison, phi_ref is not use for the GW, this parameter is for the future
	inclination = 0
	phi_ref = 0
	deltaT = 1.0/freq
	

	hp, hc = get_td_waveform(approximant="SEOBNRv4",  
			mass1=m1, 
			mass2=m2,
			delta_t=1.0/freq,
			f_lower=f_lower,
			distance=dist)

	h_tot = F_cross(1.26,1.26) * hc + F_plus(1.26,1.26) * hp
	time = h_tot.sample_times
	

	h_tot = np.asarray(h_tot)
	time = np.asarray(time)
	hp = np.asarray(hp)
	hc = np.asarray(hc)
	
	# this if is just to be sure that all GW have a minimal length of 1sec
	# If not, this if will add some zero before the GW to have 1seconde
	if freq > len(h_tot):
		h_tot2 = np.zeros([freq])
		time2 = np.zeros([freq])
		hp2 = np.zeros([freq])
		hc2 = np.zeros([freq])
		j = 0
		for i in range(freq-len(h_tot), freq):
			h_tot2[i] = h_tot[j]
			time2[i] = time[j]
			hp2[i] = hp[j]
			hc2[i] = hc[j]
			j += 1
		h_tot = h_tot2
		time = time2
		hp = hp2
		hc = hc2
	
	waveform = [h_tot, time, hp, hc, theta, phi, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, dist, inclination, phi_ref, deltaT, f_lower]
	waveform = np.asarray(waveform)
	#tab_waveform = tab_waveform + [waveform]
	
	#tab_waveform = np.asarray(tab_waveform)
	
	
	return waveform



########################################################################
#                      pick one random seconde                         #
########################################################################
#/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\#
# cette fonction n est pas utiliser dans cette version du code #
#/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\#

# create a tab with nbr_sec secondes of 0 fellow by the signal and nbr_sec of 0
# pick nbr_sec secondes in this big tab (the nbr_sec secondes pick is name h_randsec)
# The nbr_sec secondes pick can contain : 1) some 0 fellow by beginning of GW
#                                         2) a part of the GW
#                                         3) the end of the GW fellow by some 0
#                                         4) only 0 
# the percent_beginning correspond to "the fraction" of 1)
# the percent_full correspond to "the fraction" of 2)
# the percent_end correspond to "the fraction" of 3)
# the percent_nothing correspond to "the fraction" of 4)
# create also an array zero_GW which is a copy of h_randsec exept ithe signal is replace by one
def RandomSecond(tab_waveform, nbr_sec, run, percent_beginning, percent_full, percent_end, percent_nothing, seed):
	print("######################")
	print("pick one random seconde")
	tab_wf_randsec = []
	random.seed(seed)
	# the if is just to be sure thatercent_beginning + percent_full + percent_end + percent_nothing is egal to 100
	if percent_beginning + percent_full + percent_end + percent_nothing != 100:
		print("percent_beginning + percent_full + percent_end + percent_nothing different from 100")
		return 0
	else:
		print("pick random seconde waveform")
		for i in tqdm(range(0, len(tab_waveform))):
			taille =  len(tab_waveform[i][0])
			choose_cat = random.randint(1,100)
			h_randsec = []
			h_plus_randsec = []
			h_cross_randsec = []
			zero_GW = [] # copy of h_randsec exept signal replace by 1 (usefull for the loss and target in deep learning algorithm)
			freq = int(1/tab_waveform[i][17])*nbr_sec
			# We choose one of the 4 categories
			
			# categorie 1) some 0 fellow by beginning of GW
			if choose_cat <= percent_beginning:
				indice_start = random.randint(0,freq)
				itr = 0
				for k in range(0, indice_start):
					h_randsec = h_randsec + [0]
					h_plus_randsec = h_plus_randsec + [0]
					h_cross_randsec = h_cross_randsec +[0]
					zero_GW = zero_GW +[0]
				for k in range(indice_start, freq):
					h_randsec = h_randsec + [tab_waveform[i][0][itr]]
					h_plus_randsec = h_plus_randsec + [tab_waveform[i][2][itr]]
					h_cross_randsec =h_cross_randsec + [tab_waveform[i][3][itr]]
					zero_GW = zero_GW +[1]
					itr = itr + 1

			# categorie 2) a part of the GW
			elif percent_beginning < choose_cat <= percent_beginning+percent_full:
				indice_start = random.randint(0, taille-freq)
				for k in range(0,freq):
					h_randsec = h_randsec + [tab_waveform[i][0][indice_start+k]]
					h_plus_randsec = h_plus_randsec + [tab_waveform[i][2][indice_start+k]]
					h_cross_randsec =h_cross_randsec + [tab_waveform[i][3][indice_start+k]]
					zero_GW = zero_GW +[1]
				
			# categorie 3) the end of the GW fellow by some 0
			elif percent_beginning+percent_full < choose_cat <= percent_beginning+percent_full+percent_end:
				indice_start = random.randint(taille-freq, taille)
				itr = 0
				for k in range(0, taille - indice_start):
					h_randsec = h_randsec + [tab_waveform[i][0][indice_start+k]]
					h_plus_randsec = h_plus_randsec + [tab_waveform[i][2][indice_start+k]]
					h_cross_randsec =h_cross_randsec + [tab_waveform[i][3][indice_start+k]]
					zero_GW = zero_GW +[1]
				
				for k in range(taille - indice_start, freq):
					h_randsec = h_randsec + [0]
					h_plus_randsec = h_plus_randsec + [0]
					h_cross_randsec = h_cross_randsec +[0]
					zero_GW = zero_GW +[0]

			# categorie 4) the only 0
			elif  percent_beginning+percent_full+percent_end < choose_cat <= percent_beginning+percent_full+percent_end+percent_nothing:
				for k in range(0,freq):
					h_randsec = h_randsec + [0]
					h_plus_randsec = h_plus_randsec + [0]
					h_cross_randsec = h_cross_randsec +[0]
					zero_GW = zero_GW +[0]
			
			time_cut = [x*tab_waveform[i][17] for x in range(0,freq)]
			h_randsec = np.asarray(h_randsec)
			time_cut=np.asarray(time_cut)
			h_plus_randsec = np.asarray(h_plus_randsec)
			h_cross_randsec = np.asarray(h_cross_randsec)
			zero_GW = np.asarray(zero_GW)
			
			tab_wf_randsec = tab_wf_randsec +[[h_randsec, time_cut, h_plus_randsec, h_cross_randsec, tab_waveform[i][4], tab_waveform[i][5], tab_waveform[i][6], tab_waveform[i][7], tab_waveform[i][8], tab_waveform[i][9], tab_waveform[i][10], tab_waveform[i][11], tab_waveform[i][12], tab_waveform[i][13], tab_waveform[i][14], tab_waveform[i][15], tab_waveform[i][16], tab_waveform[i][17], tab_waveform[i][18], zero_GW]]
		
		tab_wf_randsec = np.asarray(tab_wf_randsec)
		return tab_wf_randsec



########################################################################
#            We cut GW to just have the last second                    #
########################################################################
def Last_second(tab_waveform, nbr_sec, run):
        print("#######################")
        print("cut last second")
        tab_wflastsec=[]
        for i in tqdm(range(0, len(tab_waveform))):
	       	taille=len(tab_waveform[i][0])
	       	freq=int(1/tab_waveform[i][17])*nbr_sec
	       	h_cut=[]
	       	time_cut=[x*tab_waveform[i][17] for x in range(0,freq)]
	       	h_plus_cut=[]
	       	h_cross_cut=[]
	       	#The if is just in case the waveform lasts less than 1 sec 
	       	if freq>taille:
	       	        for k in range(0,freq-taille):
	       	                step=freq-k
	       	                h_cut=h_cut+[0]
	       	                h_plus_cut=h_plus_cut+[0]
	       	                h_cross_cut=h_cross_cut+[0]
	       	        for k in range(freq-taille, freq):
	       	                step=freq-k
	       	                h_cut=h_cut+[tab_waveform[i][0][taille-step]]
	       	                h_plus_cut=h_plus_cut+[tab_waveform[i][2][taille-step]]
	       	                h_cross_cut=h_cross_cut+[tab_waveform[i][3][taille-step]]
	       	else:
	       	        for k in range(0,freq):
	       	                step=freq-k
	       	                h_cut=h_cut+[tab_waveform[i][0][taille-step]]
	       	                h_plus_cut=h_plus_cut+[tab_waveform[i][2][taille-step]]
	       	                h_cross_cut=h_cross_cut+[tab_waveform[i][3][taille-step]]
	       	h_cut=np.asarray(h_cut)
	       	time_cut=np.asarray(time_cut)
	       	h_plus_cut=np.asarray(h_plus_cut)
	       	h_cross_cut=np.asarray(h_cross_cut)
	       	tab_wflastsec=tab_wflastsec+[[h_cut, time_cut, h_plus_cut, h_cross_cut, tab_waveform[i][4], tab_waveform[i][5], tab_waveform[i][6], tab_waveform[i][7], tab_waveform[i][8], tab_waveform[i][9], tab_waveform[i][10], tab_waveform[i][11], tab_waveform[i][12], tab_waveform[i][13], tab_waveform[i][14], tab_waveform[i][15], tab_waveform[i][16], tab_waveform[i][17], tab_waveform[i][18]]]
		
        tab_wflastsec=np.asarray(tab_wflastsec)
        return tab_wflastsec



######################################################################
#                    data for WaveNet and normalisation              #
######################################################################
def Data_wavenet(tab_waveform, tab_snr, presenceGW, chemin, nbr_initial_file,  run, fs):
	# this fonction save the waveform under different file
	# first file is a numpy file with the GW (or whitend GW) and also the chirp mass, and the SNR
	# second file is a .txt and have all the information about the GW
	# third file is a .png and is a graph of the GW
	# this file also normalize the GW (dividing by the biggest element) 
	# difference between tab_snr and tab_SNRcomplet:
	# 	In this file we select only a short interval from the entire GW
	#	tab_snr represent the SNR on this short interval
	#	tab_SNRcomplet represent the SNR on the whole GW (in fact not really the whole GW but the last second of this GW)
	print("#######################")
	print("data wavenet")
	m_solaire=1.9884*10**30
        
	#we normalize by dividing all the elements by the biggest element
	tab_wf=[]
	for i in range (0,len(tab_waveform)):
		if tab_waveform[i][0].max(0) == 0:
			normalisation = 1
		else:
			normalisation=tab_waveform[i][0][np.argmax(tab_waveform[i][0])]
		wf_normalisation=tab_waveform[i][0]/normalisation
		tab_wf=tab_wf+[wf_normalisation]
	tab_wf=np.asarray(tab_wf)


	chirp_mass_tab=[]
	#we save our data in different type of file
	for i in tqdm(range (0,len(tab_waveform))):
		
		####the file numpy####
		m_1 = tab_waveform[i][6] #/ m_solaire (already in solar mass)
		m_2 = tab_waveform[i][7] #/ m_solaire (already in solar mass)
		
		# we create chirp masse here, but it is save in the .txt file
		chirp_mass = ((m_1*m_2)**(3./5.)) / ((m_1 + m_2)**(1./5.))
		chirp_mass_tab = chirp_mass_tab + [chirp_mass]
		
		tab_data = [tab_wf[i], chirp_mass_tab[i], tab_snr[i], presenceGW[i]]
		nbr_name = i + nbr_initial_file
		with open("%s" %chemin + "numpy/wf_run%s" %run + "_%s" %nbr_name, 'wb') as fichier_numpy:
			mon_pickler = pickle.Pickler(fichier_numpy)
			mon_pickler.dump(tab_data)

		####the file txt###
		nbr_name = i + nbr_initial_file
		with open("%s" %chemin + "txt/wf_run%s" %run+ "_%s.txt" %nbr_name, 'w') as fichier_txt:
			freq=int(1/tab_waveform[i][17])	
			fichier_txt.write("waveform %s" %nbr_name + "\n" + "\n")
			fichier_txt.write("angle theta =  %s" %tab_waveform[i][4] + "\n")
			fichier_txt.write("angle phi =  %s" %tab_waveform[i][5] + "\n")
			fichier_txt.write("masse 1 =  %s" %m_1 + " M_solaire" + "\n")
			fichier_txt.write("masse 2 =  %s" %m_2 + " M_solaire" + "\n")
			fichier_txt.write("spin 1 composante x =  %s" %tab_waveform[i][8] + "\n")
			fichier_txt.write("spin 1 composante y =  %s" %tab_waveform[i][9] + "\n")
			fichier_txt.write("spin 1 composante z =  %s" %tab_waveform[i][10] + "\n")
			fichier_txt.write("spin 2 composante x =  %s" %tab_waveform[i][11] + "\n")
			fichier_txt.write("spin 2 composante y =  %s" %tab_waveform[i][12] + "\n")
			fichier_txt.write("spin 2 composante z =  %s" %tab_waveform[i][13] + "\n")
			fichier_txt.write("distance de la source =  %s" %tab_waveform[i][14] + "\n")
			fichier_txt.write("inclinaison du plan de la source =  %s" %tab_waveform[i][15] + "\n")
			fichier_txt.write("phi de reference =  %s" %tab_waveform[i][16] + "\n")
			fichier_txt.write("frequence (nbr de point par seconde) =  %s" %freq + "\n")
			fichier_txt.write("frequence minimum =  %s" %tab_waveform[i][18] + "\n")
			fichier_txt.write("chirp_mass =  %s" %chirp_mass_tab[i] + "\n")
			fichier_txt.write("the SNR is = %s"%tab_snr[i] + "\n")
			if presenceGW[i] == 0:
				fichier_txt.write("There is no GW, only Gaussian noise" + "\n")
			else:
				fichier_txt.write("There is a GW" + "\n")
		####the file png####
		tab_time=[]
		deltaT=1. / fs
		tab_time=[x*deltaT for x in range(len(tab_wf[i]))]
		tab=np.asarray(tab_time)
		nbr_name = i + nbr_initial_file
		plt.figure('%s' %chemin + 'png/wf_run%s'%run + '_%s.png'%nbr_name)
		plt.plot(tab_time, tab_wf[i], 'g-')
		plt.savefig('%s' %chemin + 'png/wf_run%s'%run + '_%s.png'%nbr_name)
                
	#the README#
	#with open("%s" %chemin_clean + "README_xml%s.txt"%xml, 'w') as README:
	#	README.write("The folder numpy/txt/png contains %s" %len(tab_waveform) + " files wf for the  xml%s"%xml + " files")







######################################################################
#                   FAKE GAUSSIAN NOISE                              #
######################################################################
def fake_gaussian_noise(nbr_sec, seed_local, fs, flow):
        # The color of the noise matches a PSD which you provide
        delta_f= 1.0 / 16
        flen = int (fs / delta_f) + 1
        psd = pycbc.psd.aLIGOZeroDetHighPower(flen, delta_f, flow)

        frequence=np.empty([len(psd)])
        psd_numpy=np.empty([len(psd)])
        for i in range(0, len(psd)):
                frequence[i]=i
                psd_numpy[i]=psd[i]

        # Generate 1 seconds of noise at fs Hz
        freq=fs
        delta_t = 1.0 / freq
        tsamples = int(nbr_sec / delta_t)
        ts = pycbc.noise.noise_from_psd(tsamples, delta_t, psd, seed_local)

        #we transform into np.array
        ts_numpy=np.empty([len(ts)])
        for i in range(0, len(ts_numpy)):
                ts_numpy[i]=ts[i]
        return ts_numpy



################################################################
#                       Whitening                              #
################################################################
def whitening(strain, interp_psd, dt):
    Nt = len(strain)
    freqs = np.fft.rfftfreq(Nt, dt)

    # whitening: transform to freq domain, divide by asd, then transform back, 
    # taking care to get normalization right.
    hf = np.fft.rfft(strain)
    white_hf = hf / (np.sqrt(interp_psd(freqs) /dt/2.))
    white_ht = np.fft.irfft(white_hf, n=Nt)
    return white_ht


######################################################################
#                FAKE GW WITH GAUSSIAN NOISE                         #
######################################################################
def fake_noisyGW(tab_waveform, seed, run, fs, wanted_ratio):
	nbr_sec = 15
	nbr_sec_begin = 7
	nbr_sec_end = 8
	NFFT = 1*fs
	fmin = 15
	fmax = 200000
	print("########################")
	print("creation of noisy signal")
	tab_whiten=[]
	tab_snr = []
	random.seed(seed)
	for i in tqdm(range(0,len(tab_waveform))):
		#creation of the fake noise
		seed_local=random.randint(1,1000)
		fake_noise = fake_gaussian_noise(nbr_sec, seed_local, fs, fmin)
		# injection of GW into the noise
		noisyGW = np.empty([len(fake_noise)])
		# creation of the template (for SNR)
		template = np.zeros([len(fake_noise)]) # idem pure_GW but with GW in the end(for SNR)
		# creation of pureGW (no use exept plot)
		pureGW = np.zeros([len(fake_noise)]) #pure GW in nbr_sec of 0 
		
		# Research of the maximum
		max_fake_noise = np.amax(fake_noise)
		max_waveform =  np.amax(tab_waveform[i][0])
		if max_waveform > 0 : 
			#wanted_ratio = 55 #it s an argument now
			ajustement = max_fake_noise/(max_waveform*wanted_ratio)
		else: 
			ajustement = 1
		
		# injection of GW into noise
		k=0
		for j in range(0, len(fake_noise)):
			if j >= fs*nbr_sec_begin and j < fs*nbr_sec_end:
				noisyGW[j] = fake_noise[j] + tab_waveform[i][0][k]*ajustement
				pureGW[j] = tab_waveform[i][0][k]
				k=k+1
			else:
				noisyGW[j] = fake_noise[j]
				pureGW[j] = 0
		# creation of the template
		k = 0
		for j in range(0, len(fake_noise)):
			if j < len(fake_noise) - len(tab_waveform[i][0]):
				template[j] = 0
			else:
				template[j] = tab_waveform[i][0][k]
				k += 1
		
		########## the calcule of the SNR ##############
		# First we have to calculate the SNR (the following step is neccessary)
		# we calculate the SNR of noisy GW
		noisyGW_timeseries = TimeSeries(noisyGW,  delta_t=tab_waveform[i][17])
		# we need to calculate again the psd 
		# we use just 4 sec of data to calulate psd (weltch method)
		psd  = noisyGW_timeseries.psd(4) 
		# we need to interpolate our psd on all data
		psd = interpolate(psd, noisyGW_timeseries.delta_f)
		# we need to inform that the noise was creat with a lowfrequency cutoff
		psd = inverse_spectrum_truncation(psd, 4 * noisyGW_timeseries.sample_rate, low_frequency_cutoff=fmin)
		
		#Second we calcule the SNR 
		# we calculate the SNR 
		template_timeseries = TimeSeries(template, delta_t=tab_waveform[i][17])
		noisyGW_timeseries =  TimeSeries(noisyGW, delta_t=tab_waveform[i][17])
		if template.max(0) == 0:
			snrp = 0
			tab_snr = tab_snr + [snrp]
		else:
			snr = matched_filter(template_timeseries, noisyGW_timeseries, psd=psd, low_frequency_cutoff=fmin)
			#  4+1 sec remove at the begining of SNR and 4 seconde remove at the end
			# +1 is because the length of the template is one sec
			snr = snr.crop(4 + 1, 4)
			peak = abs(snr).numpy().argmax()
			snrp = snr[peak]
			tab_snr = tab_snr + [abs(snrp)]

		
		# on fait le withening
		# I) methode pycbc
		# data and then flattening the frequency response.
    		# (1) The first option sets the duration in seconds of each
    		#     sample of the data used as part of the PSD estimate.
    		# (2) The second option sets the duration of the filter to apply
		#whitened = noisyGW_timeseries.whiten(15, 4)
		# on fait aussi un bandpass
                #whitened_bp = whitened.highpass_fir(30, 512).lowpass_fir(250, 512)
		# le whitening pycbc donne meme chose mais modifie la longueur
		# par facilite je reprends l autre whitening
		
		# II) method numpy
		Pxx, freqs = mlab.psd(noisyGW, Fs = freq, NFFT=NFFT, 
                      noverlap=NFFT/2, window=np.blackman(NFFT))
		# We will use interpolations of the ASDs computed above for whitening:
		psd_numpy = interp1d(freqs, Pxx)
		whitened = whitening(noisyGW, psd_numpy, tab_waveform[i][17])
		# We need to suppress the high frequencies with some bandpassing:
		high_freq = 600.
		low_freq  = fmin
		bb, ab = butter(4, [low_freq*2./fs, high_freq*2./fs], btype='band')
		whitened_bp = filtfilt(bb, ab, whitened)
		
		# Select the secon containing the GW
		withen_bp_onesec = np.empty([fs])
		for j in range(0, len(withen_bp_onesec)):
			withen_bp_onesec[j] = whitened_bp[nbr_sec_begin*fs + j]
		tab_whiten = tab_whiten + [withen_bp_onesec]
		'''
		######################################################################
		### PLOT this is no necessary is just to be sure everithing is OK ###
		# just some  plot to see if everything looks fine
		temps = np.empty([len(fake_noise)])
		for tps in range(0, len(fake_noise)):
			temps[tps] = tps

		temps_onesec = np.empty([len(tab_waveform[i][0])])
		for tps in range(0, len(tab_waveform[i][0])):
			temps_onesec[tps] = tps
		
		# noisy GW
		plt.figure("noisyGW")
		plt.plot(temps, noisyGW, 'b')
		plt.plot(temps, pureGW, 'r')
		plt.show()
		# template 
		plt.figure("template")
		plt.plot(temps, template, 'r')
		plt.show()
		# the GW
		plt.figure("template")
		plt.plot(temps_onesec, tab_waveform[i][0], 'r')
		plt.show()	
		# psd 
		plt.figure("psd")
		plt.loglog(psd.sample_frequencies, psd)
		pylab.xlim(10, 2048)
		pylab.show()
		# SNR
		plt.figure("SNR")
		pylab.plot(snr.sample_times, abs(snr))
		pylab.grid()
		pylab.show()
		# Whitening
		whitened_bp = TimeSeries(whitened_bp, delta_t=tab_waveform[i][17])
		plt.figure("whitening")
		pylab.plot(whitened_bp.sample_times, whitened_bp)
		#pylab.xlim(7,8)
		pylab.grid()
		pylab.show()
		#whitening onesec
		whiten_bp_onesec_timeseries = TimeSeries(withen_bp_onesec, delta_t=tab_waveform[i][17])
		plt.figure("whitening_one_sec")
		pylab.plot(whiten_bp_onesec_timeseries.sample_times, whiten_bp_onesec_timeseries)
		pylab.xlim(0.6, 1.0)
		pylab.grid()
		pylab.show()
		###################################################################
		'''
	tab_snr = np.asarray(tab_snr)
	tab_whiten = np.asarray(tab_whiten)
	# we create the big array (with all the parameters,...)
	tab_waveformnoisy=copy.deepcopy(tab_waveform)
	for i in range(0,len(tab_waveformnoisy)):
		tab_waveformnoisy[i][0]=tab_whiten[i]

	return tab_waveformnoisy, tab_snr



####################################################################
#           distance in respct of the mass for good SNR            #
####################################################################
def Dist_goodSNR(m):
	# this function is rudly made, and serve to conserve a certain SNR
	# to do that we give the mass and this function give a distance
	# the SNR is between 5 and 18
	# work only if m1=m2

	distance = 0
	
	if m < 15:
		distance = 400
	elif m < 35 and m > 15 :
		#distance = random.randint(400, 1000)
		distance = 400
	elif m >35 : 
		#distance = random.randint(1000, 2000)
		distance = 400
	return distance

######################################################################
#                 Randomly no GW                                     #
######################################################################
# This function will make a wf = 0 with a certain probability
def KillGW(tab_waveform, seed, run, proba):
	tab_wf = []
	presenceGW = []
	print("####################")
	print("Proba to kill the GW")
	for i in tqdm(range(0, len(tab_waveform))):
		new_waveform = np.empty([len(tab_waveform[i][0])])
		new_hplus = np.empty([len(tab_waveform[i][2])])
		new_hcross = np.empty([len(tab_waveform[i][3])])
		choix = np.random.rand(1)
		if choix[0] <= proba :
			new_waveform = np.zeros([len(tab_waveform[i][0])])
			new_hplus = np.zeros([len(tab_waveform[i][2])])
			new_hcross = np.zeros([len(tab_waveform[i][3])])
			presenceGW = presenceGW + [0]
		else : 
			new_waveform = copy.deepcopy(tab_waveform[i][0])
			new_hplus = copy.deepcopy(tab_waveform[i][2])
			new_hcross = copy.deepcopy(tab_waveform[i][3])
			presenceGW = presenceGW + [1]

		tab_wf = tab_wf + [[new_waveform, tab_waveform[i][1],tab_waveform[i][2], tab_waveform[i][3],tab_waveform[i][4], tab_waveform[i][5], tab_waveform[i][6], tab_waveform[i][7], tab_waveform[i][8], tab_waveform[i][9], tab_waveform[i][10], tab_waveform[i][11], tab_waveform[i][12], tab_waveform[i][13], tab_waveform[i][14], tab_waveform[i][15], tab_waveform[i][16], tab_waveform[i][17], tab_waveform[i][18]]]

	tab_wf = np.asarray(tab_wf)
	presenceGW = np.asarray(presenceGW)
	
	return tab_wf, presenceGW

#####################################################################
#                           Main                                    #
#####################################################################

args,filenames = parse_command_line()

#the maximal number of file by xml files is 300 (read_injections.py creates 300 GW by xml files)
#nbr_file = args.nbr_file
#print(" ")
#print("number of files is ", nbr_file)
nbr_initial_file = args.nbr_initial_file
print("number of the first file is ", nbr_initial_file)
nbr_run= args.nbr_run #number of run by xml
print("nbr of run by file is ", nbr_run)
# number of run is there if you want to use the same GW multiple time with different noise

# parameter for data GW
distance_min = 500
distance_max = 1500
m1_min = args.mass_min_1
m1_max = args.mass_max_1
m2_min = args.mass_min_2
m2_max = args.mass_max_2
f_lower = 25
freq = 4096
theta_min = 1.26
theta_max = 1.26
phi_min = 1.26
phi_max = 1.26
jump_mass = args.jump_mass

wanted_ratio = args.wanted_ratio

tab_m1 = []
tab_m2= []
m1=m1_min
while m1 <= m1_max:
	tab_m1 = tab_m1 + [m1]
	m1 = m1 + jump_mass
m2=m2_min
while m2 <= m2_max:
        tab_m2 = tab_m2 + [m2]
        m2 = m2 + jump_mass

for run in range(0,nbr_run):
	print(" ")
	print("########RUN",run,"########")
	if run==0:
		# creation of the waveform
		print("############### data generation #####################")
		tab_waveform= []
		for j in tqdm(tab_m2):
			# Step 1: Init multiprocessing.Pool()
			pool = mp.Pool(mp.cpu_count())
			tab_waveform_m1 =  pool.starmap(Generate_GW_data, [(distance_min, 
						distance_max, 
						i, 
						i, 
						j, 
						j, 
						f_lower, 
						freq, 
						theta_min, 
						theta_max, 
						phi_min, 
						phi_max) for i in tqdm(tab_m1)])
			# Step 3: Don't forget to close
			pool.close()
			tab_waveform = tab_waveform + tab_waveform_m1
		tab_waveform = np.asarray(tab_waveform)
		
		'''	
		#plt.figure("apres_creation")
		#plt.plot(tab_waveform[0][1], tab_waveform[0][0], 'g-')
		#plt.show()
		# cut to last second
		nbr_sec = 1
		#tab_lastsec = Last_second(tab_waveform, nbr_sec, run)
		percent_beginning = 25
		percent_full = 25
		percent_end = 25
		percent_nothing = 25
		seed = nbr_initial_file*10 + run
		tab_lastsec = RandomSecond(tab_waveform, nbr_sec, run, percent_beginning, percent_full, percent_end, percent_nothing, seed)
		# This cut is to calculate the SNR of the whole end of the GW (see description of Data_wavenet) 
		'''
		print (" ")
		print("#############")
		print("#Cut last sec#")
		print("#############")
		nbr_sec = 1
		tab_lastsec =Last_second(tab_waveform, nbr_sec, run)
		#pool = mp.Pool(mp.cpu_count())
		#tab_lastsec = pool.starmap(Last_second, [(tab_waveform, nbr_sec, run) for i in tqdm(range(0, len(tab_waveform)))])
		#pool.close()
	else:
		print("NOT RUN 0, we only change noise")
	
	# We kill a GW on two to have sometimes pure Gaussian noise
	seed = run + nbr_run * nbr_initial_file + 100*wanted_ratio # to be sure is different each time
	proba = 1.
	tab_killGW, presenceGW = KillGW(tab_lastsec, seed, run, proba)
	
	# save file under different file (only GW)
        # snr_empty is a array full of 0, just to significate that is the pure wf without SNR
	snr_empty = np.zeros([len(tab_lastsec)])
	Data_wavenet(tab_killGW, snr_empty, presenceGW ,"data/pure_GW/", nbr_initial_file,run, freq)	

	# inject into noise 
	seed = run + nbr_run * nbr_initial_file + 100*wanted_ratio # to be sure is different each time
	tab_noisyGW, tab_snr = fake_noisyGW(tab_killGW, seed, run, freq, wanted_ratio)
	
	# save file under different file (whitned GW)
	Data_wavenet(tab_noisyGW, tab_snr, presenceGW ,"data/noisy_GW/", nbr_initial_file, run, freq) 
	
