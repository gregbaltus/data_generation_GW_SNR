###########################################
#Explanation of the generation of the data#
###########################################

I) Quick generation of data
_____________________________


python3 genrate_data.py --mass_min_1=20 --mass_max_1=30  --mass_min_2=28 --mass_max_2=30 --jump_mass=2 --wanted_ratio=5

This command will generate Roundup[((30-20)+1)/2]*Roundup[((30-28)+1)/2]=12 files
In fact you have to choose the minimal masses (1 and 2) and maximal masses (1 and 2) and the increment of masses (=jump_mass) for each file
It begin by generate a file with mass1=20 and mass2=28, then a file with  mass1=22 and mass2=28, then a file with  mass1=24 and mass2=28, ..., then a file with  mass1=20 and mass2=30, then a file with  mass1=22 and mass2=30, ..., then finally then a file with  mass1=30 and mass2=30.

Note you can add --nbr_initial_file=100
If you already have 100 file and you don't want to loss them. Then the first file create by the command will be wf_run0_100

Note the wanted ratio is a nimber wich will define the range of SNR of the file generate. 
For exemple --wanted_ratio=5 create file with SNR between ~[10,13]
(You need to test the value to see which range of SNR correspond to wich value or use the table here:
--wanted_ratio=3 => ~[20,26]
--wanted_ratio=3.5 => ~[16,20]
--wanted_ratio=4 => ~[13,16]
--wanted_ratio=5 => ~[10,13]
--wanted_ratio=6 => ~[9,11]
--wanted_ratio=7 => ~[6,11]
--wanted_ratio=8 => ~[4,10]
--wanted_ratio=9 => ~[4,8]
--wanted_ratio=10 => ~[4,8]
--wanted_ratio=13 => ~[4,6]


The data is store under different format : 
numpy : which is th GW, with th correspondant chirp masse and spin
.txt : where all the carateristique of the GW is write
.png :  graph of the GW (pure GW and whitened GW)


If you want to remove all the GW from the folders numpy txt png you can use the command :
./clean_data.sh



In fact genrate_data.py just call data_generation.py Roundup[((mass_max_2=30-mass_min_2)+1)/jump_mass]=2 in this exemple. 
So alternatively you can use the command
 
python3 data_generation.py --mass_min_1=20 --mass_max_1=30  --mass_min_2=28
 --mass_max_2=30 --jump_mass=2

It will do exatly the same thing. But if you try to generate to many files it will crash, because a out of memory problem. So it's recommanded to use genrate_data.py


Explanation of SNR_data.py
--------------------------
python3 SNR_data.py --mass_min_1=20 --mass_max_1=30  --mass_min_2=28 --mass_max_2=30 --jump_mass=2 --ratio_min=3 --ratio_max=4 --ratio_jump=0.5

This python scrit call genrate_data.py multiple time. The main of this code is to creat different folder contain data which each forder contain GW with a certain range of SNR (define by a ratio).
So you have to tell the minial ration (--ratio_min) maximal ratio (--ratio_max) and the increment ratio (--ratio_jump)
In our exemple SNR_data.py will create 3 folder first with ratio=3 second ratio=3.5 and last ratio=4



Explaination of data_generation.py 
-----------------------------------

parameter : --nbr_run (default 1) --mass_min_1 (default 20) --mass_max_1 (default 30) -mass_min_2 (default 20) --mass_max_2 (default 30) --jump_mass (default 1) --nbr_initial_file (default 0):

Note : --nbr_run is not relevant anymore
Note : --nbr_initial_file is important when you want to generate file without remove the other wich is present. 

This file does the data processing. It cuts the GW to the last second, genrates the Gaussian noise, injects the GW in the Gaussian noise, does the withning, normalizes.
It saved all the GW (clean and noisy) in different format located in data/pure_GW/ (data/noisy_GW/). 

Note on the SNR :  This file also calculate the SNR.
Note on the SNR and distance : here we want to obtain a certain SNR so the GW is multiply by a constante to obtain the right (range of) SNR. So the distance here have not really a true meaning.


The numpy files
----------------

Each numpy files contains the following element : [wf, chirp_mass, snr, presenceGW]
wf = th waveform
chirp_mass = chirp masse of the GW
snr = snr of the GW
presence = 1=there is a GW, 0=this is pure noise
