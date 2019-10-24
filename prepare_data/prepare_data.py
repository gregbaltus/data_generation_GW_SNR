#!/usr/bin/env python3


from optparse import OptionParser
import math
import subprocess


#######################################################################
#                        Les arguments                                #
#######################################################################
def parse_command_line():
        parser = OptionParser(description = __doc__)
        parser.set_defaults(debut=3.0)
        parser.set_defaults(fin=4.0)
        parser.set_defaults(jump=0.5)
        parser.set_defaults(trainset=0) #0 mean no
        parser.set_defaults(testset=0) #0 mean no

        parser.add_option("--debut", metavar = "1", type="float", nargs=1,  help="1er ratio du nom des dossiers")
        parser.add_option("--fin", metavar = "1", type="float", nargs=1,  help="dernier ratio du nom des dossiers")
        parser.add_option("--jump", metavar = "1", type="float", nargs=1,  help="increment du nom des dossiers")
        parser.add_option("--trainset", metavar = "1", type="float", nargs=1,  help="if you want that create into the folder a folder trainset with the data in it --trainset=1 else --trainset=0")
        parser.add_option("--testset", metavar = "1", type="float", nargs=1,  help="if you want that create into the folder a folder testset with the data in it --testset=1 else --testset=0")
        args,filenames=parser.parse_args()
        
        return args,filenames


args,filenames = parse_command_line()

debut = args.debut
fin = args.fin
jump = args.jump
trainset = args.trainset
testset = args.testset

increment=debut
while increment <= fin:
    subprocess.run(["rm","-r", "data_ratio%s"%increment])
    if trainset == 1:
        subprocess.run(["mkdir", "data_ratio%s"%increment])
        subprocess.run(["cp", "-r", "../data_ratio%s/noisy_GW/numpy/"%increment, "data_ratio%s/trainset/"%increment])
        # petite boucle pour renomer les fichier quand on a généré le bruit
        #for i in range(0,2116):
            #new_name = i+2116
            #subprocess.run(["mv", "data_ratio%s/trainset/"%increment + "wf_run0_%s"%i, "data_ratio%s/trainset/"%increment +"wf_run0_%s"%new_name])
    elif testset == 1:
        subprocess.run(["mkdir", "data_ratio%s"%increment])
        subprocess.run(["cp", "-r", "../data_ratio%s/noisy_GW/numpy/"%increment, "data_ratio%s/testset/"%increment])
    else:
        subprocess.run(["cp", "-r", "../data_ratio%s/noisy_GW/numpy/"%increment, "data_ratio%s/"%increment])
    increment+=jump
