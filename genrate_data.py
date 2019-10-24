#!/usr/bin/env python3


from optparse import OptionParser
import math
import subprocess


#######################################################################
#                        Les arguments                                #
#######################################################################
def parse_command_line():
	parser = OptionParser(description = __doc__)
	parser.set_defaults(mass_min_1=20)
	parser.set_defaults(mass_max_1=30)
	parser.set_defaults(mass_min_2=20)
	parser.set_defaults(mass_max_2=30)
	parser.set_defaults(jump_mass=1)
	parser.set_defaults(nbr_initial_file=0)
	parser.set_defaults(wanted_ratio=5)

	parser.add_option("--nbr_initial_file", metavar = "1", type="int", nargs=1,  help="nbr of GW initial. If you want more GW but conserve what you have")
	parser.add_option("--mass_min_1", metavar = "1", type="float", nargs=1,  help="mass minimal 1")
	parser.add_option("--mass_max_1", metavar = "1", type="float", nargs=1,  help="mass maximal 1")
	parser.add_option("--mass_min_2", metavar = "1", type="float", nargs=1,  help="mass minimal 1")
	parser.add_option("--mass_max_2", metavar = "1", type="float", nargs=1,  help="mass maximal 2")
	parser.add_option("--jump_mass", metavar = "1", type="float", nargs=1,  help="increment de la masse a chaque iteration")
	parser.add_option("--wanted_ratio", metavar = "1", type="float", nargs=1,  help="This ratio will define the range of SNR of the creating files")

	args,filenames=parser.parse_args()

	return args,filenames




args,filenames = parse_command_line()


m1_min = args.mass_min_1
m1_max = args.mass_max_1
m2_min = args.mass_min_2
m2_max = args.mass_max_2
jump_mass = args.jump_mass
m2_index = m2_min
initial_file = args.nbr_initial_file
wanted_ratio = args.wanted_ratio

nbr_m1 = math.ceil((m1_max-m1_min)/jump_mass)+1
nbr_m2 = math.ceil((m2_max-m2_min)/jump_mass)+1
nbr_total = int(nbr_m1 * nbr_m2)

index = 0
while m2_index <= m2_max:
	print(" ")
	print("################################################################################")
	print("#    We generate ",nbr_m1 , "files over the total of ",  nbr_total, "files   #")
	print("#      It remain ", nbr_total-(nbr_m1*index), "to generate                #")
	print("################################################################################")
	print(" ")
	
	subprocess.run(["./data_generation.py", "--nbr_initial_file=%s"%int(initial_file+index*nbr_m1), "--mass_min_1=%s"%m1_min, "--mass_max_1=%s"%m1_max, "--mass_min_2=%s"%m2_index, "--mass_max_2=%s"%m2_index, "--jump_mass=%s"%jump_mass, "--wanted_ratio=%s"%wanted_ratio])

	m2_index = m2_index + jump_mass
	index = index +1
