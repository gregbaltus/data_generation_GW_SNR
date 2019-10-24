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
	#parser.set_defaults(wanted_ratio=5)
	parser.set_defaults(ratio_min=3)
	parser.set_defaults(ratio_max=4)
	parser.set_defaults(ratio_jump=0.5)
	
	parser.add_option("--nbr_initial_file", metavar = "1", type="int", nargs=1,  help="nbr of GW initial. If you want more GW but conserve what you have")
	parser.add_option("--mass_min_1", metavar = "1", type="float", nargs=1,  help="mass minimal 1")
	parser.add_option("--mass_max_1", metavar = "1", type="float", nargs=1,  help="mass maximal 1")
	parser.add_option("--mass_min_2", metavar = "1", type="float", nargs=1,  help="mass minimal 1")
	parser.add_option("--mass_max_2", metavar = "1", type="float", nargs=1,  help="mass maximal 2")
	parser.add_option("--jump_mass", metavar = "1", type="float", nargs=1,  help="increment de la masse a chaque iteration")
	#parser.add_option("--wanted_ratio", metavar = "1", type="float", nargs=1,  help="This ratio will define the range of SNR of the creating files")
	parser.add_option("--ratio_min", metavar = "1", type="float", nargs=1,  help="minimal ratio")
	parser.add_option("--ratio_max", metavar = "1", type="float", nargs=1,  help="maximal ratio")
	parser.add_option("--ratio_jump", metavar = "1", type="float", nargs=1,  help="increment of the ratio")

	args,filenames=parser.parse_args()
	
	return args,filenames




args,filenames = parse_command_line()


m1_min = args.mass_min_1
m1_max = args.mass_max_1
m2_min = args.mass_min_2
m2_max = args.mass_max_2
jump_mass = args.jump_mass
initial_file = args.nbr_initial_file

tab_ratio = []
icr = 0.5
ratio_min = args.ratio_min
ratio_max = args.ratio_max
ration_jump = args.ratio_jump
nbr_ratio = int(math.ceil(ratio_max-ratio_min)/0.5+1)
ratio = ratio_min
for i in range(0,nbr_ratio):
	tab_ratio = tab_ratio + [ratio]
	ratio = ratio+icr
numero_folder=0
for i in tab_ratio:
        print(" ")
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print("@   creation of the forlder ",numero_folder, "on ", nbr_ratio, "folders   @")
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print(" ")
        numero_folder+=1
        subprocess.run(["./genrate_data.py", "--mass_min_1=%s"%m1_min, "--mass_max_1=%s"%m1_max, "--mass_min_2=%s"%m2_min, "--mass_max_2=%s"%m2_max, "--jump_mass=%s"%jump_mass, "--wanted_ratio=%s"%i])
        subprocess.run(["cp", "-r","data/", "data_ratio%s"%i])
        subprocess.run(["./clean_data.sh"])
