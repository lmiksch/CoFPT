from path_to_seq import *
from nussi2 import *
from path_generator import * 
from convert_functions import *
import argparse

parser = argparse.ArgumentParser(
                    prog = 'Path_to_Seq_test',
                    description = 'This programm test the path_to_seq progam. By calculating the folding path via nussinov algorithm and the comparing it to the input path',
                    epilog = 'Text at the bottom of help')

parser.add_argument("-length", metavar = "length of the test sequence", default= 4, help="Input for the length of the sequence which will be tested",type = int)
parser.add_argument("-t", metavar = "Number of Random Paths", default= 10, help="Input for how many random paths are generated",type = int)

args = parser.parse_args()

len = args.length / 2

len = int(len)


if args.length%2 == 0:
    seq = "aA" * len

else:
    seq = "aA" * len + "a"


calculated_sequences = []
possible_paths = generate_foldingpaths(seq)

paths = []
print("Generating Random Paths")
for k in range(args.t):
    paths.append(get_random_path(possible_paths))
    
for path in paths:
    print(path)

i = 0
print("Checking if paths are equal")
for path in paths:
    pertable_path = path_to_pertablepath(path)
    calc_seq = fold_to_seq(pertable_path)
    conv_seq = convert_pts_out_to_nussi_in(calc_seq)
    
    calculated_sequences.append(conv_seq)
    nussinov_path = nussinov(conv_seq)
    calculated_path = module_folding_path(nussinov_path,convert(conv_seq))
    if calculated_path == path:

        print("works")
    elif calculated_path != nussinov_path:
        print("Error: Following Path does not match calulated path")
        print(path,"Input")
        print(calculated_path,"calculated from nussinov")
        print("calculated sequence:",calc_seq)

    i += 1    

print("done")


