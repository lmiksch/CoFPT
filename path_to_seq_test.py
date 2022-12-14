from path_to_seq import *
from nussi2 import *
from path_generator import * 
from convert_functions import *
import argparse
# -*- coding: utf-8 -*-
"""This module was created to test the path_to_seq function.

"""



parser = argparse.ArgumentParser(
                    prog = 'Path_to_Seq_test',
                    description = 'This programm test the path_to_seq progam. By calculating the folding path via nussinov algorithm and the comparing it to the input path',
                    epilog = 'Text at the bottom of help')

parser.add_argument("-length", metavar = "length of the test sequence", default= 5, help="Input for the length of the sequence which will be tested",type = int)
parser.add_argument("-t", metavar = "Number of Random Paths", default= 10, help="Input for how many random paths are generated",type = int)

args = parser.parse_args()

len = args.length / 2

len = int(len)


if args.length%2 == 0:
    seq = "aA" * len

else:
    seq = "aA" * len + "a"

def check_paths(seq):
    succesfull = 0
    calculated_sequences = []
    possible_paths = generate_foldingpaths(seq)

    paths = []  
    print("Generating Random Paths")
    for k in range(args.t):
        paths.append(get_random_path(possible_paths))

    print("Checking if paths are equal")
    for path in paths:
        pertable_path = path_to_pertablepath(path)
        calc_seq = fold_to_seq(pertable_path)
        conv_seq = convert_pts_out_to_nussi_in(calc_seq)
    
        calculated_sequences.append(conv_seq)
        nussinov_path = nussinov(conv_seq)
        calculated_path = module_folding_path(nussinov_path,convert(conv_seq))
        if calculated_path == path:
            print("works:",pertable_path)
            print(calculated_path,"calculated from nussinov")
            print("calculated sequence:",calc_seq) 
            succesfull += 1

        else:
            print("Error: Following Path does not match calulated path")
            print(path,"Input",pertable_path)
            print(calculated_path,"calculated from nussinov")
            print("calculated sequence:",calc_seq)
        print("\n")
    print("Number of Succesfull attempts = ", succesfull,"/",args.t)

def check_path(path): 
    """Takes a dot-bracket notation path and calculates the domain level sequence and then checks with nussi if the calculated sequence is equal to the input sequence"""
    pertable_path = path_to_pertablepath(path)
    calc_seq = fold_to_seq(pertable_path)
    conv_seq = convert_pts_out_to_nussi_in(calc_seq)
    
    nussinov_path = nussinov(conv_seq)
    calculated_path = module_folding_path(nussinov_path,convert(conv_seq))
    if calculated_path == path:
        print("works:",pertable_path)
        print(calculated_path,"calculated from nussinov")
        print("calculated sequence:",calc_seq) 
       
    elif calculated_path != path:
        print("Error: Following Path does not match calulated path")
        print(path,"Input",pertable_path)
        print(calculated_path,"calculated from nussinov")
        print("calculated sequence:",calc_seq)
    
check_paths(seq)


#print(nussinov("ba*b*c*ecbadg*d*a*b*c*e*f*fecbadg"))
#nussinov_path = nussinov("bla*b*c*lecbadlg*d*a*b*c*e*f*lfecbadg")
#print(nussinov_path)
#print(module_folding_path(nussinov_path,convert("bla*b*c*lecbadlg*d*a*b*c*e*f*lfecbadg")))
