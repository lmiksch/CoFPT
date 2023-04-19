from functions import path_to_seq 
from functions import nussinov 
from functions import  path_generator 
from functions import convert_functions
import argparse
# -*- coding: utf-8 -*-
"""This module was created to test the path_to_seq function. 

"""



parser = argparse.ArgumentParser(
                    prog = 'Path_to_Seq_test',
                    description = 'This programm test the path_to_seq progam. By calculating the folding path via nussinov algorithm and the comparing it to the input path',
                    epilog = 'Text at the bottom of help')

parser.add_argument("-length", metavar = "length of the test sequence", default= 6, help="Input for the length of the sequence which will be tested",type = int)
parser.add_argument("-t", metavar = "Number of Random Paths", default= 10, help="Input for how many random paths are generated",type = int)

args = parser.parse_args()

length = args.length / 2

length = int(length)


if args.length%2 == 0:
    seq = "aA" * length

else:
    seq = "aA" * length + "a"

def check_paths(seq):
    succesfull = 0
    calculated_sequences = []
    possible_paths = path_generator.generate_foldingpaths(seq)

    paths = []  
    print("Generating Random Paths")
    for k in range(args.t):
        paths.append(path_generator.get_random_path(possible_paths))

    print("Checking if paths are equal")
    for path in paths:
        pairtable_path = convert_functions.path_to_pairtablepath(path)
        calc_seq = path_to_seq.fold_to_seq(pairtable_path)
        print("Calculated Sequence = ",calc_seq)
        

        nussinov_path = nussinov.nussinov_modules(calc_seq)
        if len(nussinov_path) != 1: 
            print("Warning: There are more than one optimal folding paths for this sequence. Reconsider your input")
        if nussinov_path[0] == path and len(nussinov_path) == 1:
            print("works:",pairtable_path)
            print(nussinov_path,"calculated from nussinov")
            print("calculated sequence:",calc_seq) 
            succesfull += 1

        else:
            print("Error: Following Path does not match calulated path")
            print(path,"Input",pairtable_path)
            print("Following Paths were calculated from the Nussinov Algorithm")
            for path in nussinov_path:
                print(path)
            print("calculated sequence:",calc_seq)
        print("\n")
    print("Number of Succesfull attempts = ", succesfull,"/",args.t)

def check_path(path): 
    """Takes a dot-bracket notation path and calculates the domain level sequence and then checks with nussi if the calculated sequence is equal to the input sequence"""
    pairtable_path = convert_functions.path_to_pairtablepath(path)
    calc_seq = path_to_seq.fold_to_seq(pairtable_path)
    
    nussinov_path = nussinov.nussinov(calc_seq)
    converted_seq = convert_functions.convert_UL_list(calc_seq.split())
    calculated_path = convert_functions.only_b_domainfp(converted_seq,nussinov_path)
    if calculated_path == path:
        print("works:",pairtable_path)
        print(calculated_path,"calculated from nussinov")
        print("calculated sequence:",calc_seq) 
       
    elif calculated_path != path:
        print("Error: Following Path does not match calulated path")
        print(path,"Input",pairtable_path)
        print(calculated_path,"calculated from nussinov")
        print("calculated sequence:",calc_seq)
    
check_paths(seq)


#print(nussinov("ba*b*c*ecbadg*d*a*b*c*e*f*fecbadg"))
#nussinov_path = nussinov(" blf*a*b*c*g*ldcbaele*a*b*c*d*lhgcbafili*f*a*b*c*g*h*")
#print(nussinov_path)
#print(nussinov_modules(nussinov_path))
