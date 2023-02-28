from functions import ir_domain_translator
from functions import nussinov
from functions import convert_functions
import argparse

# -*- coding: utf-8 -*-
"""This module converts the domain level sequence into a nucleotide level sequence.

"""
#input parsing prolly just from commandline a domain seq or input output file of domain_seq 


parser = argparse.ArgumentParser(
                    prog = 'domain_translator',
                    description = 'Translates a domain level sequence into a nucleotide sequence.',
                )

parser.add_argument("-i", metavar = "input ",  help="Input must be the output of doman_seq_generator otherwise a domain level seq needs to be specified.",type = str)


args = parser.parse_args()

input_file = args.i
try:
    with open(input_file) as file:
        input = [line.rstrip() for line in file]

    print("Domain Seq from input file: ",input[3])
    input = input[3]
except:
    args.i = None

if args.i == None: 
	input = input("Input here a domain level sequence:")
	print("input =", input )




domain_seq = input


print("Calculating nussinov")
nussi_output = nussinov.nussinov(domain_seq)



# Nussinov to get extended folding path 
print("\n", "Running Translation")
print(input)
ir_domain_translator.rna_design(input,nussi_output)









