from functions import ir_domain_translator
from functions import nussinov
from functions import convert_functions
import argparse
import subprocess
import drt_out_parser 

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
seq,obj_fun = ir_domain_translator.rna_design(input,nussi_output)

print("----------",seq)
#seq = "CAAUGUUUCACCCUUUCAUUGAUGUCUACUAACUUCAUCAAUGAAAACCUUCCUCUUGGUUUUCAUUGAUGAAGCACAAAAACCGAGACAUCAAUGAAAGGGCACCCUGUGUGCCCUUUCAUUGAUGUCUCGGCCAGCGAGGUGCUUCAUCAAUGAAAACCAAGUCCUCUAUGGGCUUGGUUUUCAUUGAUGAAGCACUUC"

command = "echo " + seq + "| DrTransformer --name test"
result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


# Print the output of the command
# print(result.stdout.decode('utf-8'))

# Check for any errors
if result.stderr:
    print(result.stderr.decode('utf-8'))

domain_seq = input 

with open("domain_seq_out.txt") as f:
        path_out = [line.rstrip() for line in f]

with open('test.drf') as f:
    dr_out= [line.rstrip().split() for line in f]


fp = eval(path_out[-1])


pop,ext_fp = drt_out_parser.parse_drt_out(domain_seq,fp,dr_out)

with open("obj_fun_test.out","a") as out: 
     out.write("Objecitve function \n")
     out.write(obj_fun)
     out.write("\n")
     out.write("Sequence produced:\n")
     out.write(seq)
     out.write("\n")

     out.write("Following populations were achieved: \n")
     for x in range(len(ext_fp)):
        out.write(pop[x])
        out.write("  ")
        out.write(ext_fp[x])
        out.write("\n")