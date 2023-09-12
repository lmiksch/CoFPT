from functions import ir_domain_translator
from functions import nussinov
from functions import convert_functions
import argparse
import subprocess
import drt_out_parser 
import csv 
import os

""" Pipeline which takes the domain level sequence and translates it into nucleotide sequence. DrTransformer with nt sequence and the output gets evaluated by drt_out_parser. 

    Input: 
        -i domain_seq_out.txt
"""


"""This module converts the domain level sequence into a nucleotide level sequence.

"""

#generates files and folder structure 
try: 
    with open("testing.tsv","r") as test:
        lines = test.readlines()
    used_values = []
    for line in lines:
        values = line.strip().split("\t")
        used_values.append(values[0])

    current_index = int(used_values[-1]) + 1

except: 
    with open('testing.tsv', "w") as file:
        file.write("index\taverage_pop\tobjective_function\tscore\n")
        

    current_index = 0

if not os.path.exists("testing"):
     os.makedirs("testing")


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
seq,score = ir_domain_translator.rna_design(input,nussi_output,"testing/" + str(current_index))


command = "echo " + seq + "| DrTransformer --name " + str(current_index) + "_test"
result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

if result.stderr:
    print(result.stderr.decode('utf-8'))

domain_seq = input 

with open("domain_seq_out.txt") as f:
        path_out = [line.rstrip() for line in f]

with open(str(current_index) + '_test.drf') as f:
    dr_out= [line.rstrip().split() for line in f]


fp = eval(path_out[-1])


pop,ext_fp = drt_out_parser.parse_drt_out(domain_seq,fp,dr_out)



avg_pop = round(sum(float(x) for x in pop)/len(pop),4)
s,obj_function = ir_domain_translator.objective_function(0,0,0,0)
with open("testing/" + str(current_index) + "_test.out","a") as out: 
     out.write("Objecitve function \n")
     out.write(obj_function)
     out.write("\n")
     out.write("Sequence produced:\n")
     out.write(seq)
     out.write("\n")

     out.write("Following populations were achieved with an average population of:")
     out.write(str(avg_pop))
     out.write("\n")
     for x in range(len(ext_fp)):
        out.write(pop[x])
        out.write("  ")
        out.write(ext_fp[x])
        out.write("\n")
     out.write("Score = ")
     out.write(str(score))
     out.write("\n")
           

#add current solutions to tsv file 
with open("testing.tsv", "a", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow([current_index, avg_pop, obj_function,score])