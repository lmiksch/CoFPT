from functions import nussinov
from functions import convert_functions
import argparse
from functions import path_to_seq
import os 


"""Takes an input file and generates a domain based sequence and checks if it matches the input path given. 

"""
#input parsing
parser = argparse.ArgumentParser(
                    prog = 'domain_seq_generator',
                    description = 'Translates an abstact folding path into a domain level sequence and applies nussinov to the domain based sequence',
                )

parser.add_argument("-i", metavar = "input ",  help="Input file needs to be a txt document where each line corresponds to a transcription step",type = str)


args = parser.parse_args()

input_file = args.i

if input_file == None: 
    print("No input file specified. Please specify with the -i command an input file.")
    exit()


with open(input_file) as file:
    input = [line.rstrip() for line in file]


print("Following abstract folding path was given as input: ",input)

afp = input

if not int in input: 
    input = convert_functions.path_to_pairtablepath(input)

print("converted input: ", input)


#calculating domain based sequence
domain_seq = path_to_seq.fold_to_seq(input)

print("Domain Seq:", domain_seq)

#nussinov for calculting domain level folding path
mfp = nussinov.nussinov_modules(convert_functions.convert_pts_out_to_nussi_in(domain_seq))

print(mfp)
mfp = mfp[0]
#checking if calculated path matches input path
match = True
for x in range(len(mfp)):
    if mfp[x] != afp[x]:
        match = False
if match == False:
    print("Calulated Path does not match input path")
    print("Calculated Path: ", mfp)
    print("Given Path: ",afp)
else:
    print("Calculated path matches input path.")


with open("domain_seq_out.txt", "a") as output:

    output.write("Input abstract folding path \n")
    output.write(str(afp))
    output.write("\n")
    output.write("Domain Level Seq: \n")
    output.write(domain_seq)
    output.write("\n\n")
    output.write("Calculated Path from Nussinov:\n")
    output.write(str(mfp))



