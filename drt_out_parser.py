from functions import convert_functions as cv
import argparse


parser = argparse.ArgumentParser(
                    prog='DrT_parser',
                    description='Parses the output of DrTransformer and checks if the desired structures at each transcription step are present',
                    )
parser.add_argument("-drt_out",type = str, help="DrTransformer output file you want analyze", default="NoName.drf")
parser.add_argument("-domain",type = str, help="domain_seq_generator output file you want analyze", default="domain_seq_out.txt")


#reading files
with open('NoName.drf') as f:
    dr_out= [line.rstrip().split() for line in f]

with open("domain_seq_out.txt") as f:
    path_out = [line.rstrip() for line in f]

#dr_out[][2] = occupancy 
#dr_out[][3] = structure
#dr_out[4] = energy


domain_seq = path_out[3]
fp = eval(path_out[-1])


#write function which extends fp for comparison to drt output
extended_fp = cv.extended_fp_path(fp,domain_seq)

populations = []
i = 0
for x in range(1,len(dr_out)):
    if len(dr_out[x][3]) > len(extended_fp[i]):
        populations.append("0")
        i += 1
    if dr_out[x][3] == extended_fp[i]:
        populations.append(dr_out[x][2])
        i += 1
        if i == len(extended_fp):
            break

print("Populations", "Structure")
for x in range(len(extended_fp)):
    print(populations[x],extended_fp[x])
