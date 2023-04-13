from functions import ir_domain_translator
import argparse
"""This function is designed to check the score of previously generated nt seq against the current implementation of the objective function
"""

parser = argparse.ArgumentParser(
                    prog = 'scorer',
                    description = 'Calculates the score of the given sequence against the most recent objective function',
                )

parser.add_argument("-i", metavar = "input ",  help="Input must be the output of ir_domain_translator",type = str)


args = parser.parse_args()

input_file = args.i

with open(input_file) as file:
    input = [line.rstrip() for line in file]

global d_seq
d_seq = input[4]

nt_seq = input[1]

steps = d_seq.count("l") + 1

nt_seqfp = []
extended_fp = []

for x in range(0,steps*2,2):
    nt_seqfp.append(input[x+6])
    extended_fp.append(input[x+7].split()[0])

print(d_seq)
print(nt_seq)
print(steps)
print(nt_seqfp)
print(extended_fp)





scores = ir_domain_translator.current_scores(nt_seqfp,extended_fp,nt_seq,d_seq)

print(sum(scores))