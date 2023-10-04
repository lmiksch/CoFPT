from functions import ir_domain_translator 
from functions import convert_functions as cv
import drt_out_parser
import argparse
import subprocess
import csv
import os


d_seq_1 = "b l1 b* a* l2 a b c d l3 d* c* b*"

path_1 = [[".."],["(.).."],["..((.))..."],["(.)...(((.)))"]]

#extendet_path = cv.extended_fp_path(path,d_seq)



d_seq_2 = "b l b* a* l a b c d l d* c* b* l a b e f l f* e* b*" #bbbbblllllBBBBBAAADDDlllllaaabbbbbccclllllCCCBBBBBllllldddaaabbbbbeeelllllEEEBBBBB

d_seq_3 = "b l b* a* d* l a b c l c* b* l d a b e l e* b* " 
path_2 = [[".."],["(.)..."],["..((..)).."],["(.)....((.))."],["..(((..((.)).))).."],["(.)....((.))...((.))"]]



d_seq_4 = "b l b* a* l a b l b* c* l c b d e l e* d* b*"

path_4 = [[".."],["(.).."],["..((.))."],["(.((.)).).."],["..((.)).((.))..."],["(.((.)).)...(((.)))"]]


d_seq_5 = "b a l b* c* l b d l a* b* l c b l d* b*"
path_5 = [["..."],["(..).."],["(..)....."],["((.(..)..))."],["...((.(...).))."],["(.)..((..(..).))"]]


d_seq_6 = "b d l b* a* l a b c l c* b* a* e* l e a b c l d* b*"

path_6 = [["..."],["(..).."],["...((.)).."],["(..)..(((.))).."],["(..)..(((.)))......."],["((.((.))..((((.)))).))"]]


path_7 = [['..'], ['(.)..'], ['..((.))..'], ['(.)..(((.)))..'], ['..((.))..((((.))))..'], ['(.)..(((.)))..(((((.)))))..'], ['..((.))..((((.))))..((((((.))))))..'], ['(.)..(((.)))..(((((.)))))..(((((((.)))))))..'], ['..((.))..((((.))))..((((((.))))))..((((((((.))))))))..'], ['(.)..(((.)))..(((((.)))))..(((((((.)))))))..(((((((((.)))))))))..'], ['..((.))..((((.))))..((((((.))))))..((((((((.))))))))..((((((((((.)))))))))).'], ['(.)..(((.)))..(((((.)))))..(((((((.)))))))..(((((((((.))))))))).......(.....)']]

d_seq_7 = "b   l   b* a*  l  a b c  l  c* b* a* d*  l  d a b c e  l  e* c* b* a* d* f*  l  f d a b c   e g  l  g* e* c* b* a* d* f* h*  l  h f d a b c   e   g i  l  i* g* e* c* b* a* d* f* h* j*  l  j h f d a b c   e   g   i  l   b* "





d_seq_6_steps_new = "b   l   b* a* g* l g a b c  l  c* b* a* d*  l  d a b c e  l  e* c* b* a* d* f*"


"""
For testing nothing step 

ex: 
[["."],["()"],["()."],["(())"]]

seq1: a b c b* a* b c* b* a* 
seq2: b a c a* b* b c* a* b* 
"""

d_seq_8 = "a b c l b* a* l b l c* b* a*" #seq1 
d_path_8 = [["...."],["((..))."],["((..))..."],["(((.(..).)))"]]


d_seq_9 = "b a c l a* b* l b l c* a* b*" #seq2  
d_path_9 = [["...."],["((..))."],["((..))..."],["(((..(.).)))"]]


d_seq_10 = "c b a l a* b* l b l a* b* c*"
d_path_10 = [["...."],[".((.))."],[".((.))..."],["(((..(.).)))"]]

def custom_path(d_seq,path):
    """ Pipeline which executes following steps:
        - designs a nt-sequence based on a given domain level sequence and folding path
        - Simulates CoTranscriptional folding path using DrTransformer
        - Evalutes the output and saves it 

    Args:
        d_seq(str): Domain level sequence in star annotation
        path(list): List where each sublist corresponds to the domain level structure at the corresponding folding step    
    """

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

    seq,score =  ir_domain_translator.rna_design(seq= d_seq, path=path, out= current_index)

    command = "echo " + seq + "| DrTransformer --name " + str(current_index) + "_test"
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


    if result.stderr:
        print(result.stderr.decode('utf-8'))


    with open(str(current_index) + '_test.drf') as f:
        dr_out= [line.rstrip().split() for line in f]


    pop,ext_fp = drt_out_parser.parse_drt_out(d_seq,path,dr_out)


    avg_pop = round(sum(float(x) for x in pop)/len(pop),4)
    s,obj_function = ir_domain_translator.objective_function(0,0,0,0)
    with open("testing/" + str(current_index) + "_test.out","a") as out: 
        out.write("Objective function \n")
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


    with open("testing.tsv", "a", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([current_index, avg_pop, obj_function,score])

d_seq_6_steps = "z l b   l   b* a*  l  a b c  l  c* b* a* d*  l  d a b c e  l  e* c* b* a* d* f*"
path_6_steps = [[".."],['(.).'], ['..(.)..'], ['(.).((.))..'], ['..(.)..(((.)))..'], ['(.).((.))..((((.))))..'], ['..(.)..(((.)))..(((((.))))).']]

path_8_steps = [['..'],["(.)."], ['..(.)..'], ['(.).((.))..'], ['..(.)..(((.)))..'], ['(.).((.))..((((.))))..'], ['..(.)..(((.)))..(((((.)))))..'], ['(.).((.))..((((.))))..((((((.))))))..'], ['..(.)..(((.)))..(((((.)))))..(((((((.))))))).']]

d_seq_8_steps = "z l b   l   b* a*  l  a b c  l  c* b* a* d*  l  d a b c e  l  e* c* b* a* d* f*  l  f d a b c   e g  l  g* e* c* b* a* d* f* h*"

custom_path(d_seq_6_steps,path_6_steps)