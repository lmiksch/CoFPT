from functions import ir_domain_translator
import argparse
import RNA as RNA 
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
print(" ")





#scores = ir_domain_translator.rstd_objective(nt_seq,True)

def score_info(sequence):   
        """print("Domain Seq from input file: ",input[3])
    input = input[3]
        split_nt_sequence = []
        #creates a list of list where each sublist i corresponds to the sequence at transcription step i 
        l_pointer = 0
        for z in split_seq:
            r_pointer = l_pointer + cv.d_length(z)
            
            split_nt_sequence.append(sequence[l_pointer:r_pointer])
            l_pointer = r_pointer
        
        
        nt_path = []

        for x in range(len(UL_liste)):
            if UL_liste[x][0] == "l":
            
                
                nt_path.append("".join(split_nt_sequence[:x+1]))
        nt_path.append(sequence)
        """
        nt_path = nt_seqfp
        total = []
        print("extended fp",extended_fp)
        for x in range(1,len(extended_fp)):
            
            #prepare input for finpath 
            
            ss1 = extended_fp[x-1] + ("." * (len(nt_path[x])-len(nt_path[x-1])))
         
            efe = ir_domain_translator.constrained_efe(nt_path[x],extended_fp[x])
            fc = RNA.fold_compound(nt_path[x])
            
            fe = fc.eval_structure(extended_fp[x].replace(".","x"))
            mypath, barrier = ir_domain_translator.call_findpath(nt_path[x],ss1,extended_fp[x],0,30)

           
            if mypath != None:
                deltaE = abs(mypath[-1][1]) - abs(mypath[0][1])
                
            else: 
                deltaE = 99
                barrier = 99

    
            mse = 0
            global factor
            print("\n")
            print(extended_fp[x])
            print("fe", fe)
            print("efe",efe)
            print("barrier",barrier)
            obj_score = ir_domain_translator.objective_function(fe,efe,mse,barrier)[0]
            
            total.append(obj_score) 

        
        #calculate MSE of scores to keep all scores equal 


        total_mean = sum(total)/ len(total)

        squared_error = [(x - total_mean) ** 2 for x in total]
        #print("total after", total)   
        for i,score in enumerate(total):
            total[i] = score + squared_error[i]
        print("\n Squared deviations from mean")    
        print(squared_error)

        #print("total after", total)    
        return total

print("single scores")
print(score_info(nt_seq))
