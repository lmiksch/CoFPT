import infrared as ir 
from infrared import rna
import RNA 
from functions import convert_functions as cv
from functions import nussinov
import random
import math

def objective_function(fe,efe,mse,barrier):
    """
    Function exists to only change one thing when changing the obj function
    """

    obj_fun = "(fe - efe)**2 + barrier*0.00001"
    score = eval(obj_fun.format(fe=fe,efe=efe,barrier=barrier))
    return score,obj_fun

def couple(pair):
	if pair[0].upper() == pair[1].upper() and pair[0] != pair[1]:
		return True

def identicaldomains(x):
    if x == True:

        ir.def_constraint_class(
            "IdenticalDomains",
            lambda i,j: [i,j],
            lambda x,y: x == y 
        )
    print("Imported IdenticalDomains constraint")



#could be changed to have custom input for the different domains
def d_length(domain):
    if domain[0] == "b" or domain[0] == "B":
        return 5
    elif domain == "l":
        return 5
    return 3

def identical_domains_constraint(domain,UL_seq,model):
    """Identical_domains_constraint

    Applies the constraint to our model, that the value of indices of the same domain must be equal.

    Args:
        domain: domain name
        UL_seq: upper lower case domain sequence
        model: model created in the infrared sampler
    """    

    for x in range(len(UL_seq)):
        if UL_seq[x] == domain:
            i_pointer = 0
            j_pointer = i_pointer

            for u in UL_seq[:x]:
                
                i_pointer += d_length(u)

            j_pointer = i_pointer    

            for q in range(len(UL_seq[:x]),len(UL_seq)):
           
                if domain not in UL_seq[x+1:]:
                    return    
                if UL_seq[q] == domain and j_pointer != i_pointer:
                   
                    break
                j_pointer += d_length(UL_seq[q])

        
            if i_pointer > j_pointer:
                return
            for z in range(d_length(domain)):
                if i_pointer != j_pointer:
                    
                    model.add_constraints(IdenticalDomains(i_pointer,j_pointer))
                    
                i_pointer += 1
                j_pointer += 1


def one_domains(UL_seq):
    split_seq = "".join(set(UL_seq))
    return split_seq

def domain_path_to_nt_path(path,UL_seq):
    """ Domain_path_to_nt_path

    Takes the domain level path and extends the number of the base pairings corresponding to their length: ( with length 5 --> (((((

    args:
        path (list): domain level path 
        UL_seq (list): upper lower case sequence 
    
    Returns:
        ext_path (list): where each sublist corresponds to the extended domain path
    """

    ext_path = [[] for x in path]
    
    for x in range(len(path)):
        
        
        for z in range(len(path[x][0])):
       
            ext_path[x].append(path[x][0][z] * d_length(UL_seq[z]))
        ext_path[x] = "".join(ext_path[x])


    return ext_path

def call_findpath(seq, ss1, ss2, md, fpw, mxb = float('inf')):
    """ Call ViennaRNA findpath. Modified from DrTransformer
    """
    fc = RNA.fold_compound(seq)
    
    if mxb == float('inf'):
        path = fc.path_findpath(ss1, ss2, width = fpw)
        
    else:
        e1 = round(fc.eval_structure(ss1), 4)
        dcal_bound = int(round((mxb + e1)))
        path = fc.path_findpath(ss1, ss2, maxE = dcal_bound, width = fpw)
        
    del fc

    if len(path):
        mypath = []
        barrier = None
        for step in path:
            struct = step.s
            energy = float(round(step.en,4))
            mypath.append((struct, energy))
            if barrier is None or barrier < energy:
                barrier = energy
        barrier -= float(round(path[0].en,4))
        del step, path # potentially a good idea.
        return mypath, barrier
    return None, None




def add_folding_path_constraint(path,UL_seq,model):
    """
    Takes the path calculated from nussinov for the sequences and adds those pair constraints to the model

    args: 
        path (list): output path of the nussinov algorithm
        UL_seq (list): upper lower case sequence
        model : model created in infrared sampler
    """


    ext_path = domain_path_to_nt_path(path,UL_seq)

    for x in ext_path: 
        cons = []
        bps = rna.parse(x)
        cons = [rna.BPComp(i,j) for (i,j) in bps]
        

        model.add_constraints(cons)





def current_scores(nt_seqfp,extended_fp,seq,d_seq):
    """Calculates the final score of each sequence at each transcription step 
    
    """
    seq_path = cv.split_ntseq_to_domainfp(seq,d_seq)
    
    scores = [0,]
    for x in range(1,len(nt_seqfp)):

        #prepare input for finpath 
        fc = RNA.fold_compound(seq_path[x])
        mfe_curr, e_curr = fc.mfe()

        fc = RNA.fold_compound(seq_path[x-1])
        mfe_prev , e_prev = fc.mfe()
        ss1 = mfe_prev + ("." * (len(seq_path[x])-len(seq_path[x-1])))
       
        ss1 = extended_fp[x-1] + ("." * (len(seq_path[x])-len(seq_path[x-1])))
       
        efe = constrained_efe(seq_path[x],extended_fp[x])
        
        fc = RNA.fold_compound(seq_path[x])
        
        fe = fc.eval_structure(extended_fp[x].replace(".","x"))
        
        mypath, barrier = call_findpath(seq_path[x],ss1,extended_fp[x],0,30)
        
        
        #print(mypath)
        if mypath != None:
            deltaE = abs(mypath[-1][1]) - abs(mypath[0][1])
        else: 
            deltaE = 99
            barrier = 99
            

        E_1_fc = RNA.fold_compound(seq_path[x-1])
        E_1 = E_1_fc.eval_structure(extended_fp[x-1].replace(".","x"))
       

        E_i_fc = RNA.fold_compound(seq_path[x])
        E_i = E_i_fc.eval_structure(extended_fp[x].replace(".","x"))

        E_i1_fc = RNA.fold_compound(seq_path[x])
        E_i1 = E_i1_fc.eval_structure(extended_fp[x].replace(".","x"))

        deltaE = E_i - E_i1


        
        mse = (deltaE - (fe - E_1)/x)**2
        """
        print("fe",fe)
        print("efe",efe)
        print("barrier",barrier)
        """

        scores.append(objective_function(fe,efe,mse,barrier)[0])
        
        
    return scores
        



def find_domain_length(domain):

    global domain_length
    
    length = domain_length[domain]
    
    return length

def mc_optimize(model, objective, steps, temp, start=None):
        sampler = ir.Sampler(model)
        cur = sampler.sample() if start is None else start
        curval = objective(cur)
        best, bestval = cur, curval
        
        ccs = model.connected_components()
        weights = [1/len(cc) for cc in ccs]
        
        for i in range(steps):
            cc = random.choices(ccs,weights)[0]
            new = sampler.resample(cc, cur)
            newval = objective(new)
            print("\rCurrent Score: ",newval, end="")
            if (newval >= curval
                or random.random() <= math.exp((newval-curval)/temp)):
                cur, curval = new, newval
                if curval > bestval:
                    best, bestval = cur, curval
        
        return (best, bestval)

def constrained_efe(sequence,c):
        """Calculates the ensemble free energy of sequence
        
        """
        fc = RNA.fold_compound(sequence)
        # fc.hc_add_from_db(c) # not sure if it's necessary to add constraint
        return fc.pf()[1]



def rna_design(seq,path,out):
    """
    Takes the calculated sequence(With space annotation) and the calulated nussinovpath as input and designs a RNA sequence
    
    args: 
        seq: sequence generated by path_to_seq (space annotation)
        path: extended folding path 


    """
    print("Given sequence: ",seq)
    print("Given path: ",path)
    split_seq = seq.split()
    global d_seq
    d_seq = seq
  
    seqlen = 0
    no_space_seq = "".join(split_seq)

    UL_liste = cv.UL_list(split_seq)

    #identicaldomains(True)

    ir.def_constraint_class(
            "IdenticalDomains",
            lambda i,j: [i,j],
            lambda x,y: x == y,
            module  = __name__ 
        )

    

    #calculates Sequence length
    for x in range(len(UL_liste)):
        seqlen += d_length(UL_liste[x])
       

    print("Running Calculations")
    model = ir.Model(seqlen,4)


    

    
    add_folding_path_constraint(path,UL_liste,model)


    #applies constraint, that same domains should have the same sequence
    unique_domains = "".join(set(UL_liste))
    for domain in  unique_domains:

        if domain != "l":
            identical_domains_constraint(domain,UL_liste,model)




    
    #optimization

    extended_fp = domain_path_to_nt_path(path,UL_liste)

    #addition of energy to our model optional since this does not lead to better results
   

    
    """for x in extended_fp: 
            ss = rna.parse(x)
            model.add_functions([rna.BPEnergy(i, j, (i-1, j+1) not in ss)
        for (i,j) in ss], 'energy')
    energy = -0.5
    model.set_feature_weight(energy, 'energy')
    """
    
    
    #controlling gc content since adding energy will generate a gc bias
    GCcont = -0 
    gc_funs = [rna.GCCont(i) for i in range(seqlen)]
    model.add_functions(gc_funs, "gc")
    model.set_feature_weight(GCcont, "gc")
    






    def rstd_objective(sequence):
    
        split_nt_sequence = []
        
        
        #creates a list of list where each sublist i corresponds to the sequence at transcription step i 
        l_pointer = 0
        for z in split_seq:
            r_pointer = l_pointer + d_length(z)
            
            split_nt_sequence.append(sequence[l_pointer:r_pointer])
            l_pointer = r_pointer
        
        
        nt_path = []

        for  x in range(len(UL_liste)):
            if UL_liste[x] == "l":
            
                
                nt_path.append("".join(split_nt_sequence[:x+1]))
        nt_path.append(sequence)
        

        total = []
        
        for x in range(1,len(extended_fp)):

            #prepare input for finpath 
            
            ss1 = extended_fp[x-1] + ("." * (len(nt_path[x])-len(nt_path[x-1])))
         
            efe = constrained_efe(nt_path[x],extended_fp[x])
            fc = RNA.fold_compound(nt_path[x])
            
            fe = fc.eval_structure(extended_fp[x].replace(".","x"))
            mypath, barrier = call_findpath(nt_path[x],ss1,extended_fp[x],0,30)

           
            if mypath != None:
                deltaE = abs(mypath[-1][1]) - abs(mypath[0][1])
                
            else: 
                deltaE = 99
                barrier = 99
            
            E_1_fc = RNA.fold_compound(nt_path[0])
            E_1 = E_1_fc.eval_structure(extended_fp[0].replace(".","x"))
            

            E_i_fc = RNA.fold_compound(nt_path[x])
            E_i = E_i_fc.eval_structure(extended_fp[x].replace(".","x"))

            E_i1_fc = RNA.fold_compound(nt_path[x-1])
            E_i1 = E_i1_fc.eval_structure(extended_fp[x-1].replace(".","x"))

            deltaE = E_i - E_i1


            
            mse = (deltaE - (fe - E_i)/x)**2
            global factor


            #print("fe", fe)


            
            #print("efe",efe)
            #print("barrier",barrier)
            obj_score = objective_function(fe,efe,mse,barrier)[0]
            
            total.append(obj_score) 

        
        #calculate MSE of scores to keep all scores equal 


        total_mean = sum(total)/ len(total)

        squared_error = [(x - total_mean) ** 2 for x in total]
        #print("total after", total)   
        for i,score in enumerate(total):
            total[i] = score + squared_error[i]


        #print("total after", total)    
        return sum(total)


    #[rstd-optimize-call]
    objective = lambda x: -rstd_objective(rna.ass_to_seq(x))

    best, best_val = mc_optimize(model, objective,steps = 3000, temp = 0.04)



    print("\n")#
    #Output generation
    f = open(str(out) + "_IR_out.txt", "a")
    print("Calculated NT sequence:")
    f.write("Calculated NT sequence: \n")
    print(rna.ass_to_seq(best), -best_val)
    f.write(rna.ass_to_seq(best))
    f.write("\n")
    f.write("best score: ")
    f.write(str(-best_val))
    f.write(" \n")
    print("length = ",len(rna.ass_to_seq(best)))
    f.write("Length = ")
    f.write(str(len(rna.ass_to_seq(best))))
    f.write(" \n")
    ntseq_fp = cv.split_ntseq_to_domainfp(rna.ass_to_seq(best),d_seq)


    #Visualization
    
    
    final_scores = current_scores(ntseq_fp,extended_fp,str(rna.ass_to_seq(best)),d_seq)
    print("error 1")
    i = 0
    print("\n")   
    f.write(seq)
    f.write("\n")
    print(" ")
    print("Resulting Module Folding Path with calculated structure using RNAfold")   
    f.write("Resulting Module Folding Path with calculated structure using RNAfold") 
    f.write(" \n")
    for z in ntseq_fp:
        fc = RNA.fold_compound(z)
        (mfe_struct, mfe) = fc.mfe()
        print(z)
        print(mfe_struct)
        f.write(z)
        f.write(" \n")
        f.write(mfe_struct)
        f.write("  Score: ")
        f.write(str(final_scores[i]))
        f.write(" \n")
        print("")
        i += 1
    f.write(" \n")
    print("")
        
    #output include all folding path structures so (),().,()(), with rnafold
    print(cv.extended_domain_path(UL_liste))
    f.write(cv.extended_domain_path(UL_liste))
    
    return rna.ass_to_seq(best), best_val

def split_ntseq_to_domainfp(nt_seq,domain_seq):
    """split_ntseq_to_domainfp converts the nucleotide output into a list with sublists where each sublist corresponds to the sequence at the respective transcription step 

    Args:
        nt_seq (str): nucleotide sequnce
        domain_seq (str): domain level sequence with space annotation 

    Returns: 
        nt_path (list): list with sublists where each sublist corresponds to the sequence at the respective transcription step 

    """

    split_seq = domain_seq.split()
        
    split_nt_sequence = []
    UL_list = cv.UL_list(domain_seq)
        
    l_pointer = 0
    for z in split_seq:
        r_pointer = l_pointer + d_length(z)
            
        split_nt_sequence.append(nt_seq[l_pointer:r_pointer])
        l_pointer = r_pointer
        
  
        
    nt_path = []

    for  x in range(len(UL_list)):
        if UL_list[x] == "l":
                
            nt_path.append("".join(split_nt_sequence[:x+1]))
    nt_path.append(nt_seq)

    return nt_path







if __name__ == "__main__":
    #example
    seq = " b   l  f* a* b* c* g*  l  d c  b  a e  l  e* a* b* c* d*  l  h g c  b  a f i  l  i* f* a* b* c* g* h*"
    path = [['..'], ['(...)...'], ['...(((...)))..'], ['(...)...(((((.))))).'], ['..(((((.(((((.)))))..)))))..'], ['(...)...(((((.))))).(((((((.)))))))']]
    



    rna_design(seq,path)


    