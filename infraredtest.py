import infrared as ir 
from infrared import rna
import RNA 
from convert_functions import *
from nussinov import *
import random
import math






def couple(pair):
	if pair[0].upper() == pair[1].upper() and pair[0] != pair[1]:
		return True


def show_td_info(sampler,width=600):
    td = sampler.td
    print("tree width =", td.treewidth())
    print("bags =", td.bags)
    print("edges =", td.edges)
    
    tmpfile="tmp_out.png"
    sampler.plot_td(tmpfile,'png')
    from IPython.display import Image
    return Image(filename=tmpfile,width=width)


class domain():
    def __init__(self,name,star,length):
        self.name = name
        self.star = star
        self.length = length

def identicaldomains(x):
    if x == True:

        ir.def_constraint_class(
            "IdenticalDomains",
            lambda i,j: [i,j],
            lambda x,y: x == y 
        )
    print("Imported IdenticalDomains constraint")

identicaldomains


#could be changed to have custom input for the different domains
def d_length(domain):
    if domain[0] == "b" or domain[0] == "B":
        return 5
    elif domain == "l":
        return 5
    return 3 

def identical_domains_constraint(domain,UL_seq,model):
    

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

    ext_path = [[] for x in path]
    
    for x in range(len(path)):
        
        
        for z in range(len(path[x][0])):
       
            ext_path[x].append(path[x][0][z] * d_length(UL_seq[z]))
        ext_path[x] = "".join(ext_path[x])


    return ext_path






def add_folding_path_constraint(path,UL_seq,model):
    """
    Takes the path calculated from nussinov for the sequences and adds those pair constraints to the model
    
    """
    ext_path = domain_path_to_nt_path(path,UL_seq)

    for x in ext_path: 
        cons = []
        bps = rna.parse(x)
        cons = [rna.BPComp(i,j) for (i,j) in bps]
        

        model.add_constraints(cons)





def find_domain_length(domain):

    global domain_length
    
    length = domain_length[domain]
    
    return length



def rna_design(seq,path):
    """
    Takes the calculated sequence(With space annotation) and the calulated nussinovpath as input and designs a RNA sequence
    
    """
    print("Given sequence: ",seq)
    print("Given path: ",path)
    split_seq = seq.split()
  
    seqlen = 0
    no_space_seq = "".join(split_seq)

    UL_seq = convert_to_UL(no_space_seq)

    identicaldomains(True)
    

    #calculates Sequence length
    for x in range(len(split_seq)):
        seqlen += d_length(split_seq[x])
        #print("seqlen   ",seqlen,d_length(split_seq[x]),split_seq[x])

    #print("seqlen   ",seqlen)

    print("Running Calculations")
    model = ir.Model(seqlen,4)


    add_folding_path_constraint(path,UL_seq,model)


    #applies constraint, that same domains should have the same sequence
    unique_domains = "".join(set(UL_seq))
    for domain in  unique_domains:
        #print("Domain",domain)
        identical_domains_constraint(domain,UL_seq,model)




    
    #optimization

    extended_fp = domain_path_to_nt_path(path,UL_seq)

 
    print("extended fp",extended_fp)
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
            if (newval >= curval
                or random.random() <= math.exp((newval-curval)/temp)):
                cur, curval = new, newval
                if curval > bestval:
                    best, bestval = cur, curval
    
        return (best, bestval)
    
    def constrained_efe(sequence,c):
        fc = RNA.fold_compound(sequence)
        fc.hc_add_from_db(c)
        return fc.pf()[1]
 
    def rstd_objective(sequence):
    
        split_nt_sequence = []
        
        
        l_pointer = 0
        for z in split_seq:
            r_pointer = l_pointer + d_length(z)
            
            split_nt_sequence.append(sequence[l_pointer:r_pointer])
            l_pointer = r_pointer
        
        
        nt_path = []

        for  x in range(len(UL_seq)):
            if UL_seq[x] == "l":
            
                
                nt_path.append("".join(split_nt_sequence[:x+1]))
        nt_path.append(sequence)
        
        total = 0
        

        for x in range(1,len(extended_fp)):
            #print(nt_path[x])
            #print(extended_fp[x])
            efe = constrained_efe(nt_path[x],extended_fp[x])
            fc = RNA.fold_compound(nt_path[x])
            
            fe = fc.eval_structure(extended_fp[x])
            
            total += fe - efe
        
        
        return total


    #[rstd-optimize-call]
    objective = lambda x: -rstd_objective(rna.ass_to_seq(x))

    best, best_val = mc_optimize(model, objective,steps = 500, temp = 0.04)

    print("done")
    print(" ")
    print("Calculated NT sequence:")
    print(rna.ass_to_seq(best), -best_val)
    print("length = ",len(rna.ass_to_seq(best)))
    
    ntseq_fp = split_ntseq_to_domainfp(rna.ass_to_seq(best),seq)


    #Visualization
   

    print("\b")   
    
    print("")
    print("Resulting Module Folding Path with calculated structure using RNAfold")    
    for z in ntseq_fp:
        fc = RNA.fold_compound(z)
        (mfe_struct, mfe) = fc.mfe()
        print(z)
        print(mfe_struct)
    #output include all folding path structures so (),().,()(), with rnafold


def split_ntseq_to_domainfp(nt_seq,domain_seq):

        split_seq = domain_seq.split()
        
        split_nt_sequence = []
        UL_seq = convert_to_UL(domain_seq)
        
        l_pointer = 0
        for z in split_seq:
            r_pointer = l_pointer + d_length(z)
            
            split_nt_sequence.append(nt_seq[l_pointer:r_pointer])
            l_pointer = r_pointer
        
        print("Sequences Split up based on domains",split_nt_sequence)
        
        nt_path = []

        for  x in range(len(UL_seq)):
            if UL_seq[x] == "l":
                
                nt_path.append("".join(split_nt_sequence[:x+1]))
        nt_path.append(nt_seq)

        return nt_path

if __name__ == "__main__":
    
    #example1 for path = [['.'],['()'],['.()'],['()()'],['.(())'],['()()()'],['.()(())'],['(()(()))']]
    seq = "b   l  f* a* b* c* g*  l  d c  b  a e  l  j* e* a* b* c* d* k*  l  h g c  b  a f i  l  i* f* a* b* c* g* h*  l  k d c  b  a e j  l   b*"
    path = [['..'], ['(...)...'], ['...(((...)))..'], ['(...)...(((((..)))))..'], ['..(((((.(((((..)))))...)))))..'], ['(...)...(((((..)))))..(((((((.))))))).'], ['...(((...)))..(((((((.(((((((.))))))).))))))).'], ['(..(((...)))..(((((((.(((((((.))))))).))))))).)']]
    #rna_design(seq,path)

    #"CACGCAUCCCACCCGCGUGACUUCAAUCCAUAGUCACGCGGGCCCAUCGAAGGGCCCGCGUGACUAUGCCUAUCUCGUGAAGUCACGCGGGUGGUCCAUCGGGCCACCCGCGUGACUUCACGAAUCGGGCAUAGUCACGCGGGCCCUUCAUCGCGUG"


    #example [['.', '()', '.()', '(())', '.()()', '(()())']]
    #seq = "b   l  a* b* c*  l  c  b  a  l  d* b* e*  l  e  b  d  l   b* "
    #path = [['..'], ['(..)..'], ['..(((.))).'], ['(.(((.)))..)..'], ['..(((.))).(((.))).'], ['(.(((.))).(((.))).)']]


    rna_design(seq,path)


    