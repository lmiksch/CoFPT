import infrared as ir 
from infrared import rna
import RNA 
from convert_functions import *




def couple(pair):
	if pair[0].upper() == pair[1].upper() and pair[0] != pair[1]:
		return True


#def same_domains(seq):

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


ir.def_constraint_class(
    "IdenticalDomains",
    lambda i,j: [i,j],
    lambda x,y: x == y 
)

#could be changed to have custom input for the different domains
def d_length(domain):
    if domain[0] == "b" or domain[0] == "B":
        return 5
    elif domain == "l":
        return 3 
    return 3 

def identical_domains(domain,UL_seq,model):
     print("Domain",domain,"Seq:",UL_seq)
     for x in range(len(UL_seq)):
        if UL_seq[x] == domain:
            i_pointer = 0
            j_pointer = i_pointer

            for u in UL_seq[:x]:
                print("u = ",u)
                i_pointer += d_length(u)
            j_pointer = i_pointer    
            for q in UL_seq[x+1:]:
                print("q = ",q,d_length(q))
                j_pointer += d_length(q)
                if q == domain and j_pointer != i_pointer:
                   break
                if q == UL_seq[-1] and q != domain:
                    return
            for z in range(d_length(domain)):
                if i_pointer != j_pointer:
                    model.add_constraints(IdenticalDomains(i_pointer,j_pointer))
                    print(i_pointer,j_pointer)
                i_pointer += 1
                j_pointer += 1


def one_domains(UL_seq):
    split_seq = "".join(set(UL_seq))
    return split_seq
    




def find_domain_length(domain):

    global domain_length
    
    length = domain_length[domain]
    print("length = ", length)
    return length



def rna_design(seq,path):
    """
    Takes the calculated sequence(With space annotation) and the calulated nussinovpath as input and designs a RNA sequence
    
    """

    split_seq = seq.split()
  
    seqlen = 0
    no_space_seq = "".join(split_seq)

    UL_seq = convert_to_UL(no_space_seq)

    print(no_space_seq)

    #calculates Sequence length
    for x in range(len(split_seq)):
        seqlen += d_length(split_seq[x])
        print("seqlen   ",seqlen,d_length(split_seq[x]),split_seq[x])

    print("seqlen   ",seqlen)


    model = ir.Model(seqlen,4)
    bps = rna.parse(path)
    print(bps)
    cons = [rna.BPComp(i,j) for (i,j) in bps]

    model.add_constraints(cons)
    #applies constraint, that same domains should have the same sequence
    unique_domains = "".join(set(UL_seq))
    for domain in  unique_domains:
        print("Domain",domain)
        identical_domains(domain,UL_seq,model)

    print(path)
  

   
    sampler = ir.Sampler(model)

    show_td_info(sampler)


    samples = [sampler.sample().values() for i in range(10)]

    print("samples",samples)


    sequences = [rna.values_to_seq(s) for s in samples]
    print(sequences)
    
    #Visualization
    split_sequences = [ [] for x in sequences]
    s_pointer = 0
    for x in sequences:
        l_pointer = 0
        for z in split_seq:
            r_pointer = l_pointer + d_length(z)
            
            split_sequences[s_pointer].append(x[l_pointer:r_pointer])
            l_pointer = r_pointer
        s_pointer += 1
    print(split_sequences)
    for z in sequences:
        fc = RNA.fold_compound(z)
        (mfe_struct, mfe) = fc.mfe()
        print(z)
        print(mfe_struct)


if __name__ == "__main__":

    rna_design("a b c l a* b* c* l b","(((((((((((...)))))))))))........")   

