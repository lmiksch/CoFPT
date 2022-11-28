#transforms motif seq into simple structural sequence 


def motifs_to_seq(nt,motif_seq):

    a = "("*nt
    A = ")"*nt
    seq = ""

    for x in motif_seq:
        if x == "a": 
            seq += a
        elif x == "A":
            seq += A    
        else: 
            seq += "."*nt
        
    return(seq)

print(motifs_to_seq(5,"abAcadAeafAgahAiajAkalAmanAo"))
