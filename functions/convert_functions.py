# -*- coding: utf-8 -*-
"""
This module contains all the different convert functions which can be used to convert outputs into different annotations. 
"""


def path_to_pairtablepath(path):
    """ Function which turns a structural path in dot bracket annotation into a pairtable path

    args: 
        path(list): List of the folding path in dot bracket
    output: 
        pairtable_path(List): List of the folding path using the pairtable annotation

    """
    
    pairtable_path = [[] for x in path]

    for x in range(len(path)):
        pairtable_path[x] = [0 for i in range(len(path[x])+1)]
        pairtable_path[x][0] = len(path[x]) 
        for z in range(len(path[x])):
            if path[x] == ".":
                pairtable_path[x][z] = 0
            o = 1
            i = 0
            match = False
            if path[x][z] == "(":
                while match == False:
                    i += 1
                    if path[x][z+i] == "(":
                        o += 1
                    elif path[x][z+i] == ")" and o != 0:
                        o -= 1

                    if path[x][z+i] == ")" and o == 0:
                        
                        pairtable_path[x][z+1] = z + i + 1
                        pairtable_path[x][z+i+1] = z + 1
                        match = True
                    
    return pairtable_path

def pairtable_to_path(pairtable):
    """Gives out a dot-bracket path 

        Args:
            pairtable (list): list of lists of a pairtable resembling a folding path
        
        Returns:
            path (list): sublits correspond to the given path but in pairtable format
    """
    path = []
    for x in range(len(pairtable)):
        struct = []
        
        for i in range(1,pairtable[x][0]+1):
            if pairtable[x][i] > i:
                struct.append("(")
            elif pairtable[x][i] < i and pairtable[x][i] != 0:
                struct.append(")")
            elif pairtable[x][i] == 0:
                struct.append(".")
        
        path.append("".join(struct))

    return path

def pairtable_to_struct(pairtable):
    """ Given a pairtable as input returns the structure in dot-bracket annotation

        Args:
            pairtable (list): pairtable
        Returns:
            struct (list): dot-bracket annotated structure in a list

    """

    struct = []
    for i in range(1,pairtable[0]+1):
        if pairtable[i] > i:
            struct.append("(")
        elif pairtable[i] < i and pairtable[i] != 0:
            struct.append(")")
        elif pairtable[i] == 0:
            struct.append(".")
    for x in range(1,len(struct)-1):
        if struct[x] == ".":
            return    
    struct = ("".join(struct))
    
    return struct

def convert_pts_out_to_nussi_in(string):
    """Removes spaces from string
    
    """


    new_string = string.replace(" ","")
    return new_string

def convert_to_UL(string):
    """Converts sting annotated domain sequence into a upper lower case annotated string

    Args: 
        string (string): space annotated domain level sequence

    Returns: 
        UL_seq (str): Upper lower case annotated domain level sequence
    
    """


    counts = string.count("*")
    c_string = list(string)
	
    for x in range(1,len(c_string)):
        if string[x] == "*":
            c_string[x-1] = c_string[x-1].upper()
    for x in range(counts):		
        c_string.remove("*")

    c_string = "".join(c_string)
    return c_string.replace(" ","")


def module_folding_path(structures,seq):
	"""Takes output of nussinov algorithm and converts it just to the folding path of the b domain for comparison against the input path. 
		It will only consider adding the strucure after a module has fully transcribed. 

		Args: 
			Structure(list): Output of nussinov algorithm 
			seq(str): sequences which was used in the nussinov algorithm
		
		Returns: 
			module_structure(list): 
	"""
	modules = seq.split("l")
	seq = convert(seq)
	module_path = []

	lengths = [0 for x in modules]

	for x in range(len(modules)):
		for z in range(x+1):
			lengths[x] = len(modules[z]) + lengths[x] 
	t = 1 
	for x in range(len(seq)):
		if seq[x] == "l":
			lengths[t] += t 
			t += 1

	for x in range(len(lengths)):
		module_path.append(structures[lengths[x]])
	liste = [[] for x in module_path]


	for x in range(len(module_path)):
		for z in range(len(module_path[x])):
			if seq[z] == "b" or seq[z] == "B":
				liste[x].append(module_path[x][z])
		
	for x in range(len(liste)):
		liste[x] = "".join(liste[x])
	return(liste)	  

def d_length(domain):
    if domain[0] == "b" or domain[0] == "B":
        return 5
    elif domain == "l":
        return 5
    return 3
    
  
def extended_domain_path(domain_path):
    """ Extends domain level path to the corresponding nt path but with domains: abc -> aaa bbbbb ccc

    Args:
        domain_path (list): sublist correspond to  path sequence
    
    Returns: 
        full_path (list): extended path 

    """

    

    UL_domain_path = convert_to_UL(domain_path)
    full_path = []
        
            
    for z in range(len(UL_domain_path)):
        
        full_path.append(UL_domain_path[z] * d_length(UL_domain_path[z]))
            


    return "".join(full_path)

def split_ntseq_to_domainfp(nt_seq,domain_seq):
    """Takes a nucleotide sequence and splits it up in subsequences where each sequence i corresponds to the sequence at transcription step i 

        Args: 
            nt_seq (string): nucleotide sequence
            domain_seq (string: domain level sequence 

        
        Returns:   
            nt_path (list): list where each sublist corresponds to sequence at transcription step
    
    """

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


if __name__=="__main__":
    print("convert_functions")
    print(path_to_pairtablepath(['.', '()', '.()', '(())']))
    #get_module_fp_sequences("AAABBBBBCCCLLLCCCBBBBBAAALLLBBBBBBBBBB")
    print(extended_domain_path("vbulj*d*a*b*c*e*k*ltmifcbaghnslr*o*h*g*a*b*c*f*i*p*q*lzwqpifcbaghorxyly*x*r*o*h*g*a*b*c*f*i*p*q*w*z*lkecbadjlu*b*v*lblb*"))



