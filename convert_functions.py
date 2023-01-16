# -*- coding: utf-8 -*-
"""
This module contains all the different convert functions which can be used to convert outputs into different annotations. 

"""


def path_to_pertablepath(path):
    """ Function which turns a structurral path in dot bracket annotation into a pertable path

    args: 
        path(list): List of the folding path in dot bracket
    output: 
        pertable_path(List): List of the folding path using the pertable annotation

    """
    print("Input:", path)
    pertable_path = [[] for x in path]

    for x in range(len(path)):
        pertable_path[x] = [0 for i in range(len(path[x])+1)]
        pertable_path[x][0] = len(path[x]) 
        for z in range(len(path[x])):
            if path[x] == ".":
                pertable_path[x][z] = 0
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
                        
                        pertable_path[x][z+1] = z + i + 1
                        pertable_path[x][z+i+1] = z + 1
                        match = True
                    
    return pertable_path

def pertable_to_path(pertable):
    path = []
    for x in range(len(pertable)):
        struct = []
        
        for i in range(1,pertable[x][0]+1):
            if pertable[x][i] > i:
                struct.append("(")
            elif pertable[x][i] < i and pertable[x][i] != 0:
                struct.append(")")
            elif pertable[x][i] == 0:
                struct.append(".")
        
        path.append("".join(struct))

    return path

def pertable_to_struct(pertable):
    struct = []
    for i in range(1,pertable[0]+1):
        if pertable[i] > i:
            struct.append("(")
        elif pertable[i] < i and pertable[i] != 0:
            struct.append(")")
        elif pertable[i] == 0:
            struct.append(".")
    for x in range(1,len(struct)-1):
        if struct[x] == ".":
            return    
    struct = ("".join(struct))
    
    return struct

def convert_pts_out_to_nussi_in(string):
    new_string = string.replace(" ","")
    return new_string

def convert_to_UL(string):
	counts = string.count("*")
	c_string = list(string)
	
	for x in range(1,len(c_string)):
		if string[x] == "*":
			c_string[x-1] = c_string[x-1].upper()
	for x in range(counts):		
		c_string.remove("*")

	c_string = "".join(c_string)
	return c_string.replace(" ","")

def get_module_fp_sequences(nt_seq,domain_seq):
    """
    Input NT sequence
    
    Puts out a list of the NT corresponding to the module folding path"""

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





if __name__=="__main__":
    print("convert_functions")
    print(path_to_pertablepath(['.', '()', '.()', '(())']))
    get_module_fp_sequences("AAABBBBBCCCLLLCCCBBBBBAAALLLBBBBBBBBBB")



