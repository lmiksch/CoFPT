"""
Nussinov-Jacobson python algorithm implementation
	Predicts the secondary RNA structure from an RNA sequence.
	The minimal loop length is set to a default of 0
	Argument: String: takes a string with * notation as input 
"""

import numpy as np
import itertools
import pandas as pd
from functions import convert_functions as cv 

def couple(pair):
	if pair[0].upper() == pair[1].upper() and pair[0] != pair[1]:
		return True
        
	return False

def fill(nm, rna,dir_matrix):
	"""
	Fill the matrix based on the rules of the Nussinov algorithm
	"""
	
	minimal_loop_length = 0

	


	for k in range(1, len(rna)):
		for i in range(len(rna) - k):
			j = i + k
			dir = []
			if j - i >= minimal_loop_length:
				down = nm[i + 1][j] # 1st rule
				left = nm[i][j - 1] # 2nd rule
				diag = nm[i + 1][j - 1] + couple((rna[i], rna[j])) # 3rd rule
				rc = max([nm[i][t] + nm[t + 1][j] for t in range(i, j)])  # 4th rule
				max_rc = 0
				for t in range(i,j):
					rc_score = nm[i][t] + nm[t+1][j]
					if rc_score >= max_rc:
						max_rc = rc_score
						pointer_t = t
				
					
				max_score = max(down,left,diag,rc)

				if down == max_score and max_score != 0:
					dir_matrix[i][j].append("d")
				elif left == max_score and max_score != 0:
					dir_matrix[i][j].append("l")
				elif diag == max_score and max_score != 0 and couple((rna[i],rna[j])):
					dir_matrix[i][j].append("di")
				elif rc == max_score and max_score != 0:
					
					dir_matrix[i][j].append(pointer_t)

	
				



				nm[i][j] = max(down, left, diag, rc) # max of all


			
			else:
				nm[i][j] = 0
	
	return nm, dir_matrix	


def dot_write(rna, fold,x):
	dot = ["." for i in range(x)]
	
	for s in fold:
		dot[min(s)] = "("
		dot[max(s)] = ")"
	
	return "".join(dot)

def init_matrix(rna):
	M = len(rna)

	# init matrix
	nm = np.empty([M, M])
	nm[:] = np.NAN

	# init diaganols to 0
	
	nm[range(M), range(M)] = 0
	nm[range(1, len(rna)), range(len(rna) - 1)] = 0

	return nm


def nussinov(rna):
	"""
	Takes a star notation rna sequence as input and puts out a list of lists where each sublist corresponds to the possible stuctures calculated by to nussinov algorithm for the given sequence
	"""
	
	

	split_seq = rna.split()
	
	UL_liste = cv.UL_list(split_seq)
	
	nm = init_matrix(UL_liste)
	M = len(UL_liste)
	dir_matrix = [[] for x in range(M)]
	for x in range(M):
		dir_matrix[x] = [[] for x in range(M)]  

		
	nm, dir_matrix = fill(nm, UL_liste,dir_matrix)


	structures = [[] for x in range(M+1)]

	
	module_indices = []
		
	for x in range(len(UL_liste)):
		if UL_liste[x] == "l":
			module_indices.append(x+1)
	module_indices.append(len(UL_liste))

	
	structures = []	
	
	for x in module_indices:
		
		fold = []
		
		newm = nm[:x,:x]
		
	
		new_dir_matrix = [[] for x in range(M)]
		for z in range(len(dir_matrix)):
			new_dir_matrix[z] = dir_matrix[z][:x]
		
		for z in range(len(dir_matrix)-x):
			new_dir_matrix.pop(-1)
		

		folds = [[]]
		global f_pointer
		f_pointer = 0
		global check_rc
		check_rc = 0
		global number_pairs
		number_pairs = nm[0][x-1]

		traceback(new_dir_matrix, UL_liste, fold, 0, x -1,folds,newm)
		
		final_folds = folds
		
		result = []
		for sublist in final_folds:
			if len(set([i for j in sublist for i in j])) == len([i for j in sublist for i in j]):
				result.append(sublist)
		
	
		result = [x for x in result if len(x) == number_pairs]
		
		
		cur_struct = []
		for fold in result:
			
			res = dot_write(rna, fold,x)
			
			cur_struct.append(res)
		
	
		cur_struct = list(dict.fromkeys(cur_struct))
		structures.append(cur_struct)


	return(structures)

def traceback(dir_matrix, rna, fold, i, L,folds,nm):
	"""
	Traceback through matrix and explores all possible highest scoring graphs
	"""
	
	j = L 
	global f_pointer
	global check_rc 
	global number_pairs
	
	if i < j: 
		for x in range(len(dir_matrix[i][j])):
			
			if dir_matrix[i][j][x] == "d":
				traceback(dir_matrix,rna, fold, i + 1,j,folds,nm)
			elif dir_matrix[i][j][x] == "l":
				traceback(dir_matrix,rna,fold,i,j-1,folds,nm)
			elif dir_matrix[i][j][x] == "di":
				folds[f_pointer].append((i,j))
				folds[f_pointer] = list(dict.fromkeys(folds[f_pointer]))
				if len(dir_matrix[i+1][j-1]) == 0 and check_rc == 0:
					folds.append([])
					f_pointer += 1


				elif len(dir_matrix[i+1][j-1]) == 0 and len(folds[f_pointer]) == number_pairs:
					folds.append([])
					f_pointer += 1
					check_rc = 0
					



				traceback(dir_matrix, rna, fold, i + 1, j-1,folds,nm)
			else:
				for k in range(i+1,j-1):
					if type(dir_matrix[i][j][x]) == int:
						check_rc += 1
						traceback(dir_matrix,rna,fold, i, dir_matrix[i][j][x],folds, nm)
						traceback(dir_matrix,rna,fold, dir_matrix[i][j][x] + 1, j,folds,nm)
						break
	
	
	
	return folds



def find_possible_structs(structures):
	
	
	"""
	Takes a list of dot-bracket notations and finds all possible paths
	"""
	outlist =[]; templist =[[]]
	for sublist in structures:
		outlist = templist; templist = [[]]
		for sitem in sublist:
			for oitem in outlist:
				newitem = [oitem]
				if newitem == [[]]: newitem = [sitem]
				else: newitem = [newitem[0], sitem]
				templist.append(flatten(newitem))

	outlist = list(filter(lambda x: len(x)==len(structures), templist))  
	return(outlist)
	

def flatten(B):    
    A = []
    for i in B:
        if type(i) == list: A.extend(i)
        else: A.append(i)
    return A


def convert(string):
	"""converts string in star annotation to UL_seq and keeps the spaces


	"""
	counts = string.count("*")
	c_string = list(string)
	
	for x in range(1,len(c_string)):
		if string[x] == "*":
			c_string[x-1] = c_string[x-1].upper()
	for x in range(counts):		
		c_string.remove("*")

	c_string = "".join(c_string)
	return c_string


def  convert_nussi_output(output,seq):
	struct = []

	for x in range(len(seq)):
		if seq[x] == "b":
			struct.append(output[x])
	struct = "".join(struct)
	
	return struct


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
		

def remove_non_modules(structures,seq):
	
	complete_module_path = []
	module_indices = []
	seq = convert(seq)
	
	for x in range(len(seq)):
		if seq[x] == "l":
			module_indices.append(x+1)

	for x in module_indices:
		complete_module_path.append(structures[x])

	complete_module_path.append(structures[-1])
	
	return complete_module_path

def get_domain_folds(structures,rna):
	result = [[[] for z in range(len(structures[x]))] for x in range(len(structures))]
	
	for x in range(len(structures)):
		for i in range(len(structures[x])):
			for z in range(len(structures[x][i])):
				if rna[z] == "b" or rna[z] == "B":
					result[x][i].append(structures[x][i][z])
			result[x][i] = "".join(result[x][i])
	


	return result




def nussinov_modules(rna):
	""" Takes rna sequence and applies our nussinov algorithm to it. 

	Args: module_fol
		rna (str): domain level sequence with spaces

	Returns: 
		possible_paths (list): List of all possible paths were each sublist corresponds to one path. 
	
	"""
	

	nussi_output = nussinov(rna)
	print("nussi output:",nussi_output)
	
	UL_liste = cv.UL_list(rna.split())
	domain_folds = get_domain_folds(nussi_output,UL_liste)
	
	possible_paths = find_possible_structs(domain_folds)

	return possible_paths, nussi_output





if __name__ == "__main__":
	# Test
	print("nussi output:",nussinov_modules(" b   l  d* aa* b* c* e*  l  i f c  b  aa g h  l  h* g* aa* b* c* f* i*  l  j e c  b  aa d k  l  k* d* aa* b* c* e* j*"))
	


