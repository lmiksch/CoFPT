"""
Nussinov-Jacobson python algorithm implementation
	Predicts the secondary RNA structure from an RNA sequence.
	The minimal loop length is set to a default of 0
"""

from ast import expr_context
import numpy as np
import pandas as pd # you can remove this... it's just for nice output (matrix formatting)

def couple(pair):
	
    if pair[0].upper() == pair[1].upper() and pair[0] != pair[1]:
        return True
        
    return False

def fill(nm, rna):
	"""
	Fill the matrix as per the Nussinov algorithm
	"""
	minimal_loop_length = 0

	for k in range(1, len(rna)):
		for i in range(len(rna) - k):
			j = i + k

			if j - i >= minimal_loop_length:
				down = nm[i + 1][j] # 1st rule
				left = nm[i][j - 1] # 2nd rule
				diag = nm[i + 1][j - 1] + couple((rna[i], rna[j])) # 3rd rule

				rc = max([nm[i][t] + nm[t + 1][j] for t in range(i, j)])  # 4th rule

				nm[i][j] = max(down, left, diag, rc) # max of all
			
			else:
				nm[i][j] = 0

	return nm	

def traceback(nm, rna, fold, i, L):
	"""
	Traceback through complete Nussinov matrix to find optimial RNA secondary structure solution through max base-pairs
	"""
	j = L
	if i < j:
		if nm[i][j] == nm[i + 1][j]: # 1st rule
			traceback(nm, rna, fold, i + 1, j)
		elif nm[i][j] == nm[i][j - 1]: # 2nd rule
			traceback(nm, rna, fold, i, j - 1)
		elif nm[i][j] == nm[i + 1][j - 1] + couple((rna[i], rna[j])): # 3rd rule
			fold.append((i, j))
			traceback(nm, rna, fold, i + 1, j - 1)
		else:
			for k in range(i + 1, j - 1):
				if nm[i][j] == nm[i, k] + nm[k + 1][j]: # 4th rule
					traceback(nm, rna, fold, i, k)
					traceback(nm, rna, fold, k + 1, j)
					break

	return fold

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
	# few ways to do this: np.fill_diaganol(), np.diag(), nested loop, ...
	nm[range(M), range(M)] = 0
	nm[range(1, len(rna)), range(len(rna) - 1)] = 0

	return nm




def nussinov(rna):
	rna = convert(rna)
	nm = init_matrix(rna)
	nm = fill(nm, rna)
	print(nm)
	structures = []
	for x in range(0,len(rna)+1):
		fold = []
		newm = nm[:x,:x]
		sec = traceback(newm, rna, fold, 0, x -1)
		res = dot_write(rna, fold,x)
		structures.append(res)
	return(structures)
def convert(string):
	counts = string.count("*")
	c_string = list(string)
	
	for x in range(1,len(c_string)):
		if string[x] == "*":
			c_string[x-1] = c_string[x-1].upper()
	for x in range(counts):		
		c_string.remove("*")

	c_string = "".join(c_string)
	return c_string
    

rna = "AaAaAaAaaaaaAAAA"

string = "abclb*lblc*b*a*"

test = "bla*b*c*ldcbaele*a*b*c*d*"







for struct in nussinov(string):
	print(struct)
	print(convert(test))
