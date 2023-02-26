# -*- coding: utf-8 -*-
"""
This module creates all possible folding paths given a domain sequence. This is mainly used to test our path_to_seq algortihm.
"""
import random

import math



def couple(pair):
    if pair[0].upper() == pair[1].upper() and pair[0] != pair[1]:
        return True

def find_structs(seq):
    """Finds all possible structures given a sequence

    Args: 
        seq (str): UL_seq
    Returns: 
        final_structures: all possible structures given the sequence
    """
    structs = []
    possible_pairs = [[] for x in range(len(seq))]
    all_found  = False
    for x in range(len(seq)):
            for y in range(len(seq)):
                if couple((seq[x],seq[y])):
                    possible_pairs[x].append(y)
    
    allpairs = find_all_pairs(possible_pairs)
   
    paired = [-1 for x in range(len(seq))]

    global final
    final = []
    find_all_combinations(allpairs,0,[],paired)
    





    p_list_final = []
    for x in final:
        p_list_final.append(translate(x,seq))

    p_list_final = [list(i) for i in set(map(tuple, p_list_final))]
    
    
    final_structures = []

    for x in p_list_final:
       
        final_structures.append(pairtable_to_struct(x))
   

    final_structures = [i for i in final_structures if i is not None]
    return final_structures


def find_all_combinations(pairs,x,struct,paired):
    """ Finds all possible folding paths given structures at different transcription lengths
    """
    
    if len(struct) >= math.sqrt(len(pairs)) - 0.5:
        final.append(list(struct))
        
        return 
    
    for i in range(len(pairs)):

        if pairs[i][0] not in paired and pairs[i][1] not in paired:
            struct.append(pairs[i])
            paired.append(pairs[i][0])
            paired.append(pairs[i][1])
            find_all_combinations(pairs,x+1,struct,paired)

            paired.pop()
            paired.pop()
            struct.pop()


def find_all_pairs(possible_pairs):
        pairs = []
        for x in range(len(possible_pairs)):
            for y in range(len(possible_pairs[x])):
                pairs.append((x,possible_pairs[x][y]))
        L = []
        for ll in pairs:
            if ll[::-1] not in L:
                L.append(ll)
        return L

def translate(t_list,seq):
    p_table = [0 for x in range(len(seq))]
    for x in range(len(t_list)):
            p_table[t_list[x][0]] = t_list[x][1] + 1
            p_table[t_list[x][1]] = t_list[x][0] + 1
    p_table.insert(0,len(seq))
    return p_table

def pairtable_to_struct(pairtable):
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



def generate_foldingpaths(seq):
    folding_paths = []
    for x in range(0,len(seq)+1):
        folding_paths.append(find_structs(seq[:x]))
    

    return folding_paths

def get_random_path(folding_paths):
    path = []
    for x in range(1,len(folding_paths)):
        path.append(folding_paths[x][random.randint(0,len(folding_paths[x])-1)])

    return(path)

if __name__ == "__main__":
    print("Hello")

