"""
This module creates all possible folding paths given a domain sequence. This is mainly used to test our path_to_seq algortihm.
"""
import random

import math



def couple(pair):
    if pair[0].upper() == pair[1].upper() and pair[0] != pair[1]:
        return True

def find_structs(seq):
    structs = []
    possible_pairs = [[] for x in range(len(seq))]
    all_found  = False
    for x in range(len(seq)):
            for y in range(len(seq)):
                if couple((seq[x],seq[y])):
                    possible_pairs[x].append(y)

    allpairs = find_all_pairs(possible_pairs)

    #print("allpairs:", allpairs)
    paired = [-1 for x in range(len(seq))]

    global final
    final = []
    find_all_combinations(allpairs,0,[],paired)
    #print("final",final)

    p_list_final = []
    for x in final:
        p_list_final.append(translate(x,seq))

    p_list_final = [list(i) for i in set(map(tuple, p_list_final))]
    #print(p_list_final)

    
    final_structures = []

    for x in p_list_final:
        #print("x = :",x)
        final_structures.append(pertable_to_struct(x))
    
    #print(final_structures)

    return final_structures


def find_all_combinations(pairs,x,struct,paired):
    
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

def pertable_to_struct(pertable):
    struct = []
    #print(pertable)
    for i in range(1,pertable[0]+1):
        if pertable[i] > i:
            struct.append("(")
        elif pertable[i] < i and pertable[i] != 0:
            struct.append(")")
        elif pertable[i] == 0:
            struct.append(".")
        
    struct = ("".join(struct))
    
    return struct



def generate_foldingpaths(seq):
    folding_paths = []
    for x in range(0,len(seq)+1):
        folding_paths.append(find_structs(seq[:x]))

    return folding_paths

def get_random_path(folding_paths):
    path = []
    for x in range(len(folding_paths)):
        print("x",x)
        path.append(folding_paths[x][random.randint(0,len(folding_paths[x])-1)])

    return(path)


seq = "aAaAaAaAaAaA"


paths = generate_foldingpaths(seq)

print("paths")
print(len(paths))

for x in paths:
    print(x)


print(get_random_path(paths))


#find_structs("aAaAaAaAaAaA")



