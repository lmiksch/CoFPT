# -*- coding: utf-8 -*-
"""
Takes a folding path in form a list of pairtables as input and puts out a domain level sequence which has the same folding path as the input
Folding path should be at any times satured and  free from pseudo-knots

Example: 
input: [[1,0],[2,2,1],[3,2,1,0],[4,2,1,4,3]]
output: [[.],[()],[().],[()()]]


"""
import  string
class module():
    def __init__(self,name,praefix,suffix,star):
        self.name = name
        self.praefix = praefix
        self.suffix = suffix
        self.star = star


    def __str__(self):
            return f"{self.praefix}{self.name}{self.suffix}"

def pairtable_to_path(pairtable):
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
    

def fold_to_seq(fp_length):
    """ Main function of fold_to_seq.globals()["block%s" % x].suffix = ""n 

    Return:
        new_string (str): A string of domains which should follow the given folding path. 

    """
    
    domains = list(string.ascii_lowercase)
    for x in range(len(string.ascii_lowercase)):
        for  z in range(len(string.ascii_lowercase)):
            domains.append(string.ascii_lowercase[x] + string.ascii_lowercase[z])
    domains.remove("b")
    domains.remove("l")
    domain_pointer = 0   
    #creates a list with the corresponding number of blocks  used in the folding path 
    seq = []
    used_domains = ["b"]
    path = [[] for x in range(fp_length)]
    path[-1].append(fp_length)
    for x in range(1, path[-1][0] + 1):
        
        globals()["block%s" %  x] = module("b","","",False)
        globals()["block%s" % x].praefix = ""
        globals()["block%s" % x].suffix = ""
        
        if x%2 != 0:
            globals()["block%s" % x].star = True


        seq.append(globals()["block%s" % x])


    for x in range(3,fp_length):

        #second to last block \
        #print(type(globals()["block%s" % int(x-1)].suffix)) 
        #print("suffix",globals()["block%s" % int(x-1)].suffix)


        #print(type(globals()["block%s" % int(x-1)].praefix)) 
        #print(globals()["block%s" % int(x-1)].praefix)

        if len(globals()["block%s" % int(x-1)].suffix ) == 0:
            globals()["block%s" % int(x-1)].suffix = domains[domain_pointer]

        else:
            globals()["block%s" % int(x-1)].suffix = globals()["block%s" % int(x-1)].suffix + ' ' +  domains[domain_pointer]

        #print("globals suffix", globals()["block%s" % int(x-1)].suffix)

        #last block
        suffix_split = globals()["block%s"  % int(x-1)].suffix.split()
        suffix_rev = " ".join(suffix_split[::-1])


        praefix_split = globals()["block%s"  % int(x-1)].praefix.split()
        praefix_rev = " ".join(praefix_split[::-1])
        
        globals()["block%s" % x].praefix = suffix_rev
        #globals()["block%s" % x].suffix = praefix_split
        
        globals()["block%s" % x].suffix = " ".join(praefix_rev)

        
        domain_pointer += 1

    #part which adds * notation for complementary

    for x in range(1,len(seq)+1):
        if globals()["block%s" %x].star == False: 
            
            globals()["block%s" %x].praefix = stars(globals()["block%s" %x].praefix)
            globals()["block%s" %x].suffix = stars(globals()["block%s" %x].suffix)
            globals()["block%s" %x].name = stars(globals()["block%s" %x].name)  

    #part which adds extenstion block l 

    comb_modules = []
    for x in range(1,len(seq)+1):
        comb_module = str(globals()["block%s" %x].praefix) + " " + str(globals()["block%s" %x].name) + " " + str(globals()["block%s" %x].suffix)
        comb_modules.append(comb_module)

    
    final_seq = "  l  ".join(comb_modules)
    

    #print("done")
    return final_seq

def stars(string):
    string = string.split()
    for x in range(len(string)):
        string[x] = string[x] + "*"
  
    new_string = ' '.join(string)
    return new_string

def allready_paired(step,x,y,path):
    for i in range(1,step):
        try:
            if path[i][x] == y:
                return True
        except:
            pass         
    return False



if __name__== "__main__":
    print("Imported path_to_seq")


    path = [[0, 0], [2, 2, 1], [3, 0, 3, 2], [4, 2, 1, 4, 3], [5, 0,5,4,3,2],[6,2,1,4,3,6,5]]
    
    #Out:
    #b l0 b* a* l1 a b c l2 c* b* a* d* l3 d a b c e l4 e* c* b* d* a* f* l5
    #b l0 b* a* l1 a b c l2 c* b* a* d* l3 d a b c e l4 e* c* b* a* d* 
    
    print(fold_to_seq(12))
        
   
