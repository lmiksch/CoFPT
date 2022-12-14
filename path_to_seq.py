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
    

def fold_to_seq(path):
    """ Main function of fold_to_seq.

    Agrs:
        path (list): A list of the desired folding path using the pairtable annotation 

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
    for x in range(1, path[-1][0] + 1):
        
        globals()["block%s" %  x] = module("b","","",False)
        
        
        if x%2 != 0:
            globals()["block%s" % x].star = True


        seq.append(globals()["block%s" % x])

    for x in range(1,len(path)): 
        for y in range(1,len(path[x])-1):
                z = False
                try:
                    if path[x-1][path[x][y]] == 0:
                        z = True
                except:
                    pass
                if path[x][y] > path[x-1][y] and path[x-1][y] > 0 and z == False and allready_paired(x,y,path[x][y],path) == False:
                    globals()["block%s" % y].praefix  = domains[domain_pointer] + " " + globals()["block%s"  % y].praefix
                    domain_pointer += 1
                    globals()["block%s" % y].suffix = globals()["block%s" % y].suffix+ " " + domains[domain_pointer]


                    globals()["block%s" % path[x][y]].suffix =  globals()["block%s"  % y].praefix[::-1]
                    globals()["block%s" % path[x][y]].praefix = globals()["block%s"  % y].suffix[::-1]
                    domain_pointer += 1
                   
                if path[x][y+1] == 0:
                    t = y -1
                    globals()["block%s" % t].praefix  = domains[domain_pointer] + " " + globals()["block%s"  % t].praefix
                    domain_pointer += 1
                    globals()["block%s" % t].suffix = globals()["block%s" % t].suffix+ " " + domains[domain_pointer]


                    globals()["block%s" % path[x][t]].suffix = globals()["block%s" % t].praefix[::-1]
                    globals()["block%s" % path[x][t]].praefix = globals()["block%s" % t].suffix[::-1]
                    domain_pointer += 1
                
    #part which adds * notation for complementary

    for x in range(1,len(seq)+1):
        if globals()["block%s" %x].star == False: 
            #globals()["block%s" %x] = "*".join(re.findall("",globals()["block%s" %x].praefix))
            globals()["block%s" %x].praefix = stars(globals()["block%s" %x].praefix)
            globals()["block%s" %x].suffix = stars(globals()["block%s" %x].suffix)
            globals()["block%s" %x].name = stars(globals()["block%s" %x].name)  

    #part which adds extenstion block l 

    comb_modules = []
    for x in range(1,len(seq)+1):
        comb_module = globals()["block%s" %x].praefix + " " + globals()["block%s" %x].name + " " + globals()["block%s" %x].suffix
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
                print("paired")
                return True
        except:
            pass         
    return False



if __name__=="__main__":
    print("Imported path_to_seq")


    path2 = [[0, 0], [2, 2, 1], [3, 0, 3, 2], [4, 2, 1, 4, 3], [5, 0, 3, 2, 5, 4]]
    path = [[1,0],[2,2,1],[3,2,1,0],[4,2,1,4,3]]

    
    print(fold_to_seq(path2))
        
   
#print(fold_to_seq(path2))
