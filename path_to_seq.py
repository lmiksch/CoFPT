"""
Takes a folding path in form a list of pertables as input and puts out a domain level sequence which has the same folding path as the input
Folding path should be at any times satured and  free from pseudo-knots
"""
import  string
class domain():
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
    domains = list(string.ascii_lowercase)
    domains.remove("b")
    domains.remove("l")
    domain_pointer = 0
   
    #creates a list with the corresponding number of blocks  used in the folding path 
    seq = []
    used_domains = ["b"]
    for x in range(1, path[-1][0] + 1):
        
        globals()["block%s" %  x] = domain("b","","","F")
        
        
        if x%2 != 0:
            globals()["block%s" % x].star = "T"


        seq.append(globals()["block%s" % x])

    

    for x in range(1,len(path)): 
        
        for y in range(1,len(path[x])-1):
                z = False
                try:
                    if path[x-1][path[x][y]] == 0:
                        z = True
                except:
                    pass
                if path[x][y] > path[x-1][y] and path[x-1][y] > 0 and z == False:
                    globals()["block%s" % y].praefix = domains[domain_pointer] + globals()["block%s" % y].praefix
                    domain_pointer += 1
                    globals()["block%s" % y].suffix = globals()["block%s" % y].suffix + domains[domain_pointer]

                    globals()["block%s" % path[x][y]].suffix = globals()["block%s" % y].praefix[::-1]
                    globals()["block%s" % path[x][y]].praefix = globals()["block%s" % y].suffix[::-1]
                    domain_pointer += 1

    #part which adds * notation for complementary

    for x in range(1,len(seq)+1):
        
        if x % 2 == 0: 
            #globals()["block%s" %x] = "*".join(re.findall("",globals()["block%s" %x].praefix))

            globals()["block%s" %x].praefix = stars(globals()["block%s" %x].praefix)
            globals()["block%s" %x].suffix = stars(globals()["block%s" %x].suffix)
            globals()["block%s" %x].name = stars(globals()["block%s" %x].name)  

    #part which adds extenstion domain l 

    modules = []
    for x in range(1,len(seq)+1):
        module = "".join(globals()["block%s" %x].praefix + globals()["block%s" %x].name + globals()["block%s" %x].suffix)
        modules.append(module)

    #print(modules)
    final_seq = "l".join(modules)
    

    #print("done")
    return final_seq

def stars(string):
    
    new_string = '*'.join(string[i:i+1] for i in range(0, len(string)+1, 1))
    return new_string

path2 = [[1,0],[2,2,1],[3,0,3,2],[4,4,3,2,1],[5,2,1,4,3,0],[6,2,1,6,5,4,3]]
path = [[1,0],[2,2,1],[3,0,3,2],[4,2,1,4,3]]

#print(fold_to_seq(path))
        


#print(fold_to_seq(path2))