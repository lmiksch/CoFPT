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



if __name__=="__main__":
    print("convert_functions")
    print(path_to_pertablepath(['.', '()', '.()', '(())']))



