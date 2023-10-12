


class module():
    def __init__(self,name,praefix ="",suffix ="",star= False,highest_weight = 1):
        self.name = name
        self.praefix = praefix
        self.middle = ""
        self.suffix = suffix
        self.star = star
        self.highest_weight = highest_weight
        

    def __str__(self):
            return f"{self.praefix}{self.name}{self.suffix}"
    
    # dont know if necessary
    def set_middle_domain(self,middle_domain): 
         self.middle = middle_domain
    

class graph():

    def __init__(self):
         self.graph = {}
    
    def add_node(self,node):
         if node not in self.graph:
              self.graph[node] = []

    def add_edge(self,node1,node2):
        if node1 in self.graph:
            self.graph[node1].append(node2)
        else:
            self.graph[node1] = [node2]

        if node2 in self.graph:
            self.graph[node2].append(node1)
        else:
            self.graph[node2] = [node1]

    def create_nodes_from_pairtable(self,pairtable):
        for x in range(len(1,pairtable[-1][0])):
            node = module(x)
            self.add_node(node)

    def dfs(self, start_node):
        visited = set()

        def dfs_recursive(node):
            visited.add(node.name)
            print(node.name)

            for neighbor in node.neighbors:
                if neighbor.data not in visited:
                    dfs_recursive(neighbor)

        dfs_recursive(start_node)
        
"""Modules in Graph as nodes first initialize 


Args:
    afp(list): each entry in the list corresponds to one step in the abstract folding path the input can either be dot-bracket annotated or as a pairtable should be in list 


"""

afp = [[0, 0], [2, 2, 1], [3, 2, 1, 0], [4, 4, 3, 2, 1]]

#create a module/node for each module in the folding path present at the last step and append it to the graph


