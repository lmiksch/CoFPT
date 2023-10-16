


class node():
    def __init__(self,name,praefix ="",suffix ="",star= False,highest_weight = 1,neighbours = None):
        self.name = name
        self.praefix = praefix
        self.middle = ""
        self.suffix = suffix
        self.star = star
        self.highest_weight = highest_weight
        self.neighbours = set() if neighbours is None else neighbours
        

    def __str__(self):
        neighbor_names = ', '.join(str(neighbor.name) for neighbor in self.neighbours)
        return f"Name: {self.name} Sequence: {self.praefix}{self.middle}{self.suffix}, Neighbors: {neighbor_names}"

    # dont know if necessary
    def set_middle_domain(self,middle_domain): 
         self.middle = middle_domain
    

class graph():

    def __init__(self):
        self.graph = {}
        self.edges = set()
    
    def add_node(self,node):
         if node not in self.graph:
              self.graph[node] = []

    def add_edge(self,node1,node2):    
        if node2 not in node1.neighbours:
            node1.neighbours.add(node2)  
            node2.neighbours.add(node1)  


    def create_nodes_from_pairtable(self,pairtable):
        for x in range(1,pairtable[-1][0] + 1):
            new_node = node(x)
            self.add_node(new_node)

    def dfs(self, start_node):
        visited = set()

        def dfs_recursive(node):
            visited.add(node.name)
            print(node.name)

            for neighbor in node.neighbors:
                if neighbor.data not in visited:
                    dfs_recursive(neighbor)

        dfs_recursive(start_node)

    def print_nodes(self):
        for node in self.graph:
            print(node)
    
    def get_edges(self):
        
        for node in self.graph:
            for neighbour in node.neighbours:
                self.edges.add(str(node.name) + str(neighbour.name))


        #Remove duplicate edges
        self.edges = {_ for _ in self.edges if int(_[1]) >= int(_[0])} 
        


"""Modules in Graph as nodes first initialize 


Args:
    afp(list): each entry in the list corresponds to one step in the abstract folding path the input can either be dot-bracket annotated or as a pairtable should be in list 

    
"""

afp = [[0, 0], [2, 2, 1], [3, 2, 1, 0], [4, 4, 3, 2, 1]]

#create a module/node for each module in the folding path present at the last step and append it to the graph

def build_graph(afp):

    afp_graph = graph()
    afp_graph.create_nodes_from_pairtable(afp)

    #print(afp_graph)

    print("Graph after adding all Nodes")
    afp_graph.print_nodes()

    nodes = list(afp_graph.graph.keys())
    nodes.insert(0,0)

    print("begin \n")
    for step in afp[::-1]:
        for x in range(1,step[0]):
            if step[x] == "0":
                continue
            else:
                node1 = nodes[x]
                node2 = nodes[step[x]]
                
                afp_graph.add_edge(node1, node2)
        

    # Print the graph after adding edges
    print("Graph after adding edges:")
    afp_graph.print_nodes()

    afp_graph.get_edges()
    
    print("Edges: ", afp_graph.edges)

build_graph(afp)