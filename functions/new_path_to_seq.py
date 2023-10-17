import convert_functions as cv 
from collections import deque
import inequality_solver

class node():
    def __init__(self,name,praefix ="",suffix ="",complement= False,highest_weight = 1,neighbors = None):
        self.name = name
        self.praefix = praefix
        self.middle = ""
        self.suffix = suffix
        self.complement = complement
        self.highest_weight = highest_weight
        self.neighbors = set() if neighbors is None else neighbors
        

    def __str__(self):
        neighbor_names = ', '.join(str(neighbor.name) for neighbor in self.neighbors)
        return f"Name: {self.name} Sequence: {self.praefix}{self.middle}{self.suffix}, neighbors: {neighbor_names}, Complement: {self.complement}"

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
        if node2 not in node1.neighbors:
            node1.neighbors.add(node2)  
            node2.neighbors.add(node1)  


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
            for neighbor in node.neighbors:
                self.edges.add(str(node.name) + str(neighbor.name))


        #Remove duplicate edges
        self.edges = {_ for _ in self.edges if int(_[1]) >= int(_[0])} 

    
    
    def bipartite_check(self,connected_components):
        

        visited_nodes = {}
        stack = [node for node in self.graph if node.name in connected_components]
        #print("stack",stack)
        complement = True
        while stack:
            current_node = stack.pop()
            
            print("\nnodes in stack:")
            for node in stack:
                print(node)
            print(" ")
        
            #complement = False  # Initialize as False, indicating it's in the first set
            if current_node not in visited_nodes:
                visited_nodes[current_node] = True
                current_node.complement = complement
            
                complement = not complement 

            for neighbor in current_node.neighbors:
                if neighbor not in visited_nodes:
                    neighbor.complement = complement
                    stack.append(neighbor)
                    #visited_nodes[neighbor] = True
                print("\ncurrent node:", current_node)
                print("neighbor:",neighbor)
                
                # Check if the neighbor's complement is the same as the current node
                if current_node in visited_nodes:
                    if neighbor.complement == current_node.complement:
                        print("\nhello there")
                        self.print_nodes()
                        return False  # Graph is not bipartite
                
                neighbor.complement = not current_node.complement

        return True 
    
    def get_current_edges(self,step_number):
        current_edges = []
        for edge in self.edges:
            if int(edge[1]) <= step_number + 1:
                current_edges.append(edge)
        return current_edges


    def get_inequalities(self):
        
        inequalities = []

        for x,step in enumerate(afp): 
            collected_edges = self.get_current_edges(x)
            collected_edges = sorted(collected_edges, key=lambda x: -(int(x) // 10 % 10))
            print("\nCollected Edges",collected_edges)
            print("Current Step:",step)
            active_edges = []
            inactive_edges = []
            
            for edge in collected_edges:
                if step[int(edge[0])] == int(edge[1]) and step[int(edge[1])] == int(edge[0]):
                    active_edges.append(edge)
                else:
                    inactive_edges.append(edge)
            print("Active Edges",active_edges)
            print("Inactive Edges",inactive_edges)  


            #Remove inactive edges which were allready defined
            if len(inequalities) > 0:    
                for ineq in inequalities:
                    result = [[act, inact] for inact in inactive_edges for act in active_edges if ineq[0] == act and ineq[1] == inact]
                    print(result)
                    if result:
                        print(f"Remove following edge {result[0][1]} du to result beeing {result}")
                        inactive_edges.remove(result[0][1])

            while len(inactive_edges) != 0: 
                print("\nBegin ineq inactive edges:",inactive_edges)
                current_edge = active_edges.pop()
                print(f"Current Edge: {current_edge}")
                l_node = current_edge[0]
                r_node = current_edge[1]
                
                

                neighbor_edge = [edge for edge in inactive_edges if edge[0] == r_node or edge[1] == r_node]

                print(f"R neighbor edge {neighbor_edge}")
                for edge in neighbor_edge:
                    inequalities.append([current_edge,edge])
                    inactive_edges.remove(edge)

                neighbor_edge = [edge for edge in inactive_edges if edge[0] == l_node or edge[1] == l_node]

                print(f"L neighbor edge: {neighbor_edge}")
                for edge in neighbor_edge:
                    inequalities.append([current_edge,edge])
                    inactive_edges.remove(edge)

                

            print(f"\nInequalities {inequalities}")

        final_inequalities = []
        print("\nResulting Inequalties: ",inequalities)
        for ineq in inequalities:
            left = ineq[0]
            right = ineq[1]
            

            final_inequalities.append(f"{left} > {right}")

        return(final_inequalities)
        
        







"""Modules in Graph as nodes first initialize 


Args:
    afp(list): each entry in the list corresponds to one step in the abstract folding path the input can either be dot-bracket annotated or as a pairtable should be in list 

    
"""

afp = [[1, 0], [2, 2, 1], [3, 0, 3, 2], [4, 4, 3, 2, 1]]

afp = [[1,0],[2,2,1],[3,0,3,2],[4,4,3,2,1],[5,0,3,2,5,4],[6,6,3,2,5,4,1]]
#afp = [[1,0],[2,2,1],[3,0,3,2],[4,4,3,2,1],[5,5,0,4,3,1]]
afp = [".","()",".()","(())","(.())","()()()","()(.())"]
#create a module/node for each module in the folding path present at the last step and append it to the graph

def build_graph(afp):

    if afp[0][0] != int:
        afp = cv.path_to_pairtablepath(afp)

    print(afp)


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
            if step[x] != 0:
                node1 = nodes[x]
                node2 = nodes[step[x]]
                
                afp_graph.add_edge(node1, node2)
            else:
                continue

    # Print the graph after adding edges
    print("Graph after adding edges:")
    afp_graph.print_nodes()
    afp_graph.get_edges()
    
    print("Edges: ", afp_graph.edges)

    connected_components = cv.find_connected_modules(afp)
    print("\nFollowing Nodes are connected:")
    for component in connected_components:
        nodes_in_component = set()
        for node_index in component:
            nodes_in_component.add(nodes[node_index])
        print(component)

        if not afp_graph.bipartite_check(component):
            raise ImportError("Graph not Bipartite can't design domain level sequence. Check your input")

    print("\nBipartite Check Complete:")
    for x in afp_graph.graph:
        print(x)

    inequalities = afp_graph.get_inequalities()

    final_edges = afp_graph.get_current_edges(len(afp))
    print("\n Final Inequalities",inequalities)
    inequality_solver.inequality_solver(inequalities,final_edges)



    
build_graph(afp)