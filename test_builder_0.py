# algorithm
# To build a grahene like
'''
    c  c  c  c  c  c
   / \/ \/ \/ \/ \/ \
  c  c  c  c  c  c  c
  |  |  |  |  |  |  |
  c  c  c  c  c  c  c
   \ /\ /\ /\ /\ /\ /
    c  c  c  c  c  c
    |  |  |  |  |  |
   / \/ \/ \/ \/ \/ \
  c  c  c  c  c  c  c
  |  |  |  |  |  |  |
  c  c  c  c  c  c  c
   \ /\ /\ /\ /\ /\ /
    c  c  c  c  c  c

We need to input two lines of structure as repeat
units, where we also need to decide how many layers
we need
'''
import networkx as nx
import numpy as np
import copy

def print_note(note):
    element = note["type"]
    x, y, z = note["coords"][:]
    print(f"{element} {x} {y} {z}\n")

def Change_vector_to_bondlength(x, y, bondlength):
    norm = x**2 + y**2
    norm = np.sqrt(norm)
    x = x/norm*bondlength
    y = y/norm*bondlength
    return x, y




length_C_C = 1.414
length_C_H = 1.14
G = nx.Graph()
# Create a benzene
G.add_node(0, type=1, coords=np.array([0,               length_C_C, 0]), pos=0)
G.add_node(1, type=1, coords=np.array([np.sqrt(3)/2*length_C_C, length_C_C/2, 0]), pos=0)
G.add_node(2, type=1, coords=np.array([np.sqrt(3)/2*length_C_C, -length_C_C/2, 0]), pos=0)
G.add_node(3, type=1, coords=np.array([0,               -length_C_C, 0]),  pos=0)
G.add_node(4, type=1, coords=np.array([-np.sqrt(3)/2*length_C_C, -length_C_C/2, 0]), pos=0)
G.add_node(5, type=1, coords=np.array([-np.sqrt(3)/2*length_C_C,length_C_C/2, 0]), pos=0)
# Add bond
for i in range(5):
    G.add_edge(i, i+1)
G.add_edge(5, 0)

# Trans direction add molecule

def Trans_add(trans_num, Graph):
    # trans_num should be an integer>0
    # means how many benzene in a line
    # Graph means nx.graph
    '''
    Each time add 4 more atoms at right side
    '''
    if trans_num<1:
        print("Trans direction is less than 0")
    # At boundary
    Graph.add_node(6, type=1, coords=np.array([np.sqrt(3)*length_C_C,     length_C_C, 0]), pos=0)
    Graph.add_edge(6, 1)
    Graph.add_node(7, type=1, coords=np.array([np.sqrt(3)*3/2*length_C_C, length_C_C/2, 0]), pos=0)
    Graph.add_edge(6, 7)
    Graph.add_node(8, type=1, coords=np.array([np.sqrt(3)*3/2*length_C_C, -length_C_C/2, 0]), pos=0)
    Graph.add_edge(8, 7)
    Graph.add_node(9, type=1, coords=np.array([np.sqrt(3)*length_C_C,     -length_C_C, 0]), pos=0)
    
    Graph.add_edge(9, 8)
    Graph.add_edge(9, 2)
    # In loop
    for i in range(0, trans_num-2):
        Graph.add_node(9 + 4*i+1, type=1, coords=np.array([(i+2)*np.sqrt(3)*length_C_C, length_C_C, 0]), pos=0)
        # The new ring is always 3 different
        Graph.add_edge(9 + 4*i+1, 9 + 4*i-2)
        Graph.add_node(9 + 4*i+2, type=1, coords=np.array([(5/2+i)*np.sqrt(3)*length_C_C, length_C_C/2, 0]), pos=0)
        Graph.add_edge(9 + 4*i+2, 9 + 4*i+1)
        Graph.add_node(9 + 4*i+3, type=1, coords=np.array([(5/2+i)*np.sqrt(3)*length_C_C, -length_C_C/2, 0]),
                pos=0)
        Graph.add_edge(9 + 4*i+3, 9 + 4*i+2)
        Graph.add_node(9 + 4*i+4, type=1, coords=np.array([(i+2)*np.sqrt(3)*length_C_C, -length_C_C, 0]), pos=0)
        Graph.add_edge(9 + 4*i+4, 9 + 4*i+3)
        Graph.add_edge(9 + 4*i+4, 9 + 4*i-1)

    print("The total node is: ", Graph.number_of_nodes())
    return Graph


# Add trans benzene
G=Trans_add(6, G)

# Vertically add benzene
'''
            c  c  c  c  c  c
           / \/ \/ \/ \/ \/ \
          c  c  c  c  c  c  c
          |  |  |  |  |  |  |
          c  c  c  c  c  c  c
     /     \ /\ /\ /\ /\ /\ /        \
    C       c  c  c  c  c  c          C
    |  +    |  |  |  |  |  |     +    |
    C      / \/ \/ \/ \/ \/ \         C
     \    c  c  c  c  c  c  c        /
          |  |  |  |  |  |  |
          c  c  c  c  c  c  c
           \ /\ /\ /\ /\ /\ /
            c  c  c  c  c  c

'''


def Vert_add(vert_add, G):
    layer = int((vert_add-1)/2)
    print("layer is :", layer)
    num_benzene = 0
    num_benzene = int((G.number_of_nodes()-2)/4)
    mapping = {}
    for i in range(layer):
        print("i is ", i)
        Graph = copy.deepcopy(G)
        for node in Graph:
            Graph.nodes[node]["coords"][1] += 3*length_C_C*layer
        print("Graph is: ", Graph.number_of_nodes())
        print("G is: ", G.number_of_nodes())
        # Save the index for later use
        line_num = copy.deepcopy(max(G.nodes))
        Graph = nx.relabel_nodes(Graph, {n: n + max(G.nodes)+1 for n in Graph.nodes})
        G = nx.compose(G, Graph)
        print("Total number: ", G.number_of_nodes())
        # How many benzene
        #G.add_edge(0+(i-1)*line_num, Graph.number_of_nodes()+3)
        G.add_edge(0, Graph.number_of_nodes()+3)
        
        for num in range(0, num_benzene-1):
            # G.add_edge(9+Graph.number_of_nodes()+4*num, 6+4*num+line_num*(i-1))
            G.add_edge(9+Graph.number_of_nodes()+4*num, 6+4*num)

    if layer%2 ==1:
        # Waiting for fix, time out here.
        print("Please change to even")

    # Add boundary atom.
    if vert_add==3:
        print("Adding boundary")
        # For the beginning
        G.add_node(G.number_of_nodes()+1, type=1, coords=np.array([G.nodes[5]["coords"][0]-np.sqrt(3)/2*length_C_C, G.nodes[5]["coords"][1]+1/2*length_C_C, 0]), pos=0)
        G.add_edge(G.number_of_nodes(), 5)
        G.add_node(G.number_of_nodes()+1, type=1, coords=np.array([G.nodes[5]["coords"][0]-np.sqrt(3)/2*length_C_C, G.nodes[5]["coords"][1]+3/2*length_C_C, 0]), pos=0)
        G.add_edge(G.number_of_nodes(), G.number_of_nodes()-1)
        G.add_edge(G.number_of_nodes(), line_num+5)
        # For the end
        G.add_node(G.number_of_nodes()+1, type=1, coords=np.array([G.nodes[line_num-2]["coords"][0]+np.sqrt(3)/2*length_C_C, G.nodes[5]["coords"][1]+1/2*length_C_C, 0]), pos=0)
        G.add_edge(G.number_of_nodes(), line_num-2)
        G.add_node(G.number_of_nodes()+1, type=1,
                coords=np.array([G.nodes[line_num-2]["coords"][0]+np.sqrt(3)/2*length_C_C, G.nodes[5]["coords"][1]+3/2*length_C_C, 0]), pos=0)
        G.add_edge(G.number_of_nodes(), G.number_of_nodes()-1)
        G.add_edge(G.number_of_nodes(), line_num*2)

    return G







G = Vert_add(3, G)






x = []
y = []
x_v = 0
y_v = 0

atom_n = G.number_of_nodes()
# Add hydrogen
Orginal = copy.deepcopy(G)
for node in Orginal:
    if len(list(Orginal.neighbors(node)))<3:
        atom_n += 1
        for neighbor in Orginal.neighbors(node):
            # print_note(Orginal.nodes[neighbor])
            x.append(Orginal.nodes[neighbor]["coords"][0]-Orginal.nodes[node]["coords"][0])
            y.append(Orginal.nodes[neighbor]["coords"][1]-Orginal.nodes[node]["coords"][1])
        print(node)
        # Vector operation
        x_v = -(x[1]+x[0])/2
        y_v = -(y[1]+y[0])/2
        x_v, y_v = Change_vector_to_bondlength(x_v, y_v, length_C_H)
        x_v = Orginal.nodes[node]["coords"][0]+x_v
        y_v = Orginal.nodes[node]["coords"][1]+y_v
        # Add Hydrogen
        G.add_node(atom_n, type=2, coords=np.array([x_v, y_v, 0]), pos=0)
        G.add_edge(atom_n, node)
        # Initial
        x = []
        y = []

print("The total node is: ", G.number_of_nodes())
# Write the XYZ file
with open("molecule.xyz", "w") as f:
    f.write(f"{len(G.nodes)}\n")
    f.write("\n")
    for node in G.nodes:
        attrs = G.nodes[node]
        element = attrs["type"]
        if element==1:
            element="C"
        elif element==2:
            element="H"
        x, y, z = attrs["coords"][:]
        f.write(f"{element} {x} {y} {z}\n")

