import networkx as nx

# Create a molecular graphene structure
G = nx.Graph()

# Add nodes for each carbon atom in the hexagonal ring
num_rings = 3
ring_size = 6
for i in range(num_rings):
    for j in range(ring_size):
        G.add_node(i*ring_size+j, element="C")

# Add edges to form hexagonal ring 1
for i in range(ring_size):
    G.add_edge(i, (i+1)%ring_size)
    G.add_edge(i, i+ring_size)
    G.add_edge(i, (i-1+ring_size)%ring_size+2*ring_size)
for i in range(ring_size, 2*ring_size):
    G.add_edge(i, i+ring_size)

# Add edges to form hexagonal ring 2
for i in range(ring_size, 2*ring_size):
    G.add_edge(i, (i+1)%ring_size+2*ring_size)
    G.add_edge(i, i+ring_size)
    G.add_edge(i, (i-1+ring_size)%ring_size+3*ring_size)
for i in range(2*ring_size, 3*ring_size):
    G.add_edge(i, i+ring_size)

# Add edges to connect the rings
for i in range(ring_size):
    G.add_edge(i, i+2*ring_size)
    G.add_edge(i+ring_size, i+2*ring_size)
for i in range(ring_size):
    G.add_edge(i+2*ring_size, (i-1+ring_size)%ring_size+3*ring_size)
    G.add_edge(i+2*ring_size, (i+1)%ring_size+3*ring_size)

# Add nodes and edges for the phenyl ring
phenyl_nodes = [3*ring_size, 3*ring_size+1, 3*ring_size+2, 3*ring_size+3, 3*ring_size+4, 3*ring_size+5, 3*ring_size+6]
for node in phenyl_nodes:
    G.add_node(node, element="C")
G.add_edge(phenyl_nodes[0], phenyl_nodes[1])
G.add_edge(phenyl_nodes[1], phenyl_nodes[2])
G.add_edge(phenyl_nodes[2], phenyl_nodes[3])
G.add_edge(phenyl_nodes[3], phenyl_nodes[4])
G.add_edge(phenyl_nodes[4], phenyl_nodes[5])
G.add_edge(phenyl_nodes[5], phenyl_nodes[0])
G.add_edge(phenyl_nodes[0], phenyl_nodes[6])
G.add_edge(phenyl_nodes[2], phenyl_nodes[6])
G.add_edge(phenyl_nodes[4], phenyl_nodes[6])

# Write molecular graphene to .xyz file
with open("molecular_graphene.xyz", "w") as f:
    f.write(str(len(G.nodes())) + "\n")
    f.write("Molecular Graphene\n")
    for node in G.nodes():
        element = G.nodes[node]["element"]
        x, y, z = node // ring_size, node % ring_size, 0
        f.write(element + " " + str(x) + " " + str(y) + " " + str(z) + "\n")

