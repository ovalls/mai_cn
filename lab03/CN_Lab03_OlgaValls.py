import networkx as nx
from networkx.algorithms import community
from networkx.algorithms.community import greedy_modularity_communities
from networkx.algorithms.community import k_clique_communities
from networkx.algorithms.community import asyn_fluidc
import community

import matplotlib.pyplot as plt
import igraph

# Networks to analyze:
f01 = 'A3-networks/toy/20x2+5x2'
f02 = 'A3-networks/toy/graph3+1+3'
f03 = 'A3-networks/toy/graph4+4'
f04 = 'A3-networks/toy/star'
f05 = 'A3-networks/model/256_4_4_2_15_18_p'
f06 = 'A3-networks/model/256_4_4_4_13_18_p'
f07 = 'A3-networks/model/rb125'
f08 = 'A3-networks/real/airports_UW'
f09 = 'A3-networks/real/cat_cortex_sim'
f10 = 'A3-networks/real/dolphins'
f11 = 'A3-networks/real/football'
f12 = 'A3-networks/real/zachary_unwh'

filenames = [f01, f02, f03, f04, f05, f06, f07, f08, f09, f10, f12]

##################### NETWORKX ########################

for fname in filenames:
    print(fname.upper())
    G = nx.read_pajek(fname+'.net')
    G = nx.Graph(G)

    num_nodes_G = nx.number_of_nodes(G)  # number of nodes of graph
    nodes_G = G.nodes  # number of nodes of graph
    print('Number of Nodes: {}'.format(num_nodes_G))
    print('Nodes: {}'.format(nodes_G))
    print('Number of Edges: {}'.format(nx.number_of_edges(G)))

    # Plot original graph
    #nx.draw(G)
    #plt.show()

    ### MODULARITY Community: Clauset-Newman-Moore greedy modularity maximization
    #############################################################################
    c = list(greedy_modularity_communities(G))
    # number of clusters
    print('Number of clusters (Clauset-Newman-Moore): {}'.format(len(c)))

    colors_dic = {}
    csset = []
    for i in range(len(c)):
        c_sorted = sorted(c[i])
        csset.append(set(c_sorted))
        print('nodes in cluster (Clauset-Newman-Moore) {}: {}'.format(i, c_sorted))
        for n in c_sorted:
            #print('n {}, i {}'.format(n,i))
            colors_dic[n]=i                 # cluster where each node belongs to

    print('cluster dictionary (Clauset-Newman-Moore): {}'.format(colors_dic))
    colors = []
    for key in sorted(colors_dic.keys()):
        #print("%s: %s" % (key, mydict[key]))
        colors.append(colors_dic[key])
    print('cluster/node sorted (Clauset-Newman-Moore): {}'.format(colors))
    print('list sets Modularity (Clauset-Newman-Moore): {}'.format(csset))

    # MODULARITY value
    #mod_value = nx.algorithms.community.modularity(G, [{0, 1, 2}, {3, 4, 5}])
    mod_value = nx.algorithms.community.modularity(G, csset)
    print('Modularity (Clauset-Newman-Moore): {}'.format(mod_value))


    # Save partitions in .clu
    file = open(fname+'_clauset_part.clu', 'w')
    file.write('*Vertices {}\n'.format(len(G.nodes())))
    for i in colors:
        file.write(str(i) + '\n')
    file.close()

    # Plot the different clusters -- Clauset-Newman-Moore greedy modularity maximization
    #nx.draw(G, node_color=colors)
    #plt.show()

    # Kamada-Kawai Plot view -- Clauset-Newman-Moore greedy modularity maximization
    if (fname != f07):
        pos = nx.kamada_kawai_layout(G)
    else:  # airports amb un altre layout
        pos = nx.spring_layout(G)
    print('pos: {}'.format(pos))
    x = [pos[i][0] for i in pos]
    y = [pos[i][1] for i in pos]
    print('x: {}'.format(x))
    print('y: {}'.format(y))
    plt.scatter(x, y, c=colors)
    plt.show()


    # FLUID Communities
    # based on the idea of fluids interacting in an environment,
    # expanding and contracting as a result of that interaction.
    # Fluid Communities is based on the propagation methodology,
    # which represents the state-of-the-art in terms of computational cost and scalability.
    # While being highly efficient, Fluid Communities is able to find communities in synthetic graphs
    # with an accuracy close to the current best alternatives.
    # Additionally, Fluid Communities is the first propagation-based algorithm capable of
    # identifying a variable number of communities in network.
    c = list(asyn_fluidc(G, 4, max_iter=100))
    print('Number of clusters (Fluid): {}'.format(len(c)))

    csset = []
    colors_dic = {}
    for i in range(len(c)):
        c_sorted = sorted(c[i])
        csset.append(set(c_sorted))
        print('nodes in cluster (Fluid) {}: {}'.format(i, c_sorted))
        for n in c_sorted:
            colors_dic[n]=i                 # cluster where each node belongs to
    print('cluster dictionary (Fluid): {}'.format(colors_dic))

    colors = []
    for key in sorted(colors_dic.keys()):
        colors.append(colors_dic[key])
    print('cluster/node sorted (Fluid): {}'.format(colors))
    print('list sets (Fluid): {}'.format(csset))

    # Modularity value
    mod_value = nx.algorithms.community.modularity(G, csset)
    print('Modularity (Fluid): {}'.format(mod_value))

    # save partitions in .clu
    file = open(fname + '_fluid_part.clu', 'w')
    file.write('*Vertices {}\n'.format(len(G.nodes())))
    for i in colors:
        file.write(str(i) + '\n')
    file.close()

    # Kamada-Kawai Plot view -- Fluid
    if (fname != f07):
        pos = nx.kamada_kawai_layout(G)
    else:               # airports amb un altre layout
        pos = nx.spring_layout(G)
    print('pos: {}'.format(pos))
    x = [pos[i][0] for i in pos]
    y = [pos[i][1] for i in pos]
    print('x: {}'.format(x))
    print('y: {}'.format(y))
    plt.scatter(x, y, c=colors)
    plt.show()

##################### IGRAPH + NetworkX ########################
# Louvain Algorithm

for fname in filenames:
    print(fname.upper())

    # load data into a graph
    G = igraph.read(fname + '.net')
    c = G.community_multilevel()
    print('num_clusters: {}'.format(len(c)))

    csset = []
    colors_dic = {}
    for i in range(len(c)):
        c_sorted = sorted(c[i])
        csset.append(set(c_sorted))
        print('nodes in cluster (Louvain) {}: {}'.format(i, c_sorted))
        for n in c_sorted:
            colors_dic[n]=i                 # cluster where each node belongs to
    print('cluster dictionary (Louvain): {}'.format(colors_dic))

    colors = []
    for key in sorted(colors_dic.keys()):
        colors.append(colors_dic[key])
    print('cluster/node sorted (Louvain): {}'.format(colors))
    print('list sets (Louvain): {}'.format(csset))

    # Modularity value
    mod_value = G.modularity(c)
    print('Modularity (Louvain): {}'.format(mod_value))

    # save partitions in .clu
    G_nx = nx.read_pajek(fname+'.net')
    G_nx = nx.Graph(G_nx)
    file = open(fname + '_louvain_part.clu', 'w')
    file.write('*Vertices {}\n'.format(len(G_nx.nodes())))
    for i in colors:
        file.write(str(i) + '\n')
    file.close()

    # Kamada-Kawai Plot view -- Louvain -- we plot it with NetworkX
    if (fname != f07):
        pos = nx.kamada_kawai_layout(G_nx)
    else:               # airports amb un altre layout
        pos = nx.spring_layout(G_nx)
    print('pos: {}'.format(pos))
    x = [pos[i][0] for i in pos]
    y = [pos[i][1] for i in pos]
    print('x: {}'.format(x))
    print('y: {}'.format(y))
    plt.scatter(x, y, c=colors)
    plt.show()
