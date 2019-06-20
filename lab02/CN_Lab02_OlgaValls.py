# Complex Networks
# LAB02: Models of complex networks
# Olga Valls
# 20190417

import networkx as nx
import math
import matplotlib.pyplot as plt
import random
import numpy as np
from scipy.special import factorial

#
# Erdös-Rényi Network (ER)
#
def ERnet(N,p):
    # Create empty graph
    G = nx.Graph()
    # Add n nodes to it. No edges yet.
    # add_nodes_from(nodes) --> nodes: iterable container
    G.add_nodes_from(range(N))

    # Add edges to the graph randomly:
    # For all possible pairs of nodes
    # 1. Take a pair of nodes
    # 2. get a random number between 0 and 1
    # 3. if r <= p --> add edge, else ignore
    list_edges = []
    for i in range(N):
        for j in range(i+1,N):  # we don't need a->a edges or a->b and b->a
            r = random.random()
            if r <= p:
                G.add_edge(i,j)
                edge = (i,j)
                list_edges.append(edge) # keep created edges
    print('list of edges: {}'.format(list_edges))

    # Plot the Network
    if (N<=50):
        nx.draw(G)
        plt.show()
    print('nodes: {}, probability: {}'.format(N,p))

    return G

# Plot PDFs for ER Network
def plot_PDF_ER(G,p):
    # Compute experimental values for PDF
    temp = nx.degree(G)
    print('temp: {}'.format(temp))

    list_degrees = []
    for l in temp:
        list_degrees.append(l[1])

    print('k: {}'.format(list_degrees))
    list_degrees.sort()
    print('k sorted: {}'.format(list_degrees))

    k_min = list_degrees[0]
    k_max = list_degrees[-1]
    print('kmin: {}, kmax: {}'.format(k_min, k_max))

    degrees = []
    count_p = []
    for i in range(k_min, k_max + 1):
        degrees.append(i)
        count_p.append(list_degrees.count(i) / len(list_degrees))
    print('degrees: {}'.format(degrees))
    print('count_p: {}'.format(count_p))

    # Compute theoretical values for PDF
    k_av = p*N              # average distribution
    h = nx.degree_histogram(G)
    poisson = []
    #np_degrees = np.array(hists)
    #for i in range(0,len(np_degrees)):
    for i in range(len(h)):
    #for i in np_degrees:
        poisson.append(np.exp(-k_av) * k_av**i / factorial(i))
    print('hists: {}'.format(h))
    print('poisson: {}'.format(poisson))

    # Plot Linear Degree distributions (experimental (blue) and theoretical (red)
    plt.plot(degrees, count_p, 'o-', color='b', label='Experimental PDF')
    plt.plot(poisson, 'o-', color='r', label='Theoretical PDF')
    plt.title('ER Degree distribution (N=%d' % N + ' & p=%f' % p + ')')
    plt.xlabel('Degree k')
    plt.ylabel('Prob(k)')
    plt.legend(loc=2)
    plt.grid()
    plt.show()

# # ER Network parameters
# N = 50            # number of nodes
# p = 0.3           # probability [0, 1]

# # Create the ER network
# G = ERnet(N,p)

# # Plot the linear PDFs (experimental and theoretical)
# plot_PDF_ER(G,p)


#
# Barabási & Albert Network (BA)
#
# N: total nodes of the network
# m0: initial nodes
# m: edges to be connected from the new nodes m<=m0
def BAnet(N,m0,m):
    # Create a path graph with the initial number of nodes
    G = nx.path_graph(m0)
    nx.draw(G)
    plt.show()
    temp = nx.degree(G)
    print('degrees initial nodes: {}'.format(temp))
    k = [i[1] for i in temp]
    print('degrees: {}'.format(k))

    for i in range(m0,N):   # add new nodes (N-m0) to the network
        # 1. add the node to the network
        G.add_node(i)
        # 2. add m edges to this node (one by one) -- following Preferential Attachment
        # 2.1. Preprocessing
        # dictionary of degrees
        temp = nx.degree(G)
        deg = [i[1] for i in temp]
        #print('deg: {}'.format(deg))

        # dictionary of probabilities (more degree more probability)
        # p = k / sum(k)
        node_probs = {}
        for n in G.nodes():
            #print(n)
            node_probs[n] = deg[n]/sum(deg)

        # order of nodes -- list cumulative probabilities
        node_probs_cum = []
        acc = 0
        for k,v in node_probs.items():
            temp = [k,acc+v]
            node_probs_cum.append(temp)
            acc = acc + v

        # 3. while the number of edges added is != m
        nodes_used = []
        num_added = 0
        edges_used = []
        new_edges = []

        while (num_added < m):
            # choose a random number btw 0 and 1
            r = random.random()
            # ens quedem amb el primer node que té cumulative >= random number
            # more window in the cumulative "bin", more probability --> assign node to more prob.
            cum = 0
            k = 0
            while(not(node_probs_cum[k][1] >= r)):
                cum = node_probs_cum[k][1]
                k = k+1
            node = node_probs_cum[k][0]
            # if the node is not used yet and the edge is not in the list of edges, we add it to the list of nodes used
            if (node_probs_cum not in nodes_used) and (node_probs not in edges_used):
                nodes_used.append(node)
                G.add_edge(i,node)                  # create new edge from old node to new node
                num_added = num_added +1
                edges_used.append((i,node))
                new_edges.append((i,node))
        print('new_edges (it #{}): {}'.format(i,new_edges))

    print('node_probs: {}'.format(node_probs))
    print('node_probs_cum: {}'.format(node_probs_cum))

    # Plot the Network
    if N<=50:
        nx.draw(G)
        plt.show()

    return G

# Plot PDFs for BA Network
def plot_PDF_BA(G,m):
    # Compute experimental values for PDF
    temp = nx.degree(G)
    print('temp: {}'.format(temp))

    list_degrees = []
    for l in temp:
        list_degrees.append(l[1])

    print('k: {}'.format(list_degrees))
    list_degrees.sort()
    print('k sorted: {}'.format(list_degrees))

    #k_min = list_degrees[0]
    # els nous node han de connectar amb m nodes de la network inicial
    # --> ho fixem com a minim degree, tot i que la network inicial és un path i alguns tenen k=1 o k=2
    k_min = m
    k_max = list_degrees[-1]
    print('kmin: {}, kmax: {}'.format(k_min, k_max))

    degrees = []
    count_p = []
    for i in range(k_min, k_max + 1):
        degrees.append(i)
        count_p.append(list_degrees.count(i) / len(list_degrees))
    print('degrees: {}'.format(degrees))
    print('count_p: {}'.format(count_p))

    # Compute theoretical values for PDF
    #p(k) = (2m*(m+1)) / (k*(k+1)*(k+2))
    h = nx.degree_histogram(G)
    theo = []
    for i in range(m,len(h)):
        theo.append((2*m*(m+1))/(i*(i+1)*(i+2)))
    print('hists: {}'.format(h))
    print('theorethical: {}'.format(theo))

    # Plot Linear Degree Distributions (experimental (blue) and theoretical (red)
    plt.plot(degrees, count_p, 'o-', color='b', label='Experimental PDF')
    plt.plot(theo, 'o-', color='r', label='Theoretical PDF')
    plt.title('BA Degree distribution (N=%d' % N + ' & m0=%d' % m0 + ' & m=%d' % m + ')')
    plt.xlabel('Degree k')
    plt.ylabel('Prob(k)')
    plt.legend(loc=1)
    plt.grid()
    plt.show()


# Plot log-log PDFs for BA Network
def plot_PDF_BA_log(G,num_bins):
    temp = nx.degree(G)
    list_degrees = []
    for l in temp:
        list_degrees.append(l[1])
    list_degrees.sort()
    print('k sorted: {}'.format(list_degrees))

    #k_min = list_degrees[0]
    # els nous node han de connectar amb m nodes de la network inicial
    # --> ho fixem com a minim degree, tot i que la network inicial és un path i alguns tenen k=1 o k=2
    k_min = m
    k_max = list_degrees[-1]

    # 2. Compute the logarithm of ki for all the data elements
    k_log = []
    for k in list_degrees:
        k_log.append(math.log10(k))
    print('log(k): {}'.format(k_log))

    # 3. Divide the interval [log(kmin), log(kmax +1)] in equal size bins, e.g. 10 bins, to get de
    # values x0 = log(kmin), x1, x2, ..., x10 = log(kmax +1)
    x_bins_log = []
    log_k_min = k_log[0]
    x_bins_log.append(log_k_min)
    log_k_max_1 = math.log10(k_max + 1)

    size_bin = (log_k_max_1 - log_k_min) / num_bins
    x_temp = log_k_min
    for i in range(1, num_bins):
        x_temp = x_temp + size_bin
        x_bins_log.append(x_temp)
    x_bins_log.append(log_k_max_1)
    print('bins: {}; log(kmin): {}; log(kmax +1): {}'.format(len(x_bins_log) - 1, log_k_min, log_k_max_1))
    print('x_bins_log: {}'.format(x_bins_log))

    # 4. Count how many elements ki have their log(ki) in each bin [x0, x1), [x1, x2), [x2, x3),
    # ..., [x9, x10)
    histo_log = [0] * num_bins  # inicialitzo comptador de cada bin a 0
    for i in k_log:
        for b, x in enumerate(x_bins_log):
            if (i >= x):
                baix = b
            else:
                pass
        histo_log[baix] = histo_log[baix] + 1
    print('histograma log: {}'.format(histo_log))

    x_bins_log.pop(0)
    print('x_bins_log (without first): {}'.format(x_bins_log))

    # 5. Dividing the number of elements in each bin by the total number of elements n we
    # get estimations for the probabilities pb of bin [xb-1, xb)
    pb_bin_log = [h / len(list_degrees) for h in histo_log]
    print('estimations probabilities log: {}'.format(pb_bin_log))

    '''
    # Compute CCDF
    ccdf_log = []
    acc = 0
    for i, p in enumerate(reversed(pb_bin_log)):
        posi = len(pb_bin_log) - i - 1
        ccdf_log.append(pb_bin_log[posi] + acc)
        acc = acc + pb_bin_log[posi]

    print('ccdf_log normal: {}'.format(ccdf_log))
    ccdf_log = list(reversed(ccdf_log))
    print('ccdf_log reversed: {}'.format(ccdf_log))
    '''

    # Arrodoneixo a 4 decimals per a que no ocupin tant els números a la gràfica
    x_bins_log = [float(round(x, 4)) for x in x_bins_log]
    pb_bin_log = [float(round(p, 4)) for p in pb_bin_log]
    #ccdf_log = [float(round(c, 4)) for c in ccdf_log]


    # Plot log-log Degree distributions (experimental (blue) and theoretical (red)
    plt.plot(x_bins_log, pb_bin_log, 'o-', color='b')
    #plt.plot(x_bins_log, ccdf_log, 'o-', color='r', label='Theoretical PDF')
    plt.title('BA log-log Degree distribution (N=%d' % N + ' & m0=%d' % m0 + ' & m=%d' % m + ')')
    plt.xlabel('log(k)')
    plt.ylabel('Estimation p(k)')
    plt.legend(loc=1)
    plt.grid()
    plt.show()


# # BA Network parameters
# N=10000            # total number of nodes
# m0=750            # initial number of nodes
# m=120             # edges to be connected from the new node m<=m0

# # Create the BA network
# G = BAnet(N,m0,m)

# # Plot the linear PDFs (experimental and theoretical)
# plot_PDF_BA(G,m)

# # Plot the log-log PDFs (experimental and theoretical)
# b = 20;                 # number of bins
# plot_PDF_BA_log(G,b)

