# Complex Networks
# LAB01: Structural descriptors of complex networks
# Olga Valls
# 20190318

import networkx as nx
import math
import matplotlib.pyplot as plt

# Networks to analyze:
f01 = 'A1-networks/toy/circle9.net'
f02 = 'A1-networks/toy/star.net'
f03 = 'A1-networks/toy/graph3+1+3.net'
f04 = 'A1-networks/toy/grid-p-6x6.net'
f05 = 'A1-networks/model/homorand_N1000_K4_0.net'
f06 = 'A1-networks/model/ER1000k8.net'
f07 = 'A1-networks/model/SF_1000_g2.7.net'
f08 = 'A1-networks/model/ws1000.net'
f09 = 'A1-networks/real/zachary_unwh.net'
f10 = 'A1-networks/real/airports_UW.net'

filenames = [f01, f02, f03, f04, f05, f06, f07, f08, f09, f10]

# Networks to plot:
filenames_plot = [f06, f07, f08, f10]

# Number of bins
num_bins = 20;

### 1. Numerical Descriptors ###
for fname in filenames:
    print(fname.upper())
    G = nx.read_pajek(fname)
    G = nx.Graph(G)

    print('Number of Nodes: {}'.format(nx.number_of_nodes(G)))
    print('Number of Edges: {}'.format(nx.number_of_edges(G)))

    temp = nx.degree(G)
    list_degrees = []
    for l in temp:
        list_degrees.append(l[1])

    #print('k: {}'.format(list_degrees))
    list_degrees.sort()
    print('k sorted: {}'.format(list_degrees))

    k_min = list_degrees[0]
    k_max = list_degrees[-1]

    print('Minimum Degree: {}'.format(k_min))
    print('Maximum Degree: {}'.format(k_max))
    print('Average Degree: {}'.format(sum(list_degrees)/len(list_degrees)))

    print('Average Clustering Coefficient: {}'.format(nx.average_clustering(G)))
    print('Assortativity: {}'.format(nx.degree_assortativity_coefficient(G)))
    print('Average Path Length: {}'.format(nx.average_shortest_path_length(G)))
    print('Diameter of Network: {}'.format(nx.diameter(G)))
    print('\n')

    ### 2. Plot PDF and CCDF (linear and log scales) ###
    if (fname in filenames_plot):
        name = fname.split('/')[-1]
        print('name of network: {}'.format(name))
        # 1. Find kmin = min(k) and kmax = max(k)
        #print('k sorted: {}'.format(list_degrees))
        print('k_min: {}'.format(k_min))
        print('k_max: {}'.format(k_max))

        # Compute Histogram
        bins_k = k_max + 1 - k_min
        degrees = []
        count_p = []
        for i in range(k_min,k_max+1):
            degrees.append(i)
            count_p.append(list_degrees.count(i)/len(list_degrees))
        print('degrees: {}'.format(degrees))
        print('count_p: {}'.format(count_p))

        # Compute CCDF
        ccdf = []
        acc = 0

        for i, p in enumerate(reversed(count_p)):
            posi = len(count_p) - i - 1
            ccdf.append(count_p[posi] + acc)
            acc = acc + count_p[posi]

        print('ccdf normal: {}'.format(ccdf))
        ccdf = list(reversed(ccdf))
        print('ccdf reversed: {}'.format(ccdf))

        # Linear-scale Plot PDF and CCDF
        fig1 = plt.figure(figsize=(12,5)) #inches
        ax1 = fig1.add_subplot(1, 2, 1)
        ax2 = fig1.add_subplot(1, 2, 2)
        ax1.set_title('Degree PDF (Linear Scale) -- ' + str(name), fontsize='medium')
        ax1.bar(degrees, count_p, 1, align='center', color='Blue', edgecolor='black')
        ax1.set_xlabel('Degree k')
        ax1.set_ylabel('Prob(k)')
        ax2.set_title('Degree CCDF (Linear Scale) -- ' + str(name), fontsize='medium')
        ax2.bar(degrees, ccdf, 1, align='center', color='Green', edgecolor='black')
        ax2.set_xlabel('Degree k')
        ax2.set_ylabel('Acc Prob(k)')
        #plt.show()
        fig1.savefig(str(name) + '_linear_plot.pdf', bbox_inches='tight')
        plt.close(fig1)

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
        log_k_max_1 = math.log10(k_max+1)

        size_bin = (log_k_max_1 - log_k_min) / num_bins
        x_temp = log_k_min
        for i in range(1,num_bins):
            x_temp = x_temp + size_bin
            x_bins_log.append(x_temp)
        x_bins_log.append(log_k_max_1)
        print('bins: {}; log(kmin): {}; log(kmax +1): {}'.format(len(x_bins_log) - 1, log_k_min, log_k_max_1))
        print('x_bins_log: {}'.format(x_bins_log))

        # 4. Count how many elements ki have their log(ki) in each bin [x0, x1), [x1, x2), [x2, x3),
        # ..., [x9, x10)
        histo_log = [0] * num_bins    # inicialitzo comptador de cada bin a 0
        for i in k_log:
            for b,x in enumerate(x_bins_log):
                if (i>=x):
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

        # Arrodoneixo a 4 decimals per a que no ocupin tant els números a la gràfica
        x_bins_log = [float(round(x, 4)) for x in x_bins_log]
        pb_bin_log = [float(round(p, 4)) for p in pb_bin_log]
        ccdf_log = [float(round(c, 4)) for c in ccdf_log]

        # Log-scale Plot PDF and CCDF
        fig2 = plt.figure(figsize=(12,5)) #inches
        ax1 = fig2.add_subplot(1, 2, 1)
        ax2 = fig2.add_subplot(1, 2, 2)
        ax1.set_title('Degree PDF (Log Scale) -- ' + str(name), fontsize='medium')
        ax2.set_title('Degree CCDF (Log Scale) -- ' + str(name), fontsize='medium')

        width = x_bins_log[1] - x_bins_log[0]
        ax1.bar(x_bins_log, pb_bin_log, align='edge', width=width, color='Blue', edgecolor='black')
        ax2.bar(x_bins_log, ccdf_log, align='edge', width=width, color='Green', edgecolor='black')

        ax1.set_xticks(x_bins_log)
        ax1.set_xticklabels(x_bins_log, rotation=70, ha='right')
        ax1.set_xlabel('log(k)')
        for tick in ax1.xaxis.get_ticklabels():
            tick.set_fontsize('small')
        ax1.set_yticks(pb_bin_log)
        ax1.set_yticklabels(pb_bin_log, rotation=40)
        ax1.set_ylabel('Estimation p(k)')
        for tick in ax1.yaxis.get_ticklabels():
            tick.set_fontsize('small')
        ax2.set_xticks(x_bins_log)
        ax2.set_xticklabels(x_bins_log, rotation=70, ha='right')
        ax2.set_xlabel('log(k)')
        for tick in ax2.xaxis.get_ticklabels():
            tick.set_fontsize('small')
        ax2.set_yticks(ccdf_log)
        ax2.set_yticklabels(ccdf_log, rotation=40)
        ax2.set_ylabel('Acc Estimation p(k)')
        for tick in ax2.yaxis.get_ticklabels():
            tick.set_fontsize('small')

        #plt.show()
        fig2.savefig(str(name) + str(num_bins) +'bins_log_plot.pdf', bbox_inches='tight')
        plt.close(fig2)
