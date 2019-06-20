## SIRinTime.py
## Simulation to PLOT S,I,R fraction of population in time for ROWS/RANDOM arrangement

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
np.set_printoptions(threshold=sys.maxsize)
import time


# create grilled graph          # defined also in ComputeDistances !!!!!!!!
def CreateGraph(num_rows, num_cols):
    all_nodes = num_cols * num_rows - 1
    G = nx.Graph()
    horizontal_edges = [(i, i + 1) for i in range(num_cols * num_rows) if (i + 1 - num_cols) % num_cols != 0]
    vertical_edges = [(i, i + num_cols) for i in range(num_cols * num_rows) if i + num_cols <= all_nodes]
    all_edges = horizontal_edges + vertical_edges
    G.add_edges_from(all_edges)
    return G


# conversion from Graph to dictionary: node: neighbours, state
def Network2dict(G):
    dic = {}
    for n in G.nodes():
        values = {'neighbours': list(G.neighbors(n)), 'crop': 'None', 'state': 'None', 'p_inf': '0.2', 'p_re': '0.0'}
        dic[n] = values
    return dic


# Initialize nodes (Crop A / Crop B) and modify status field in dictionary
# Crop A initialized as state Susceptible. Crop B initialized as Removed/Recovered
# p0 will be used in Simulation function when formulas are applied for t=0
def init_nodes(dictio, num_rows, type_arrangement):
    cols = len(dictio)/num_rows
    is_infected = 0
    crops_list = []

    if type_arrangement == 'rows':
        for n in dictio:
            y = int(n / cols)
            if y % 2 == 0:
                dictio[n]['crop'] = 'A'
                if is_infected == 0:
                    dictio[n]['state'] = 'I'     # To start infection we infect 1st node of crop 'A'
                    is_infected = 1
                else:
                    dictio[n]['state'] = 'S'
            else:
                dictio[n]['crop'] = 'B'
                dictio[n]['state'] = 'R'
                dictio[n]['p_inf'] = '0.0'
            crops_list.append(dictio[n]['crop'])
    else:
        i = 0
        randoms = []
        while i < len(dictio)/2:
            r = np.random.randint(0, len(dictio))
            if r not in randoms:
                randoms.append(r)
                i += 1
            is_infected = 0
        for n in dictio:
            if n in randoms:
                dictio[n]['crop'] = 'A'
                if is_infected == 0:
                    dictio[n]['state'] = 'I'     # To start infection we infect 1st node of crop 'A'
                    is_infected = 1
                else:
                    dictio[n]['state'] = 'S'
            else:
                dictio[n]['crop'] = 'B'
                dictio[n]['state'] = 'R'
                dictio[n]['p_inf'] = '0.0'
            crops_list.append(dictio[n]['crop'])

    return dictio, crops_list


# Simulation (for all time steps)
def Simulation(dictio, max_steps, beta, esse):
#def Simulation(dictio, max_steps):
    # At t=0 half of the crop is B, then it's Recovered.
    # All the crop A is Susceptible, except for 1 node which is Infected.
    num_recovered = len(dictio) / 2
    num_susceptible = num_recovered - 1
    num_infected = 1

    # fraction of nodes S I R at the beggining of Simulation (t=0)
    sus = num_susceptible / len(dictio)
    inf = num_infected / len(dictio)
    re = num_recovered / len(dictio)

    susceptible_simulation = [sus]
    infected_simulation = [inf]
    recovered_simulation = [re]

    # Loop number of steps
    for t in range(0, max_steps):
        #print('t: {}'.format(t))
        next_susceptible = []           # next step
        next_infected = []
        next_recovered = []

        # list of nodes
        for n in dictio:
            if dictio[n]['crop'] == 'A':         # we are only interested in primary crop A (can be affected by infection)
                prod = 1
                for j in Gdict:                 # from node n to all other nodes (possible spread of disease)
                    if n == j:
                        A = 0
                    else:
                        dist = distances[n][j]
                        A = dist ** (-esse)
                    resta = beta * A * float(dictio[j]['p_inf'])
                    prod = prod * (1 - resta)

                # Formulas Markovian formulation of the epidemiological model
                qu = 1 - prod                                     # (Eq. 2.13) -- Prob Susceptible (t) to Infected (t+1)
                inf_t = float(dictio[n]['p_inf'])
                re_t = float(dictio[n]['p_re'])
                inf_next = inf_t * (1 - mu) + (1 - inf_t - re_t) * qu  # (Eq. 2.11) -- Infected (t+1)
                re_next = re_t + mu * inf_t                            # (Eq. 2.12) -- Recovered/recovered (t+1)

                r = np.random.uniform(0.0, 1.0)
                if dictio[n]['state'] == 'I':    # infected in time step t (can become I or R in t+1)
                    if r < re_next:
                        next_recovered.append(n)
                    else:
                        next_infected.append(n)
                elif dictio[n]['state'] == 'S':    # susceptible in time step t
                    if r < inf_next:
                        next_infected.append(n)
                    else:
                        next_susceptible.append(n)
                else:                               # recovered in time step t
                    next_recovered.append(n)

                dictio[n]['p_inf'] = inf_next        # update probabilities of I and R for the node
                dictio[n]['p_re'] = re_next

        for ni in next_infected:
            dictio[ni]['state'] = 'I'
        for ns in next_susceptible:
            dictio[ns]['state'] = 'S'
        for nr in next_recovered:
            dictio[nr]['state'] = 'R'

        num_susceptible = 0
        num_infected = 0
        num_recovered = 0

        for n in dictio:
            if dictio[n]['state'] == 'S':
                num_susceptible += 1
            elif dictio[n]['state'] == 'I':
                num_infected += 1
            else:
                num_recovered += 1

        # fraction of nodes S I R at the end of time step t
        sus = num_susceptible / len(dictio)
        inf = num_infected / len(dictio)
        re = num_recovered / len(dictio)

        # fraction of nodes A that got infected at the end of time step t
        #inf = float(len(next_infected)) / len(dictio)
        infected_simulation.append(inf)
        # fraction of nodes A that got recovered/recovered at the end of time step t
        #re = float(len(next_recovered)) / len(dictio)
        recovered_simulation.append(re)
        # fraction of nodes A that got susceptible at the end of time step t
        #sus = float(len(next_susceptible)) / len(dictio)
        susceptible_simulation.append(sus)

    return infected_simulation, recovered_simulation, susceptible_simulation, dictio


# PLOT population/time ROWS/RANDOM arrangement
def plot_SIR_time(inf, rem, sus, b, mu, p0, s, maxt):
    plt.plot(inf[0], 'g-', label='Infected - Rows')  # rows
    plt.plot(rem[0], 'r-', label='Recovered - Rows')
    plt.plot(sus[0], 'b-', label='Susceptible - Rows')
    plt.plot(inf[1], 'g:', label='Infected - Random')  # random
    plt.plot(rem[1], 'r:', label='Recovered - Random')
    plt.plot(sus[1], 'b:', label='Susceptible - Random')
    plt.xlabel('t')
    plt.ylabel('% population')
    plt.legend(loc=7)  # 'center right'
    plt.title('SIR(β=%.1f, μ=%.1f, p0=%.1f, s=%.1f)' % (b, mu, p0, s))
    plt.savefig('Population_beta_' + str(b) + '_mu_' + str(mu) + '_s_' + str(s) + '_t_' + str(maxt) + '.png')
    plt.close()


### MAIN ###
# SIR parameters
# infection probability of a Susceptible individual when is contacted with an Infected one
# at least 51 values (increases of 0.02)
betas = [0.1, 0.2, 0.5, 0.8, 1.0]
# mu: death rate of the disease (spontaneous recovery????????) !!!!!!
mu = 0.5
# p_i(t): prob node i Infected at time t
# p0: initial fraction of infected nodes    --  20%
p0 = 0.2
# s: aphid mobility (near 0: high, near infinit: low)
esses = [1.0, 2.5, 10.0]
# Define the size of the grids (fields)
rows = 20
columns = 50
Tmax = 20
arrangements = ['rows', 'random']         # 'rows' or 'random' arrangement of CROPS

# Matrix of short path lengths from node to node read from CSV file (to improve performance)
df = pd.read_csv("distances_grill.csv", header=None)
distances = df.values.tolist()

for s in esses:
    for b in betas:
        infected = []
        removed = []
        susceptible = []

        for i, arrangement in enumerate(arrangements):
            # Also in ComputeDistances.py
            # I need it here to create the dictionaries from it
            G_grill = CreateGraph(rows, columns)            # grid for rows and random arrangements

            start = time.time()

            # Create dictionary from Network
            Gdict = Network2dict(G_grill)

            # Initialize nodes ROWS/RANDOM arrangement
            Gdict, crops = init_nodes(Gdict, rows, arrangement)
            print('Gdict_' + arrangement + ': {}'.format(Gdict))
            print('crops_' + arrangement + ': {}'.format(crops))

            # Simulation ROWS/RANDOM arrangement
            inf_simulation, rec_simulation, sus_simulation, Gdict = Simulation(Gdict, Tmax, b, s)
            print('FINAL Gdict_' + arrangement + ': {}'.format(Gdict))
            infected.append(inf_simulation)
            removed.append(rec_simulation)
            susceptible.append(sus_simulation)

            end = time.time()
            print('    time: %.2f s' % (end - start))

            G_grill.clear()                                 # remove graph to avoid problems in next iteration

        # PLOT population/time ROWS/RANDOM arrangement
        plot_SIR_time(infected, removed, susceptible, b, mu, p0, s, Tmax)

