import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

import time

# Networks Definition
Net1 = nx.Graph(nx.read_pajek('dolphins.net'))
Net2 = nx.Graph(nx.read_pajek('ER1000k8.net'))
Net3 = nx.Graph(nx.read_pajek('SF_500_g2.7.net'))

# SIS parameters
# spontaneous recovery probability  -- (e.g.: 0.1, 0.5, 0.9 for instance)
mus = [0.1, 0.5, 0.9]
# infection probability of a Susceptible individual when is contacted with an Infected one
# at least 51 values (increases of 0.02)
betas = [0.01*i for i in range(0, 102, 2)]    # 51 values: 0 ~ 1


# conversion from Graph to dictionary: node: neighbours, state
def Network2dict(G):
    dict = {}
    for n in G.nodes():
        values = {'neighbours': list(G.neighbors(n)), 'state': 'None'}
        dict[n] = values
    return dict


# Initialize nodes (Infected (I) / Susceptible (S)) and modify status field in dictionary
def init_nodes(Gdict, p0):
    infected = []
    susceptible = []
    for n in Gdict:
        if np.random.uniform(0.0, 1.0) < p0:
            Gdict[n]['state'] = 'I'
            infected.append(n)
        else:
            Gdict[n]['state'] = 'S'
            susceptible.append(n)
    return Gdict


# Determines which neighbours of a particular node n are infected
def get_neighbours_infected(n):
    list_neighbours_infected = []
    neighbours = Gdict[n]['neighbours']
    for ne in neighbours:
        if Gdict[ne]['state'] == 'I':
            list_neighbours_infected.append(ne)
    return list_neighbours_infected


# Monte Carlo simulation -- one repetition for all time steps (until Tmax)
def one_simulation(Gdict, Tmax):
    # Loop number of steps
    phi_simulation = []
    for t in range(0, Tmax):
        next_susceptible = []           # next step
        next_infected = []
        # llista nodes
        for n in Gdict:
            if Gdict[n]['state'] == 'I':    # infected in time step t
                r = np.random.uniform(0.0, 1.0)
                if r < mu:
                    next_susceptible.append(n)
                else:
                    next_infected.append(n)
            else:                           # susceptible in time step t
                list_neighbours_infected = get_neighbours_infected(n)
                infected = 0
                i = 0
                while i < len(list_neighbours_infected):
                    r = np.random.uniform(0.0, 1.0)
                    if r < beta:
                        infected = 1
                        break
                    i += 1
                if infected == 1:
                    next_infected.append(n)
                else:
                    next_susceptible.append(n)

        for ni in next_infected:
            Gdict[ni]['state'] = 'I'
        for ns in next_susceptible:
            Gdict[ns]['state'] = 'S'

        # fraction of nodes that got infected at the end of time step t
        phi = float(len(next_infected)) / len(Gdict)
        phi_simulation.append(phi)

        #print('Gdict step: {}'.format(Gdict))
        #print('phi step: {}'.format(phi))

    return phi_simulation, Gdict


# Monte Carlo simulation -- all simulations
def MonteCarlo(Gdict, Nrep, Tmax, Ttrans, p0):
    #Loop number of simulations
    for rep in range(0, Nrep):
        phi_final = 0

        Gdict = init_nodes(Gdict, p0)
        #print('Gdict initialized: {}'.format(Gdict))

        phi_simulation, Gdict = one_simulation(Gdict, Tmax)
        #print('Gdict final: {}'.format(Gdict))
        #print('List fraction infected: {}'.format(phi_simulation))

        stationary_phi = phi_simulation[Ttrans:]
        mean_phi = sum(stationary_phi)/len(stationary_phi)
        phi_final += mean_phi

    phi_final = phi_final/Nrep

    return phi_final


### MAIN PROGRAM ###

#Nrep = 100     # number of repetitions of the simulation     -- 50 si triga molt!!!
Nrep = 50
Tmax = 1000     # max number of time steps of each simulation
Ttrans = 900    #x number of steps of the transitory
p0 = 0.2        # initial fraction of infected nodes          -- 20%

i = 0
networks = ['dolphins.net', 'ER1000k8.net', 'SF_500_g2.7.net']
for net in [Net1, Net2, Net3]:
    print('Network:', networks[i])
    for mu in mus:
        phi_list = []
        print('- mu:', mu)
        for beta in betas:
            print(' - beta:', beta)
            start = time.time()
            G = net
            # Create dictionary from Network
            Gdict = Network2dict(G)
            phi_final = MonteCarlo(Gdict, Nrep, Tmax, Ttrans, p0)
            end = time.time()
            phi_list.append(phi_final)
            print('   p: %.2f' % phi_final)
            print('   time: %.2f s' % (end-start))

        plt.plot(betas, phi_list, 'o-')
        plt.xlabel('β')
        plt.ylabel('P')
        plt.title(networks[i] + ', SIS(μ=%.1f, P0=%.1f)' % (mu, p0))
        plt.savefig(networks[i] + '_' + str(mu) + '.png')
        plt.close()

    i += 1

