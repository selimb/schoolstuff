"""Provides necessary functions to make children according to [Soremekun2000].
The application is "composite laminates". 

To be used with DEAP.

If you plan on using this module, you should register the following:
    1. either one_point_cross or two_point_cross in toolbox.mate
    2. a ply generator in toolbox.ply

You can then call "createChildren(parents,toolbox,cxpb,mutpb,swappb)"
whenever you need to create a population of children.

This does not modify the input parents
"""

from operator import attrgetter
from itertools import izip
import random
# from composites.sim import Sim,ureg,Q_
# from composites.failure import FailureAnalysis
import numpy as np

def rank_parents(parents):
    """This ranks the parents according to a 
    LINEAR roulette wheel.
    
    Returns a list of weights and 
    a list of sorted parents.
    
    As an example, the weights will be 
    [0.5,0.33,0.17] for a list of 3 parents.
    """
    sorted_parents = sorted(parents, 
                            key=attrgetter("fitness"), 
                            reverse=True)
    P = len(parents)
    i_range = xrange(1, P + 1)
    def get_rank(i):
        return 2.0*(P-i+1)/(P*(P+1))
    weights = [get_rank(k) for k in i_range]
    
    return sorted_parents,weights


#http://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python#id2
def roulette(weights):
    """Turns roulette wheel and returns the INDEX
    of the chosen individual. To be used in conjunction
    with rank_parents()
    """
    rnd = random.random() * sum(weights)
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return i


# In[229]:

def select_parents(weights):
    """Selects a pair of parents to be crossed to
    form a child via a linear roulette.
    
    The parents must be unique!
    
    Returns the indices corresponding to the parents.
    """
    while True:
        index_1 = roulette(weights)
        index_2 = roulette(weights)
        if index_1 != index_2:
            break
    
    return index_1,index_2

def two_point_cross(ind1,ind2):
    assert(len(ind1)==len(ind2))
    size = len(ind1)
    cxpoint1 = random.randint(1, size)
    cxpoint2 = random.randint(1, size - 1)
    print cxpoint1
    print cxpoint2
    if cxpoint2 >= cxpoint1:
        cxpoint2 += 1
    else: # Swap the two cx points
        cxpoint1, cxpoint2 = cxpoint2, cxpoint1
   
    ind1[cxpoint1:cxpoint2] = ind2[cxpoint1:cxpoint2]
    return ind1

def one_point_cross(ind1,ind2):
    assert(len(ind1)==len(ind2))
    size = len(ind1)
    cxpoint = random.randint(1, size-1)

    ind1[cxpoint:] = ind2[cxpoint:]
    return ind1

def breed(sorted_parents,weights,k,toolbox,cxpb):
    """Makes *k* children from the list of parents
    with crossover-rate cxpb.
    
    The parents are first ranked w.r.t fitness. 
    If crossover occurs:
        two unique parents are selected using a linear roulette 
        wheel to form a single child.
    otherwise:
        a parent is cloned in the child string
    """
    # children = [toolbox.clone(ind) for ind in sorted_parents]
    children = [None]*k
    # Apply crossover. 
    for i in range(k):
        indices = select_parents(weights)
        if random.random() < cxpb:
            assert len(indices) == 2
            parents = [toolbox.clone(sorted_parents[j]) for j in indices]
            assert len(parents) == 2
            children[i] = toolbox.mate(*parents)
            del children[i].fitness.values
        else:
            children[i] = toolbox.clone(sorted_parents[indices[0]])
            # We don't need to delete its fitness value!
    
    assert(None not in children)
    return children

def mutAddRandom(ind,new_angle,mutpb):
    """ In this implementation, each gene is given a
    small probability to switch to any other permissible in-
    teger value except the value of the gene before ply al-
    teration occurs.
    
    new_angle : toolbox function used to generate ply angle.
    I could've required "toolbox" instead but I want to provide
    the user freedom over what he calls it.
    """
    for i,gene in enumerate(ind):
        r = random.random()
        if r < mutpb:
            while True:
                new_gene = new_angle()
                if new_gene != gene:
                    break
            ind[i]=new_gene

def mutPlySwap(ind,swappb):
    """ The gene swap operator is implemented by randomly 
    selecting two unique genes in the string and switching 
    their positions.
    """
    r = random.random()
    if r > swappb:
        return
    
    size = len(ind)-1
    while True:
        gene1 = random.randint(0,size)
        gene2 = random.randint(0,size)
        if gene1 != gene2:
            break
    
    ind[gene1],ind[gene2] = ind[gene2],ind[gene1]


def createChildren(sorted_parents,weights,k,toolbox,cxpb,mutpb,swappb):
    """Algorithm from Soremekun[2000]. 
    
    1. Rank parents
    2. Breed children (crossover)
    3. Apply mutation (random insert and replace)
    4. Gene swap
    
    Returns a population of children.
    
    The number of children is equal to the number
    of parents.

    Requires toolbox attributes:
    o ply
    o mate (either one_point_cross or two_point_cross defined in this module)
    """
    #1.
    # weights,sorted_parents = rank_parents(parents)
    #2. 
    children = breed(sorted_parents,weights,k,toolbox,cxpb)
    #3/#4
    for child in children:
        mutAddRandom(child,toolbox.ply,mutpb)
        mutPlySwap(child,swappb)
    
    return children

# import random

# from deap import base
# from deap import creator
# from deap import tools
# from deap import algorithms


# # In[167]:

# creator.create("FitnessMax", base.Fitness, weights = (1.0,))
# creator.create("Individual",list,fitness=creator.FitnessMax)


# # In[168]:

# def new_angle(a,b,coeff):
#     return random.randint(a,b)*coeff


# # In[169]:

# def evaluate_layup(individual):
#     list_of_str = ['p{0:d}'.format(k) for k in individual]
#     layup = '/'.join(list_of_str)+'s'
#     sim_1 = Sim(layup=layup,
#                 materialID=materialID)
#     sim_1.apply_N(N)
#     sim_1.solve()
#     fail = FailureAnalysis(sim_1)
#     the_mins = fail.find_min()
#     return float(min([R[0] for R in the_mins.values()])),


# # In[290]:

# toolbox = base.Toolbox()
# toolbox.register("ply",new_angle,0,9,10)
# toolbox.register("individual",tools.initRepeat,creator.Individual,toolbox.ply,n=4)
# toolbox.register("population",tools.initRepeat,list, toolbox.individual,n=3)
# toolbox.register("evaluate",evaluate_layup)
# # toolbox.register("evaluate",evalOneMax)
# toolbox.register("mate",one_point_cross)
# toolbox.register("mutate",tools.mutFlipBit,indpb = 0.1)
# toolbox.register("select",tools.selNSGA2)


# # In[291]:

# materialID = 4; 
# P = Q_(1.25,'MPa'); D = Q_(0.08,'m')
# N1 = P*D/4; N2 = P*D/2
# N = Q_([N1.magnitude,N2.magnitude,0],N1.units)


# # In[292]:

# population = toolbox.population()
# fits = toolbox.map(toolbox.evaluate,population)
# for fit,ind in zip(fits,population):
#     ind.fitness.values = fit


# # In[293]:

# population


# # In[339]:

# createChildren(population,toolbox,1.0,0.01,0.9)

