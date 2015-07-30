from composites.sim import Sim,ureg,Q_
from composites.failure import FailureAnalysis
import numpy as np
from composites.laminate import Laminate

import random

from operator import attrgetter
from deap import base,creator,tools,algorithms
from .algorithm import MEtwo
import makechildren
import sys
print "Yay!"
sys.exit()


# In[4]:

zc = 0.15/100.0
w = Q_(3.0/100.0,'m')
B = Q_(10.0/100.0,'m')
M_1 = Q_([-1000,-100,-100],'N')
N_1 = Q_([-22400,-3000,-2000],'Nperm')
M_2 = Q_([-980,-98,-110],'N')
N_2 = Q_([-20800,-2800,-2200],'Nperm')


# In[5]:

def evaluate_layup(individual,materialID):
    list_of_str = ['{0:d}'.format(k) for k in individual]
    layup = '/'.join(list_of_str)+'s'
    laminate=Laminate(layup,materialID,zc)
    sim_1 = Sim(laminate=laminate)
    sim_1.apply_N(N_1)
    sim_1.apply_M(M_1)
    sim_1.solve()
    fail = FailureAnalysis(sim_1)
    the_mins = fail.find_min()
    first_min = min([R[0] for R in the_mins.values()])
    sim_1.apply_N(N_2)
    sim_1.apply_M(M_2)
    fail = FailureAnalysis(sim_1)
    the_mins = fail.find_min()
    second_min = min([R[0] for R in the_mins.values()])
    return float(min((first_min,second_min))),


# In[6]:

creator.create("FitnessMax", base.Fitness, weights = (1.0,))
creator.create("Individual",list,fitness=creator.FitnessMax)


# In[7]:

def new_angle(a,b,coeff):
    return random.randint(a,b)*coeff


# In[8]:

materialID = 1
num_of_plies = 7
pop_size = 200
nk = 25
rand_counter = 15
rand_ind = 50
gen_numbers = 100


# In[9]:

toolbox = base.Toolbox()
toolbox.register("ply",new_angle,-9,9,10)
toolbox.register("individual",tools.initRepeat,creator.Individual,toolbox.ply,n=num_of_plies)
toolbox.register("population",tools.initRepeat,list, toolbox.individual,n=pop_size)
toolbox.register("evaluate",evaluate_layup,materialID=materialID)
toolbox.register("mate",makechildren.one_point_cross)


# In[ ]:

filename = "material%i--%iplies.txt" % (materialID,num_of_plies)
file_path = 'design3_results/'+filename
f = open(file_path,'a')
header = """Material ID : %i
Number of Plies : %i
Population Size : %i
Nk : %i\n\n""" % (materialID,num_of_plies,pop_size,nk)
print header,; f.write(header)
try:
    population = toolbox.population()
    fits = toolbox.map(toolbox.evaluate,population)
    for fit,ind in zip(fits,population):
        ind.fitness.values = fit
    
    sorted_pop = sorted(population, key=attrgetter("fitness"), reverse=True)
    last_fit = sorted_pop[0].fitness.values[0]
    counter = 0
    bumped_nk = False
    
    print len(population)
    print "Writing to %s" % file_path
    
    for gen in range(gen_numbers):
        population = MEtwo(population,nk,toolbox,1.0,0.01,0.9)
        sorted_pop = sorted(population, key=attrgetter("fitness"), reverse=True)
        new_fit = sorted_pop[0].fitness.values[0]
        if new_fit - last_fit < 10**-5:
            counter = counter + 1
        else:
            counter = 0
        last_fit = new_fit        
        results = "Generation %i (%i) \n%r \n %r \n" % (gen, 
                                                      counter,
                                                      sorted_pop[0],
                                                      sorted_pop[0].fitness.values
                                                      )
        print results,
        f.write(results)
        

        
        if counter == rand_counter:
            print "Reached %i" % rand_counter
            random_individuals = toolbox.population()[:rand_ind]
            sorted_pop[-rand_ind:] = random_individuals
            population = sorted_pop
            message = "Inserted %i random individuals\n" % rand_ind
            print message,; f.write(message)
            counter = 0
#                 
#             if bumped_nk:
#                 message = "Could not get higher"
#                 print message; f.write(message)
#                 raise AssertionError
#             else:
#                 random_individuals = toolbox.population()[:rand_ind]
#                 sorted_pop[-rand_ind:] = random_individuals
#                 population = sorted_pop
#                 message = "Inserted %i random individuals\n" % rand_ind
#                 counter = 0
#                 print message,; f.write(message)
#                 nk = pop_size
#                 bumped_nk = True
#                 message = "Bumping N_k to population size (%i)" % pop_size
#                 print message; f.write(message)
except AssertionError:
    pass
finally:
    f.close()
print "Done!"

