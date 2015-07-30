import random
from makechildren import rank_parents,createChildren
from operator import attrgetter

def MEtwo(parents,nk,toolbox,cxpb,mutpb,swappb):
    """Multiple Elitist one selection. The new generation
    will have as many parents as the input number of parents.
    
    the parent and child populations of size P are
    combined into one population, forming 2P laminates. 
    Members of the combined population are then arranged from best to worst 
    according to their fitness values. The number of top laminates
    from the combined population which will be carried over to the next 
    generation is designated N k.
    
    In ME 2 selection, laminates from the
    portion of the child population not already in the new population 
    will be randomly chosen to fill the remainder of the new population. 
    This method will allow for the possibility of some of the poorer laminates, 
    which may contain important genetic information, from the child population to
    migrate to the new population
    
    1. Make children
    2. Rank children+parents
    3. Place Nk best in new generation
    4. Randomly select remaining new_gen parents from children
        This means a children can't go in twice.
    """
    
    size = len(parents)
    if size < nk:
    	raise ValueError("Population size (%i) can't be bigger than Nk (%i)" 
    									% (size,nk))
    	
    new_gen = [None]*size
    sorted_parents,weights = rank_parents(parents)
    children = createChildren(sorted_parents,weights,size,
                              toolbox,cxpb,mutpb,swappb)
    
    fits = toolbox.map(toolbox.evaluate,children)
    for fit,ind in zip(fits,children):
        ind.fitness.values = fit
    
    sorted_children = sorted(children, key=attrgetter("fitness"), reverse=True)
    for i in xrange(nk):
        if sorted_children[0].fitness < sorted_parents[0].fitness:
            new_gen[i] = sorted_parents.pop(0)
        else:
            new_gen[i] = sorted_children.pop(0)
    
    random.shuffle(sorted_children)
    for i in xrange(nk,size):
        new_gen[i] = sorted_children.pop()
    
    assert None not in new_gen
    
    return new_gen