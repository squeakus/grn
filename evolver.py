""" 
Evolver uses a population based ES to evolve GRNS. 
It can run on a single core or distribute the evaluations.
"""

import popstrat, pp, sys, random, graph
from data import run_grn

def singlecore(filename, delta, popsize, generations):
    """Run serially"""
    #set up the evo strategy
    best_list, mut_list = [], []
    evo = popstrat.Evostrategy(5000, popsize)
    children = evo.iterate(evo.pop)

    for i in range(generations):
        for child in children:
            results, conclist = run_grn(child['genome'], delta)
            bestidx = results.index(max(results))
            child['fitness'] = results[bestidx]
            print "fitness:", child['fitness']
            
        children = evo.iterate(children)
        bestgenome = evo.pop[-1]['genome']
        results, conclist = run_grn(bestgenome, delta)
        filename = "best_gen_"+str(i)
        graph.plot_2d(conclist, filename)
    
        if evo.adaptive:
            evo.adapt_mutation()

        best_list.append(evo.pop[-1]['fitness'])
        mut_list.append(evo.mut_rate)
        
    print "best overall fitness", evo.pop[-1]['fitness']
    
    graph.plot_2d([best_list], 'bestfit')
    graph.plot_2d([mut_list], 'mutrate')

def multicore(filename, syncsize, offset, delta, popsize, generations):
    """ Uses parallel python to evaluate, PP can also be used to
    distribute evaluations to machines accross the network"""
    #set up the evo strategy
    
    evo = popstrat.Evostrategy(5000, popsize)
    children = evo.iterate(evo.pop)

    nodes = ("*",)
    job_server = pp.Server(ncpus=8, ppservers=nodes)
    print "Starting pp with", job_server.get_ncpus(), "workers"

    for _ in range(generations):
        jobs = [(child, job_server.submit(run_grn, 
                                          (child['genome'], delta,
                                           syncsize, offset),
                                           (),
                                           ("cgrn as grn","numpy","math")))
                                           for child in children]
        for child, result in jobs:
            results, conclist = result()
            bestidx = results.index(max(results))
            child['fitness'] = results[bestidx]

        print "gen:", evo.gen_count, "fitness:", evo.pop[-1]['fitness']
        children = evo.iterate(children)

        if evo.adaptive:
            evo.adapt_mutation()

        res_file = open(filename,"a")
        res_file.write(str(evo.pop[-1])+'\n')
        res_file.close()

    restopname = filename+".rest"
    #plotting the best with colors
    bestgenome = evo.pop[-1]['genome']
    bestresult, conclist = run_grn(bestgenome, delta, syncsize, offset, restopname)

    bestidx = bestresult.index(max(bestresult))
    fitness = round(max(bestresult),1)
    colors = []

    for idx, result in enumerate(bestresult):
        # draw the input in black
        if idx == len(bestresult)-1:
            colors.append('k')
        # draw the best in green
        elif idx == bestidx:
            colors.append('g')
        # draw p_genes in red
        elif result > 0:
            colors.append('r')
        # draw TFs in blue
        else: 
            colors.append('b')

    print "colors", colors
    graphname = filename.split('/')[2]
    graphname = graphname[:-4] + "-F" + str(fitness)
    graph.plot_2d(conclist, graphname, colors, (0,1))        
        
def main():
    import os
    """Run an experiment"""
    if(len(sys.argv) != 4):
        print "Usage: " + sys.argv[0] + " <rseed> <syncsize> <offset>"
        sys.exit(1)
    seed = int(sys.argv[1])
    random.seed(seed)
    syncsize = int(sys.argv[2])
    offset = int(sys.argv[3])

    #save in named folder
    name = "sync"+str(syncsize)+"offset"+str(offset)
    pathname = "results/"+name+"/"
    filename = pathname + name + "-" + str(seed) + ".dat"
    if not os.path.exists(pathname): os.makedirs(pathname)
    multicore(filename, syncsize, offset,
              delta=1, popsize=10, generations=50)
        
if __name__ == "__main__":
    main()
