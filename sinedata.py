# This file is part of Python GRN implementation.
# Architype is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# Architype is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License 
# along with GRN.  If not, see <http://www.gnu.org/licenses/>.
# Author Jonathan Byrne 2014

"""
iterate grn with offset sine input
"""
import numpy, math, graph
import cgrn as grn
#import grn

def run_grn(genome, delta, syncsize, offset):
    """Create GRN, give sine as input,
    calculate fitness for each p_gene and return best
    """
    stabiter = 10000 / delta
    if syncsize > delta:
    	runiter = syncsize / delta
    else:
      runiter = syncsize
      delta = 1	
    completed = False
    results = None

    #run until delta stays above zero
    while not completed:
        regnet = grn.GRN(genome, delta)
        regnet.build_genes()
        regnet.add_extra("EXTRA_sineval", 0.2, [0]*32)
        results = [0] * len(regnet.genes)
        regnet.precalc_matrix()

        #build list of p-genes, quit if only 1 p gene
        p_genes = [ idx for idx, gn in enumerate(regnet.genes)
                   if gn.gene_type == "P"]
        if len(p_genes) < 2:
            return [0], [0]

        # stabilise the grn
        regnet.regulate_matrix(stabiter)

        # generate complete sine revolution
        sinesamples = [(math.sin(i) + 1)/2
                    for i in numpy.arange(0,4*math.pi, 0.1)]

        for inputval in sinesamples:
            extra_vals = {'sineval':inputval}
            regnet.set_extras(extra_vals)
            regnet.regulate_matrix(runiter)

        #only finish if TF concs don't go below zero
        if regnet.below_zero == False:
            completed = True
        else:
            delta = delta - 1
            stabiter = 10000 / delta
        if syncsize > delta:
            runiter = syncsize / delta
	else:
	    delta = 1

    # Iterations complete, now calculate each p-gene fitness
    for p_gene in p_genes:
        fitness = 0
        for index in range(len(sinesamples)):
            offidx = (index - offset) % len(sinesamples)
            target = sinesamples[offidx] * 0.4
            signal = regnet.conc_list[p_gene][index]
            fitness += 1 - abs(target - signal)
        results[p_gene] += fitness
    return results, regnet.conc_list

def main():
    """comparison code"""
    import random

    random.seed(1)
    genome = [random.randint(0, 1) for _ in range(0, 5000)]
    results, conc_list = run_grn(genome, delta=1, syncsize=100, offset=10)
    minresult = 126 - max(results)
    print "bestfit:", minresult
    if len(conc_list) > 1:
        graph.plot_2d(conc_list, "sinedata"+str(1))

if __name__ == "__main__":
    main()
