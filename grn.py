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
NEW Gene Regulatory Model based on the work of Wolfgang Banzhaf
and Miguel Nicolau.
"""
#fix asynchronous update
#zero it at the very end

import random, math, graph, sys
from banzhaf_parser import reorder_gene
import numpy as np    # for matrices

class Gene:
    """
    Genes objects for GRN model. Generates protein and stores
    concentration
    """
    def __init__(self, info, gene_type):
        self.promoter =  info[0:32]
        self.enhancer = info[32:64]
        self.inhibitor = info[64:96]
        self.geneinfo = info[96:256]
        self.protein = []
        self.create_protein()
        self.concentration = 0.0
        self.gene_type = gene_type

    def create_protein(self):
        """Use geneinfo and majority rule to make protein"""
        blocks = len(self.geneinfo) / 32

        for bit in range(32):
            values = []

            for block in range(0, blocks):
                values.append(self.geneinfo[bit + block * 32])

            values.sort() #use midpoint of sort to find majority
            self.protein.append(values[len(values)/2])

class GRN:
    """ Simplified GRN model"""
    def __init__(self, genome=None, delta=1):
        """Default constructor, randomly initialised genome"""
        if genome == None:
            self.genome = self.random_init(5000)
        else:
            self.genome = genome

        self.sensible_parsing = False
        self.promoters = [[[0]*8, "TF"], [[1]*8, "P"]]
        self.delta = delta
        self.below_zero = False
        self.rest_delta = 0.00005 * delta
        self.extra_limit = 0.4
        self.tf_genes = 0
        self.p_genes = 0
        self.extras = 0
        self.genes = []
        self.conc_list = []
        self.weight_matrix = []

    def random_init(self, genome_size):
        """randomly populates a genome with ones and zeros"""
        indiv = [random.randint(0, 1) for _ in range(0, genome_size)]
        return indiv

    def read_genome(self, filename):
        """parse genome list from file"""
        gene_file = open(filename, "r")
        self.genome = eval(gene_file.readline())

    def build_genes(self):
        """Search genome for promoter sites and create gene list"""
        index = 0
        max_length = len(self.genome) - 256

        if not self.sensible_parsing:
            index = 96 - len(self.promoters[0][0])
            max_length = len(self.genome) - (160 + len(self.promoters[0][0]))

        while index <= max_length:
            for seq in self.promoters:
                found = self.genome[index:index+len(seq[0])] == seq[0]

                if found and index <= max_length:
                    if seq[1] == "TF": self.tf_genes += 1
                    elif seq[1] == "P": self.p_genes += 1

                    if self.sensible_parsing:
                        gene_segment = self.genome[index:index + 256]
                    else:
                        gene_segment = reorder_gene(index,
                                                    seq[0],
                                                    self.genome)
                    self.genes.append(Gene( gene_segment,
                                           seq[1]))
                    index += 255
            index += 1
        self.update_concentrations(True)


    def add_extra(self, gene_type, concentration, signature=None):
        self.extras +=1
        """Initialise and add a gene"""
        segment = [random.randint(0, 1) for _ in range(256)]
        extra_gene = Gene(segment, gene_type)
        extra_gene.concentration = concentration

        if not signature == None:
            extra_gene.protein = signature

        self.genes.append(extra_gene)
        self.conc_list.append([])

    def set_extras(self, extra_vals):
        """Set all EXTRA protein concentrations"""
        e_total = sum([ extra_vals[key] for key in extra_vals])
        for key in extra_vals:
            for gene in self.genes:
                extra_name = 'EXTRA_'+key
                if gene.gene_type == extra_name:
                    if len(extra_vals) == 1:
                        gene.concentration = extra_vals[key]*self.extra_limit
                    else:
                        gene.concentration = ((extra_vals[key] / e_total)
                                              * self.extra_limit)

    def precalc_matrix(self):
        """Generate concentration matrix using interdependent TF and
        input EXTRA proteins, XOR the protein with the enhancer and
        inhibitor sites on every the gene and calculate exponential.
        """
        enh_matrix = np.zeros((len(self.genes),
                               len(self.genes)-self.extras))
        inh_matrix = np.zeros((len(self.genes),
                               len(self.genes)-self.extras))

        #iterate proteins and genes to create matrix
        for idx, gene in enumerate(self.genes):
            if not gene.gene_type == "P":
                protein = gene.protein

                for idy in range(len(self.genes)-self.extras):
                    target = self.genes[idy]
                    enhancer, inhibitor = target.enhancer, target.inhibitor

                    xor_enhance = sum(protein[i] != enhancer[i]
                                      for i in range(len(protein)))
                    xor_inhibit = sum(protein[i] != inhibitor[i]
                                      for i in range(len(protein)))

                    enh_matrix[idx][idy] = xor_enhance
                    inh_matrix[idx][idy] = xor_inhibit


        #ensure it is always a negative number
        max_observed =  max(np.max(enh_matrix), np.max(inh_matrix))
        enh_matrix = enh_matrix - max_observed
        inh_matrix = inh_matrix - max_observed

        # precalculate the exponential function and generate weight matrix
        vector_exp = np.vectorize(math.exp)
        enh_matrix = vector_exp(enh_matrix)
        inh_matrix = vector_exp(inh_matrix)
        self.weight_matrix = enh_matrix - inh_matrix

    def regulate_matrix(self, iterations):
        """ Grind the concs against the weights and normalise"""
        for itr in range(iterations):
            updates = []

            for idy in range(len(self.genes)-self.extras):
                enhance, inhibit, signal = 0, 0, 0

                for idx in range(0, len(self.weight_matrix)):
                    weight = self.weight_matrix[idx][idy]
                    signal += self.genes[idx].concentration * weight

                # scale by number of proteins
                signal = signal / (self.tf_genes + self.extras)
                updates.append(signal)

            #THIS FOR LOOP UPDATES SYNCHRONOUSLY
            #for idy in range(len(self.genes)-self.extras):
                signal = updates[idy]

                if self.genes[idy].gene_type == "P":
                    self.genes[idy].concentration += self.delta * signal
                elif self.genes[idy].gene_type == "TF":
                    signal = signal * self.genes[idy].concentration
                    self.genes[idy].concentration += self.delta * signal
                    # Check if TFs(only!) fall below zero
                    if self.genes[idy].concentration < 0:
                        self.below_zero = True #delta too damn high
                # never goes below zero (computed before norm)
                if self.genes[idy].concentration < 1e-10:
                    self.genes[idy].concentration = 1e-10

            # normalise and add to the conc array
            self.update_concentrations(False)

        for idx, gene in enumerate(self.genes):
            self.conc_list[idx].append(gene.concentration)

    def update_concentrations(self, initialising):
        """Normalise the values of the TF and P proteins and add them
        to the concentration list"""
        if initialising:
            self.conc_list = []
            for gene in self.genes:
                self.conc_list.append([])
                if gene.gene_type == "TF":
                    gene.concentration = 1.0/self.tf_genes
                elif gene.gene_type == "P":
                    gene.concentration = 1.0/self.p_genes
                elif gene.gene_type.startswith("EXTRA"):
                    gene.concentration = 0

        tf_total = sum([i.concentration for i in self.genes
                        if i.gene_type == "TF" ])
        p_total = sum([i.concentration for i in self.genes
                       if i.gene_type == "P"])
        e_total = sum([i.concentration for i in self.genes
                       if i.gene_type.startswith('EXTRA')])

        # normalise concentrations separately
        for idx, gene in enumerate(self.genes):
            if gene.gene_type == "TF":
                gene.concentration *= 1.0 - e_total
                gene.concentration = gene.concentration / tf_total
            elif gene.gene_type == "P":
                gene.concentration = gene.concentration / p_total

def main():
    # random.seed(2)
    # grn = GRN(delta=1)
    # grn.build_genes()
    # grn.add_extra("EXTRA_sineval", 0.0, [0]*32)
    # grn.precalc_matrix()
    # grn.regulate_matrix(5)

    import time
    start_time = time.time()

    for seed in range(10):
        random.seed(seed)
        stabiter = 10000
        runiter = 1000
        grn = GRN(delta=1)

        grn.build_genes()
        grn.add_extra("EXTRA_sineval", 0.0, [0]*32)
        grn.precalc_matrix()
        grn.regulate_matrix(stabiter)

        onff = 0
        for i in range(0, 50):
            if i % 10 == 0:
                if onff == 1:
                    onff = 0
                else:
                    onff = 1

            inputval = onff
            extra_vals = {'sineval': inputval}
            grn.set_extras(extra_vals)
            grn.regulate_matrix(runiter)

        for conc in grn.conc_list:
            print conc[-1]
        filename = "demo"+str(seed)
        graph.plot_2d(grn.conc_list, filename)
    print "took", str(time.time() - start_time)

if __name__ == "__main__":
    main()
