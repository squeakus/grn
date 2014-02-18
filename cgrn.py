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
Gene Regulatory Model based on the work of Wolfgang Banzhaf
and Miguel Nicolau.
"""
# Bug when casting from ctype array to nparray
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

import random, math, graph, time
from banzhaf_parser import reorder_gene
import numpy as np    # for matrices
import ctypes as ct   # for grn c-library

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
            self.genome = GRN.random_init(5000)
        else:
            self.genome = genome

        ct.cdll.LoadLibrary('libgrn.so')
        self.grnlib = ct.CDLL('libgrn.so')
        self.promoters = [[[0]*8, "TF"], [[1]*8, "P"]]
        self.logfile = str(time.time())+".out"

        self.sensible_parsing = False
        self.sync_array = []
        self.below_zero = False
        self.restbreak = True
        self.delta = delta
        self.extra_limit = 0.4
        self.genes = []
        self.tf_genes = []
        self.p_genes = []
        self.extras = []
        self.conc_list = []
        self.weight_matrix = None

    @staticmethod
    def random_init(genome_size):
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
        gene_cnt = 0
        max_length = len(self.genome) - 256

        if not self.sensible_parsing:
            index = 96 - len(self.promoters[0][0])
            max_length = len(self.genome) - (160 + len(self.promoters[0][0]))

        while index <= max_length:
            for seq in self.promoters:
                found = self.genome[index:index+len(seq[0])] == seq[0]

                if found and index <= max_length:
                    if seq[1] == "TF":
                        self.tf_genes.append(gene_cnt)
                    elif seq[1] == "P":
                        self.p_genes.append(gene_cnt)

                    if self.sensible_parsing:
                        gene_segment = self.genome[index:index + 256]
                    else:
                        gene_segment = reorder_gene(index,
                                                    seq[0],
                                                    self.genome)
                    self.genes.append(Gene( gene_segment,
                                           seq[1]))
                    index += 255
                    gene_cnt += 1
            index += 1

    def initialise_concentrations(self):
        # TF genes are normalised with Extras, Pgenes are separate
        self.conc_list = []
        e_total = 0
        for idx in self.extras:
            e_total += self.genes[idx].concentration

        for gene in self.genes:
            self.conc_list.append([])
            if gene.gene_type == "TF":
                gene.concentration = (1.0-e_total)/len(self.tf_genes)
            elif gene.gene_type == "P":
                gene.concentration = 1.0/len(self.p_genes)

    def add_extra(self, gene_type, concentration, signature=None):
        """Initialise and add a gene"""
        self.extras.append(len(self.genes))
        segment = [random.randint(0, 1) for _ in range(256)]
        extra_gene = Gene(segment, gene_type)
        extra_gene.concentration = concentration

        if not signature == None:
            extra_gene.protein = signature

        self.genes.append(extra_gene)
        self.conc_list.append([])

    def set_extras(self, extra_vals):
        """Set all EXTRA protein concs"""
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
        Then generate weights by subtracting inhibitor from the enhancer
        """
        enh_matrix = np.zeros((len(self.genes),
                               len(self.genes)-len(self.extras)))
        inh_matrix = np.zeros((len(self.genes),
                               len(self.genes)-len(self.extras)))

        #cheekily initialise concs after extras are added (TF+E)
        self.initialise_concentrations()
        #iterate proteins and genes to create matrix
        for idx, gene in enumerate(self.genes):
            if not gene.gene_type == "P":
                protein = gene.protein

                for idy in range(len(self.genes)-len(self.extras)):
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

    def regulate_matrix(self, iterations, restop=False):
        """ Grind the concs against the weights and normalise.
        Pass concs, weights tfs and p genes to the
        grnlib for fast calculation"""
        weights = self.weight_matrix.astype(np.float64)
        cweights = weights.ctypes.data_as(ct.POINTER(ct.c_double))
        rows, cols = self.weight_matrix.shape

        concs = []
        for gene in self.genes:
            concs.append(gene.concentration)
        e_total = ct.c_double(sum([i.concentration for i in self.genes
                                   if i.gene_type.startswith('EXTRA')]))
        # compute array lengths
        conc_len = len(concs)
        old_concs = np.zeros(conc_len)
        tf_len = len(self.tf_genes)
        p_len = len(self.p_genes)

        # cast to c arrays
        cconcs = (ct.c_double * conc_len)(*concs)
        ctfs = (ct.c_int * tf_len)(*self.tf_genes)
        cps = (ct.c_int * p_len)(*self.p_genes)

        # call grnlib
        for current in range(iterations):
            bzero = self.grnlib.regulate(conc_len, cconcs,
                                         rows, cols, cweights,
                                         tf_len, ctfs, p_len, cps,
                                         self.delta, e_total)
            # quit if system at rest
            if self.restbreak:
                tmp_concs = np.ctypeslib.as_array((ct.c_double * conc_len).from_address(ct.addressof(cconcs)))
                difference = sum(abs(tmp_concs-old_concs))/float(conc_len)
                if difference < 1E-10:
                    break
                old_concs = np.copy(tmp_concs)

            if bzero == 1:
                self.below_zero = True


        outfile = open(self.logfile, 'a')
        outfile.write(str(current)+'\n')
        outfile.close()
        if restop == True:
            self.sync_array.append(current)
        # cast ctypes to numpy and set gene values
        final_concs = np.ctypeslib.as_array(
            (ct.c_double * conc_len).from_address(ct.addressof(cconcs)))

        for idx, gene in enumerate(self.genes):
            gene.concentration = final_concs[idx]
            self.conc_list[idx].append(final_concs[idx])

def main():
    """ Code to compare with eoins grn"""
    random.seed(2)
    grn = GRN(delta=1)
    grn.read_genome("eoinseed3.txt")
    grn.build_genes()
    grn.add_extra("EXTRA_1", 0.05, [0]*32)
    grn.add_extra("EXTRA_2", 0.05, [1]*32)
    grn.add_extra("EXTRA_3", 0.05, [0]*16 +[1]*16)
    grn.add_extra("EXTRA_4", 0.05, [1]*16 +[0]*16)

    grn.precalc_matrix()
    for _ in range(10):
        # opstring = ""
        # for geneid in grn.tf_genes:
        #     tmpstr = "%0.9f " % grn.genes[geneid].concentration
        #     opstring +=  tmpstr
        # for geneid in grn.extras:
        #     tmpstr = "%0.9f " % grn.genes[geneid].concentration
        #     opstring +=  tmpstr
        # for geneid in grn.p_genes:
        #     tmpstr = "%0.9f " % grn.genes[geneid].concentration
        #     opstring +=  tmpstr
        # print opstring
        grn.regulate_matrix(10000)


    # speed test code
    # import time
    # start_time = time.time()
    # for seed in range(10):
    #     random.seed(seed)
    #     stabiter = 10000
    #     runiter = 1000

    #     grn = GRN(delta=1)
    #     grn.build_genes()
    #     grn.add_extra("EXTRA_sineval", 0.0, [0]*32)
    #     grn.precalc_matrix()
    #     grn.regulate_matrix(stabiter)

    #     onff = 0
    #     for i in range(0, 50):
    #         if i % 10 == 0:
    #             if onff == 1:
    #                 onff = 0
    #             else:
    #                 onff = 1

    #         inputval = onff
    #         extra_vals = {'sineval': inputval}
    #         grn.set_extras(extra_vals)
    #         grn.regulate_matrix(runiter)

    #     for conc in grn.conc_list:
    #         print conc[-1]

    #     filename = "sync"+str(seed)
    #     graph.plot_2d(grn.conc_list, filename)
    # print "took", str(time.time() - start_time)

if __name__ == "__main__":
    main()
