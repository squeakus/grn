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

"""Simple Evolutionary strategy algorithm"""
import random, graph

class Evostrategy:
    """A population based evolutionary strategy using the 1/5th rule"""
    def __init__(self, genome_size, pop_size):
        self.genome_size = genome_size
        self.pop_size = pop_size
        self.pop = []
        self.initialise_pop()
        self.gen_count = 0
        self.mut_rate = 0.005
        self.success_mut = 0
        self.generational = True
        self.elite_size = 1
        self.adaptive = False
        self.target = [1 for _ in range(0, self.genome_size)]

    def initialise_pop(self):
        for _ in range(self.pop_size):
            self.pop.append({'genome':self.create_indiv(),
                             'fitness':0,
                             'testfit':0})

    def create_indiv(self):
        """initialise a binary individual"""
        indiv = [random.randint(0, 1) for _ in range(0, self.genome_size)]
        return indiv

    def mutate(self,indiv):
        """mutate each codon with a given probability"""
        mutant = list(indiv)
        for i in range(len(mutant)):
            if random.random() < self.mut_rate:
                mutant[i] = random.randint(0, 1)
        return mutant

    def adapt_mutation(self):
        success_rate = float(self.success_mut) / self.pop_size
        # 1/5th rule
        if success_rate > 0.2:
            self.mut_rate = self.mut_rate * 2
        else:
            self.mut_rate = self.mut_rate / 2
        self.success_mut = 0
        # make sure self.mut_rate is within boundaries
        if self.mut_rate < 0.0005: self.mut_rate = 0.0005
        if self.mut_rate > 0.5: self.mut_rate = 0.5

    def onemax_fitness(self, indiv):
        """sum the number of ones in the genome"""
        fitness = sum(indiv[i] == self.target[i]
                      for i in range(0, len(self.target)))
        return fitness

    def steady_state_replacment(self, children):
        for child in children:
            for idx, parent in enumerate(self.pop):
                if parent['fitness'] < child['fitness']:
                    self.success_mut += 1
                    parent['genome'] = child['genome']
                    parent['fitness'] = child['fitness']
                    parent['testfit'] = child['testfit']
                    break
        return children

    def generational_replacement(self, children):
        for idx in range(0, len(self.pop)-self.elite_size):
            self.pop[idx]['fitness'] = children[idx]['fitness']
            self.pop[idx]['genome'] = children[idx]['genome']
            self.pop[idx]['testfit'] = children[idx]['testfit']
        return children

    def iterate(self, children):
        # sort the parent and child pops
        self.gen_count += 1
        self.pop.sort(key=lambda k: k['fitness'])
        children.sort(key=lambda k: k['fitness'])
        children.reverse()

        # choose replacement strategy
        if self.generational:
            children = self.generational_replacement(children)
        else:
            children = self.steady_state_replacment(children)

        # re-sort pop and mutate
        self.pop.sort(key=lambda k: k['fitness'])
        children = []
        for parent in self.pop:
            mutated_genome = self.mutate(parent['genome'])
            children.append({'genome':mutated_genome,
                             'fitness':0,
                             'testfit':0})
        return children

def main():
    run_list = []

    for run in range(30):
        print "run ", run
        evo = Evostrategy(5000, 100)
        children = evo.iterate(evo.pop)
        best_list = []

        for i in range(100):

            for child in children:
                child['fitness'] = evo.onemax_fitness(child['genome'])
            children = evo.iterate(children)


            if evo.adaptive:
                evo.adapt_mutation()

            best_list.append(evo.pop[-1]['fitness'])
        run_list.append(best_list)

    graph.plot_ave(run_list, "steadystate")
if __name__ == "__main__":
    main()
