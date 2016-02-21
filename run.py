#generate and then calculate how far - multiplier
import matplotlib.pyplot as plt
import copy
import random
import numpy as np

MAX_INT = 4294967295.0

def to_gray(x):
    return x ^ (x >> 1)

def from_gray(x):
    x ^= x>>16
    x ^= x>>8
    x ^= x>>4
    x ^= x>>2
    x ^= x>>1
    return x

class SimpleFloatGene:

    def __init__(self, x=0):
        self.bits = to_gray(int(x*MAX_INT))

    def get_number(self):
        return from_gray(self.bits)/MAX_INT

    def mutate(self, m_rate):
        for bit in range(32):
            if random.random() < m_rate:
                self.bits ^= 1 << bit

increase number of chromosomes!
http://www.gene-expression-programming.com/GepBook/Introduction.htm
class SimpleDNA:

    def __init__(self, genes_in_numbers):
        self.genes = [SimpleFloatGene(q) for q in genes_in_numbers]

    def crossover(self, this, other):
        """ accepts deep copies """
        pivot = round(random.random()*len(this.genes))
        dna1 = this
        dna2 = other
        if random.random() > 0.5:
            for i in range(pivot, len(this.genes)):
                dna1.genes[i] = dna2.genes[i]
            return dna1
        else:
            for i in range(pivot):
                dna2.genes[i] = dna1.genes[i]
            return dna2

    def mutate(self, m_rate):
        for gene in self.genes:
            gene.mutate(m_rate)


class GenAlg:

    def __init__(self, orgs, pop_size, fitness_f, iterations, choose_best_per=0.1,
            m_rate=0.01, crossover_m_rate=0.8, plotting=False):
        self.orgs = orgs
        self.pop_size = pop_size
        self.fitness_f = fitness_f
        self.iterations = iterations
        self.plot_data = [[], []]#best org, avg fitness
        self.plotting = plotting
        self.choose_best_per = choose_best_per
        self.m_rate = m_rate
        self.crossover_m_rate = crossover_m_rate

    def choose_best(self):
        best = []
        for org in self.orgs:
            best.append((org, self.fitness_f(org)))

        best.sort(key=lambda x: x[1], reverse=True)

        if self.plotting:
            self.plot_data[0].append(copy.deepcopy(best[0][0]))
            self.plot_data[1].append(sum([f[1] for f in best])/len(best))

        return [b[0] for b in best[:int(round(len(best)*self.choose_best_per))]]

    def reproduce(self, first_parents, second_parents):
        children = []
        need_more = self.pop_size
        while need_more > 0:
            iterate_another_n = min(need_more, len(first_parents))
            need_more -= iterate_another_n
            for i in range(iterate_another_n):
                parent1 = copy.deepcopy(first_parents[i])
                if random.random() > self.crossover_m_rate:#yes > this sign not < we are testing for false
                    children.append(parent1)
                    continue
                parent2 = copy.deepcopy(random.choice(second_parents))
                parent1.mutate(self.m_rate)
                parent2.mutate(self.m_rate)
                child = parent1.crossover(parent1, parent2)
                children.append(child)

        self.orgs = children


    def iteration(self):
        best = self.choose_best()
        self.reproduce(best, self.orgs)

    def run(self):
        for i in range(self.iterations):
            self.iteration()
            print("Iteration ", i)
            if self.plotting:
                print("Average: ", self.plot_data[1][-1])


genes_n = 50
init_val = 0.0
pop_size = 10
m_rate = 0.001
crossover_m_rate = 0.9
iters = 1000
best_per = 0.2
orgs = [SimpleDNA([init_val]*genes_n) for i in range(pop_size)]
#orgs, pop_size, fitness_f, iterations, choose_best_per=0.1, m_rate=0.01, crossover_m_rate=0.8, plotting=False):

def fit1(org):
    t = 0
    for g in org.genes:
        t += g.get_number()
    return t

genAlg = GenAlg(orgs, pop_size, fit1, iters, best_per, m_rate, crossover_m_rate, True)
genAlg.run()

plt.plot(genAlg.plot_data[1])
plt.ylabel('Average 0..%d' % genes_n)
plt.show()
plt.imshow(np.transpose([[g.get_number() for g in org.genes] for org in genAlg.plot_data[0]]), cmap=plt.cm.gray, aspect='auto', interpolation='nearest')
plt.show()


#AIM
# "1/(1+exp(-(x)))"

there are operators, numbers, 1, input_number.
crossover at equal substrings.
insertions...deletions?...
Insert and delete parenthesis at the same time.
"(1)/(1+e(-(x)))"
can only delete parenthesis as a whole?

make sigma neuron a unit in the future and evolve whole cluster to achieve an LStM cell?
Fitness functions based on the correlation coefficient and R-square are also very smooth. For
regression problems, these functions work best by combining them with other measures because, by
themselves, they only tend to measure correlation, not caring for the range of values of the model
output. So by combining them with functions that work at approximating the range of the target
values, they form very efficient fitness functions for finding models with good correlation and good
fit between predicted and actual values.

So by counting the TP, TN, FP, and FN and further assigning different weights to these four types of
classifications, it is possible to create smoother and therefore more efficient fitness functions.
Some popular fitness functions based on the confusion matrix include sensitivity/specificity,
recall/precision, F-measure, Jaccard similarity, Matthews correlation coefficient, and cost/gain
matrix which combines the costs and gains assigned to the 4 different types of classifications.

stochastic Acceptence! https://en.wikipedia.org/wiki/Fitness_proportionate_selection
http://www.gene-expression-programming.com/GepBook/Introduction.htm
http://www.gene-expression-programming.com/GepBook/Introduction.htm
http://www.gene-expression-programming.com/GepBook/Introduction.htm
