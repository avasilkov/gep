import random
import copy

def compute_tail_length(head_len, max_args_number):
    """ computes tail size """
    return head_len*(max_args_number - 1) + 1

def compute_gene_length(head_len, max_args_number):
    """ computes gene length: head + tail """
    return head_len + compute_tail_length(head_len, max_args_number)

#make an additional gene that operates on other genes!
#gradually increase gene/chromosome size and mayb functions complexity p8

def coin_flip_swap(first, second):
    if random.random() < 0.5:
        return second, first
    return first, second

def roulette_wheel(weights, max_weight):
    """ from https://en.wikipedia.org/wiki/Fitness_proportionate_selection
    chooses an index based on weight values as probabilities"""
    n = len(weights)
    limit = 5#check if it's enough
    while limit > 0:
        index = int(n*random.random())#random() in [0,1) and int(1.9)==1
        if random.random() < weights[index]/max_weight:
            break
        limit -= 1
    return index

class Node:

    def __init__(self, code, abc):
        """ code is only needed for graphical repr right now.
        I don't think it's too much overhead to store it here"""
        self.f, self.args_number = abc[code]
        self.code = code
        self.children = []

    def compile_expression(self):
        return self.f.format(*[c.compile_expression() for c in self.children])

    def __repr__(self):
        return self.code

def parse_gene_into_tree(gene, abc):
    abc = abc.funs_and_terms
    layers = []
    layers.append([Node(gene[0], abc)])
    if layers[0][0].args_number == 0:
        return layers
    terminated = False
    pos = 1#len of a geneosome can't be one because then it's just a terminal
    #and zero element is not terminal either.
    #head and tail boundaries are preserved so any geneosome is a valid expression terminated at the
    #end
    while not terminated:
        next_layer = []
        terminated = True
        for node in layers[-1]:
            for j in range(node.args_number):
                next_layer.append(Node(gene[pos], abc))
                node.children.append(next_layer[-1])
                if next_layer[-1].args_number != 0:
                    terminated= False
                pos += 1
        layers.append(next_layer)

    return layers

def get_code_and_args_n_from_layers(layers):
    return [[(n.code, n.args_number) for n in l] for l in layers]


class Alphabet:

    def __init__(self):
        self.functions = {}
        self.terminals = {}
        self.funs_and_terms = None

        self.functions_lst = None
        self.terminals_lst = None
        self.funs_and_terms_lst = None

        self.max_args_number = 0

    def finished_adding_elements(self):
        self.funs_and_terms = {**self.functions, **self.terminals}
        self.functions_lst = list(self.functions.keys())
        self.terminals_lst = list(self.terminals.keys())
        self.funs_and_terms_lst = self.functions_lst + self.terminals_lst

    def add_function(self, code, str_val, args_number):
        assert not code in self.functions, "{0} function has already been added to functions, WARNING".format(code)
        assert not code in self.terminals, "{0} function has already been added to terminals, WARNING".format(code)
        if self.max_args_number < args_number:
            self.max_args_number = args_number
        self.functions[code] = (str_val, args_number)

    def add_terminal(self, code, str_val):
        assert not code in self.terminals, "{0} terminal has already been added to terminals, WARNING".format(code)
        assert not code in self.functions, "{0} terminal has already been added to functions, WARNING".format(code)
        self.terminals[code] = (str_val, 0)



class GeneObj:

    def __init__(self, head_len, max_args_number):
        self.head_len = head_len
        self.tail_len = compute_tail_length(head_len, max_args_number)
        self.gene_len = compute_gene_length(head_len, max_args_number)

class Organism:

    def __init__(self, gene_obj, gene_number, abc, linker=None, random_init=True):
        self.gene_obj = gene_obj
        self.gene_number = gene_number
        self.dna = []
        self.linker = linker
        if random_init:
            self.set_random_dna(abc)

    def set_random_dna(self, abc):
        dna = []
        for g_n in range(self.gene_number):
            for i in range(self.gene_obj.head_len):
                dna.append(random.choice(abc.funs_and_terms_lst))
            for i in range(self.gene_obj.tail_len):
                dna.append(random.choice(abc.terminals_lst))
        self.dna = dna

    def __repr__(self):
        s = ''
        for i in range(self.gene_number):
            for j in range(self.gene_obj.gene_len):
                s += str(j % 10)
            s += ' '
        s += '\n'
        for i in range(self.gene_number):
            s += ''.join(self.dna[i*self.gene_obj.gene_len:(i + 1)*self.gene_obj.gene_len]) + ' '

        return s

    def copy(self):
        return copy.deepcopy(self)

    def mutate(self, abc, m_rate):
        """ TODO calc m_rate based on chromosome size and how many points you want to mutate per
        mutation.
        TODO add tail and head separate rates and/or encode m_rate i dna as a number.
        TODO make different m_rates for mutations that change args_number and drastically change
        expr tree"""
        for g_n in range(self.gene_number):
            offset = g_n*self.gene_obj.gene_len
            for i in range(self.gene_obj.head_len):
                if random.random() < m_rate:
                    self.dna[offset + i] = random.choice(abc.funs_and_terms_lst)
            offset += self.gene_obj.head_len
            for i in range(self.gene_obj.tail_len):
                if random.random() < m_rate:
                    self.dna[offset + i] = random.choice(abc.terminals_lst)


    def insertionIS(self, max_len):
        """ 0.1 is suggested m_rate(checking it outside.
        maybe add random filter that filters less the farther it wants to insert from root"""
        if self.gene_obj.head_len < 2 or max_len > self.gene_obj.gene_len:#can't insert at root and
            #max len need to be of reasonable size
            return
        size = random.randrange(max_len) + 1
        pos = random.randrange(self.gene_number)*self.gene_obj.gene_len
        pos_in_gene = random.randrange(self.gene_obj.head_len - 1) + 1#excluding root
        transposon_pos = random.randrange(len(self.dna))
        seq_len = min(size, self.gene_obj.head_len - pos_in_gene)
        seq = self.dna[transposon_pos:transposon_pos + seq_len]#while its still intact
        seq_len = len(seq)

        self.dna[pos + pos_in_gene + seq_len:pos + self.gene_obj.head_len] = self.dna[pos + pos_in_gene:pos + self.gene_obj.head_len - seq_len]
        self.dna[pos + pos_in_gene:pos + pos_in_gene + seq_len] = seq
        #print(size, pos, pos_in_gene, seq_len, transposon_pos)

    def root_insertionRIS(self, max_len, abc):
        """ 0.1 is suggested m_rate(checking it outside.
        maybe add random filter that filters less the farther it wants to insert from root"""
        if max_len > self.gene_obj.gene_len:
            #max len need to be of reasonable size
            return
        size = random.randrange(max_len) + 1
        pos = random.randrange(self.gene_number)*self.gene_obj.gene_len
        pos_in_head = random.randrange(self.gene_obj.head_len)
        ult_pos = -1
        for i in range(pos + pos_in_head, pos + self.gene_obj.head_len):
            if self.dna[i] in abc.functions:
                ult_pos = i
                break
        if ult_pos == -1:
            return#function not found, do nothing
        pos_in_head = ult_pos - pos
        seq_len = min(size, self.gene_obj.gene_len - pos_in_head)
        seq = self.dna[ult_pos:ult_pos + seq_len]
        seq_len = len(seq)

        self.dna[pos + seq_len:pos + self.gene_obj.head_len] = self.dna[pos:pos + self.gene_obj.head_len - seq_len]
        self.dna[pos:pos + seq_len] = seq
        #print(size, pos, pos_in_gene, seq_len, transposon_pos)

    def gene_transposition(self):
        if self.gene_number < 2:
            return
        gn = random.randrange(self.gene_number - 1) + 1
        gene = self.dna[gn*self.gene_obj.gene_len:(gn + 1)*self.gene_obj.gene_len]
        self.dna[self.gene_obj.gene_len:(gn + 1)*self.gene_obj.gene_len] = self.dna[:gn*self.gene_obj.gene_len]
        self.dna[:self.gene_obj.gene_len] = gene

    @staticmethod
    def gene_transposition_any_place(self):
        #TODO ?
        pass

    @staticmethod
    def crossover_pack(org1, org2, one_p_partition, two_p_partition):
        """ typically sum of all three rates is 0.7 """
        n = random.random()
        org1, org2 = coin_flip_swap(org1, org2)
        if n < one_p_partition:
            return Organism.one_point_recombination(org1, org2)
        if n < one_p_partition + two_p_partition:
            return Organism.two_point_recombination(org1, org2)
        else:
            return Organism.one_gene_recombination(org1, org2)

    @staticmethod
    def one_point_recombination(org1, org2):
        pivot = random.randrange(1, len(org1.dna))
        org1 = org1.copy()
        org1.dna[pivot:] = copy.deepcopy(org2.dna[pivot:])
        return org1

    @staticmethod
    def two_point_recombination_with_params(org1, org2, pivot1, pivot2):
        org1 = org1.copy()
        if random.random() < 0.5:
            org1.dna[pivot1:pivot2] = copy.deepcopy(org2.dna[pivot1:pivot2])
        else:
            org1.dna[:pivot1] = copy.deepcopy(org2.dna[:pivot1])
            org1.dna[pivot2:] = copy.deepcopy(org2.dna[pivot2:])
        return org1

    @staticmethod
    def two_point_recombination(org1, org2):
        pivot1 = random.randrange(1, len(org1.dna) - 1)
        pivot2 = random.randrange(pivot1 + 1, len(org1.dna))
        return Organism.two_point_recombination_with_params(org1, org2, pivot1, pivot2)

    @staticmethod
    def one_gene_recombination(org1, org2):
        if org1.gene_number < 2:
            return
        g = random.randrange(org1.gene_number)*org1.gene_obj.gene_len
        return Organism.two_point_recombination_with_params(org1, org2, g, g + org1.gene_obj.gene_len)



#TODO speciation by a chromosome tree size!!! see NEAT
#p8 see important about fitness function in GEP paper




class MutationRates:

    def __init__(self, simple_m_rate, IS, RIS, gene_transposition,
                 one_p_recomb, two_p_recomb, one_gene_recomb,
                 max_transposition_insertion_seq_len):
        self.m_rate = simple_m_rate
        self.max_transposition_insertion_seq_len = max_transposition_insertion_seq_len
        self.IS = IS
        self.RIS = RIS
        self.gene_transposition = gene_transposition
        self.one_p_recomb = one_p_recomb
        self.two_p_recomb = two_p_recomb
        self.one_gene_recomb = one_gene_recomb
        self.crossover_total = one_p_recomb + two_p_recomb + one_gene_recomb
        self.one_two_recomb_percents = (self.one_p_recomb/self.crossover_total,
                                        self.two_p_recomb/self.crossover_total)


class EvoAlg:

    def __init__(self, orgs, pop_size, elite_partition, fitness_f, iterations, mutation_rates, abc, plotting=False):
        self.orgs = orgs
        self.pop_size = pop_size
        self.fitness_f = fitness_f
        self.iterations = iterations
        self.elite_partition = elite_partition
        self.plot_data = []#best org, avg fitness, max f
        self.plotting = plotting
        self.m_rates = mutation_rates
        self.abc = abc
        self.best = []

    def compute_fitnesses(self):
        return [(org, self.fitness_f(org)) for org in self.orgs]

    def get_elite(self, orgs_with_fit):
        best = sorted(orgs_with_fit, key=lambda x: x[1], reverse=True)
        return [b[0] for b in best[:int(len(best)*self.elite_partition)]]

    def mutate_organism(self, org):
        """ WARNING MODIFIES ORGANISM  VARIABLE"""
        org.mutate(self.abc, self.m_rates.m_rate)
        if random.random() < self.m_rates.IS:
            org.insertionIS(self.m_rates.max_transposition_insertion_seq_len)
        if random.random() < self.m_rates.RIS:
            org.root_insertionRIS(self.m_rates.max_transposition_insertion_seq_len, self.abc)
        if random.random() < self.m_rates.gene_transposition:
            org.gene_transposition()

        return org

    def sexual_reproduction(self, org1, org2):
        return self.mutate_organism(Organism.crossover_pack(org1, org2, *self.m_rates.one_two_recomb_percents))

    def next_gen(self, orgs_with_fit, max_fitness):
        elite = self.get_elite(orgs_with_fit)
        orgs, fitnesses = zip(*orgs_with_fit)#inverse order from actual

        next_gen = elite[:]

        fitnesses = list(fitnesses)
        f_min = min(fitnesses)
        for i in range(len(fitnesses)):
            fitnesses[i] -= f_min - 0.1
        max_fitness -= f_min - 0.1

        for i in range(self.pop_size - len(next_gen)):
            if random.random() < self.m_rates.crossover_total:
                next_gen.append(self.sexual_reproduction(orgs[roulette_wheel(fitnesses, max_fitness)], orgs[roulette_wheel(fitnesses, max_fitness)]))
            else:
                next_gen.append(self.mutate_organism(orgs[roulette_wheel(fitnesses, max_fitness)].copy()))#important. copy()

        self.orgs = next_gen
        return elite[0]


    def run(self):
        for i in range(self.iterations):
            orgs_with_fit = self.compute_fitnesses()
            max_fitness = orgs_with_fit[0][1]
            avg = 0
            for o in orgs_with_fit:
                if o[1] > max_fitness:
                    max_fitness = o[1]
                avg += o[1]
            avg /= len(orgs_with_fit)

            if self.plotting:
                self.plot_data.append((avg, max_fitness))


            best = self.next_gen(orgs_with_fit, max_fitness)
            self.best.append((best, max_fitness))


            print('----------------')
            print("Iteration ", i)
            if self.plotting:
                print("Max: ", self.plot_data[-1][1])
                print("Average: ", self.plot_data[-1][0])

