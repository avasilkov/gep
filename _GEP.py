import random

def compute_tail_length(head_len, max_args_number):
    """ computes tail size """
    return head_len*(max_args_number - 1) + 1

def compute_gene_length(head_len, max_args_number):
    """ computes gene length: head + tail """
    return head_len + compute_tail_length(head_len, max_args_number)

#make an additional gene that operates on other genes!
#gradually increase gene/chromosome size and mayb functions complexity p8

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

class Chromosome:

    def __init__(self, gene_obj, gene_number, abc, random_init=True):
        self.gene_obj = gene_obj
        self.gene_number = gene_number
        self.dna = ''
        if random_init:
            self.set_random_dna(abc)

    def set_random_dna(self, abc):
        dna = ''
        for g_n in range(self.gene_number):
            for i in range(self.gene_obj.head_len):
                dna += random.choice(abc.funs_and_terms_lst)
            for i in range(self.gene_obj.tail_len):
                dna += random.choice(abc.terminals_lst)
        self.dna = dna

    def __repr__(self):
        s = ''
        for i in range(self.gene_number):
            for j in range(self.gene_obj.gene_len):
                s += str(j % 10)
            s += ' '
        s += '\n'
        for i in range(self.gene_number):
            s += self.dna[i*self.gene_obj.gene_len:(i + 1)*self.gene_obj.gene_len] + ' '

        return s

#TODO speciation by a chromosome tree size!!! see NEAT
#p8 see important about fitness function in GEP paper
