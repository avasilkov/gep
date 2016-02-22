import random

def compute_tail_length(head_len, max_operands_number):
    """ computes tail size """
    return head_len*(max_operands_number - 1) + 1

def compute_gene_length(head_len, max_operands_number):
    """ computes gene length: head + tail """
    return head_len + compute_tail_length(head_len, max_operands_number)

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

    def to_string(self):
        pass
        #return to_string(self.children) + [child.to_string() for child in children]

    def __repr__(self):
        return self.code

def parse_gene_into_tree(chrom, abc):
    layers = []
    layers.append([Node(chrom[0], abc)])
    if layers[0][0].args_number == 0:
        return layers
    current_layer_size = layers[0][0].args_number
    pos = 1#len of a chromosome can't be one because then it's just a terminal
    #and zero element is not terminal either.
    #head and tail boundaries are preserved so any chromosome is a valid expression terminated at the
    #end
    while current_layer_size > 0:
        next_layer = []
        next_layer_size = 0
        for i in range(current_layer_size):
            next_layer.append(Node(chrom[pos], abc))
            next_layer_size += next_layer[-1].args_number
            pos += 1
        layers.append(next_layer)
        current_layer_size = next_layer_size

    return layers

def get_code_and_args_n_from_layers(layers):
    return [[(n.code, n.args_number) for n in l] for l in layers]


