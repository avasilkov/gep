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

weights = [0.4, 0.3, 1.2, 0.1]
max_weight = max(weights)
import collections
d = collections.defaultdict(int)
for i in range(1000):
    d[roulette_wheel(weights, max_weight)] += 1
print(d)


s = "Q*+-abcd"
ops = {"*":("%s*%s", 2),"-":("%s-%s", 2),"+":("%s+%s", 2),"/":("%s/%s", 2), "Q":("%s**0.5", 1)}

class Node:

    def __init__(self, operator, ops):
        self.op_s, self.operands_number = ops[operator]
        self.operands = []

    def to_string(self):
        pass
        #return to_string(self.children) + [child.to_string() for child in children]


added_operators = 0
tree = []
if s[0] in ops:
    tree = ;
row_terminated = False
filled_operators = 0
prev_row_st = 0
prev_row_end = 0
for c in s:
    if c in ops:
        tree.append(Node(c))
    if filled_operators == 0:
        break


