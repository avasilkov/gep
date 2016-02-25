import _GEP as gep
import random
import collections
import copy
from _GEP import Organism
import matplotlib.pyplot as plt
import numpy as np

print(gep.compute_tail_length(10, 2) == 11)
print(gep.compute_gene_length(10, 2) == 21)

weights = [0.4, 0.3, 1.2, 0.1]
max_weight = max(weights)
d = collections.defaultdict(int)
for i in range(1000):
    d[gep.roulette_wheel(weights, max_weight)] += 1
print([q[1] for q in sorted(d.items(), key=lambda x : x[0])], " random test. need something like ", [216, 135, 595, 54])




def get_alphabet():
    abc = gep.Alphabet()

    abc.add_function("*", "({0})*({1})", 2)
    abc.add_function("-", "({0})-({1})", 2)
    abc.add_function("+", "({0})+({1})", 2)
    abc.add_function("/", "({0})/({1})", 2)
    abc.add_function("Q", "({0})**0.5", 1)

    abc.add_terminal("a", "a")
    abc.add_terminal("b", "b")
    abc.add_terminal("c", "c")
    abc.add_terminal("d", "d")

    abc.finished_adding_elements()

    return abc

abc = get_alphabet()
gene = "Q*+-abcd"
expr_tree_layers = gep.parse_gene_into_tree(gene, abc)
print(expr_tree_layers)
expr_tree_layers_fo_gui = gep.get_code_and_args_n_from_layers(expr_tree_layers)
print(expr_tree_layers_fo_gui)
compiled_expr = expr_tree_layers[0][0].compile_expression()
print(compiled_expr)
print(eval(compiled_expr, globals(), {'a':1, 'b':2, 'c':4, 'd':2})==((1+2)*(4-2))**0.5)
gene_obj = gep.GeneObj(10, abc.max_args_number)
org1 = gep.Organism(gene_obj, 3, abc)
print('------mutation------')
print(org1)
org1.mutate(abc, 0.9)
print(org1)
print('------insertionIS/RIS max 3------')
org2 = org1.copy()
for i in range(100):
    org2.insertionIS(3)
    org2.root_insertionRIS(3, abc)
    assert len(org2.dna) == org2.gene_obj.gene_len*org2.gene_number, "error in IS insertion"
print(org1)
org1.root_insertionRIS(3, abc)
print(org1)
print('------gene transposition------')
print(org1)
org1.gene_transposition()
print(org1)
print('------recombinations------')
print(org1)
print(org2)
print('------res one point recombination------')
org1_dna = copy.deepcopy(org1.dna)
org2_dna = copy.deepcopy(org2.dna)
org3 = Organism.one_point_recombination(org1, org2)
assert org1.dna == org1_dna
assert org2.dna == org2_dna
print(org3)
print('------res two point recombination------')
org1_dna = copy.deepcopy(org1.dna)
org2_dna = copy.deepcopy(org2.dna)
org3 = Organism.two_point_recombination(org1, org2)
assert org1.dna == org1_dna
assert org2.dna == org2_dna
print(org3)
print('------res one gene recombination------')
org1_dna = copy.deepcopy(org1.dna)
org2_dna = copy.deepcopy(org2.dna)
org3 = Organism.one_gene_recombination(org1, org2)
assert org1.dna == org1_dna
assert org2.dna == org2_dna
print(org3)
print(Organism.crossover_pack(org1, org2, 0.5, 0.3))
m_rates = gep.MutationRates(0.01, 0.1, 0.1, 0.01, 0.4, 0.15, 0.15, 3)
print(vars(m_rates))

def get_alphabet2():
    abc = gep.Alphabet()

    abc.add_function("*", "({0})*({1})", 2)
    abc.add_function("-", "({0})-({1})", 2)
    abc.add_function("+", "({0})+({1})", 2)

    abc.add_terminal("a", "a")
    abc.add_terminal("b", "b")
    abc.add_terminal("c", "c")
    abc.add_terminal("d", "d")

    abc.finished_adding_elements()

    return abc

abc2 = get_alphabet2()

gene_obj2 = gep.GeneObj(5, abc2.max_args_number)
pop_size = 20
orgs = [Organism(gene_obj2, 1, abc2) for _ in range(pop_size)]

def fff(org):
    expr_tree_layers = gep.parse_gene_into_tree(org.dna, abc2)
    compiled_expr = expr_tree_layers[0][0].compile_expression()
    return -abs(eval(compiled_expr, globals(), {'a':1, 'b':2, 'c':4, 'd':2}) - 5)
import os
os.system('clear')
ea = gep.EvoAlg(orgs, pop_size, 0.2, fff, 100, m_rates, abc2, True)

ea.run()
print(ea.best[-1])
plot_data = list(zip(*ea.plot_data))
plt.plot(plot_data[1], 'r')
plt.plot(plot_data[0], 'g')
plt.show()


