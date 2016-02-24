import _GEP as gep
import random
import collections
import copy

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
chro1 = gep.Chromosome(gene_obj, 3, abc)
print('------mutation------')
print(chro1)
chro1.mutate(abc, 0.9)
print(chro1)
print('------insertionIS/RIS max 3------')
chro2 = chro1.copy()
for i in range(100):
    chro2.insertionIS(3)
    chro2.root_insertionRIS(3, abc)
    assert len(chro2.dna) == chro2.gene_obj.gene_len*chro2.gene_number, "error in IS insertion"
print(chro1)
chro1.root_insertionRIS(3, abc)
print(chro1)
print('------gene transposition------')
print(chro1)
chro1.gene_transposition()
print(chro1)
print('------recombinations------')
print(chro1)
print(chro2)
print('------res one point recombination------')
chro1_dna = copy.deepcopy(chro1.dna)
chro2_dna = copy.deepcopy(chro2.dna)
chro3 = chro1.one_point_recombination(chro2)
assert chro1.dna == chro1_dna
assert chro2.dna == chro2_dna
print(chro3)
print('------res two point recombination------')
chro1_dna = copy.deepcopy(chro1.dna)
chro2_dna = copy.deepcopy(chro2.dna)
chro3 = chro1.two_point_recombination(chro2)
assert chro1.dna == chro1_dna
assert chro2.dna == chro2_dna
print(chro3)
print('------res one gene recombination------')
chro1_dna = copy.deepcopy(chro1.dna)
chro2_dna = copy.deepcopy(chro2.dna)
chro3 = chro1.one_gene_recombination(chro2)
assert chro1.dna == chro1_dna
assert chro2.dna == chro2_dna
print(chro3)
print(chro1.crossover_pack(chro2, 0.5, 0.5, 0.5))
