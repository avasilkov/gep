import random
import _GEP as gep
import collections

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
chro = gep.Chromosome(gene_obj, 1, abc)
print(chro)
