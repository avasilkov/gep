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




def get_abc():
    abc = {"*":("(%s)*(%s)", 2),"-":("(%s)-(%s)", 2),"+":("(%s)+(%s)", 2),"/":("(%s)/(%s)", 2), "Q":("(%s)**0.5", 1),
            "a":("a", 0), "b":("b", 0), "c":("c", 0), "d":("d", 0)}
    return abc
abc = get_abc()
chrom = "Q*+-abcd"
print(gep.parse_gene_into_tree(chrom, abc))
print(gep.get_code_and_args_n_from_layers(gep.parse_gene_into_tree(chrom, abc)))

