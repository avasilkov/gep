from run import *
import random

for i in range(1000000):
    a = int(random.random()*4294967295)
    if from_gray(to_gray(a)) != a:
        print("error ", a)
