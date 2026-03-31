import matplotlib.pyplot
import numpy as np
import matplotlib.pyplot as plt




f1 = open("../Results/Nonzeroelements_init.txt")
f2 = open("../Results/Nonzeroelements_first_eval.txt")
f3 = open("../Results/Nonzeroelements_real_eval.txt")


line1 = f1.readline()
line2 = f2.readline()
line3 = f3.readline()

Splitted_line1 = line1.split(' ')
Splitted_line2 = line2.split(' ')
Splitted_line3 = line3.split(' ')


print("The mats have {} {} {} nonzero components ".format(len(Splitted_line1),len(Splitted_line2),len(Splitted_line3)))

x1 = []
y1 = []
x2 = []
y2 = []
x3 = []
y3 = []


for i in range(len(Splitted_line1)):
    Pair = Splitted_line1[i].split(',')
    if(len(Pair)<2):
        continue
    x1.append(float(Pair[0]))
    y1.append(float(Pair[1]))

    Pair = Splitted_line2[i].split(',')
    x2.append(float(Pair[0]))
    y2.append(float(Pair[1]))

    Pair = Splitted_line3[i].split(',')
    x3.append(float(Pair[0]))
    y3.append(float(Pair[1]))




plt.scatter(x3,y3)
# plt.scatter(x2,y2)
# plt.scatter(x1,y1)
plt.axis('equal')
plt.show()