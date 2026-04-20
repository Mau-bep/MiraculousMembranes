import matplotlib.pyplot
import numpy as np
import matplotlib.pyplot as plt




f1 = open("../Results/Nonzeroelements_init_2.txt")
f2 = open("../Results/Nonzeroelements_first_eval_2.txt")
f3 = open("../Results/Nonzeroelements_real_eval_2.txt")
f4 = open("../Results/Nonzeroelements_init.txt")

f5 = open("../Results/Hessian_structure.txt")
f6 = open("../Results/hessian_debug.txt")

line1 = f1.readline()
line2 = f2.readline()
line3 = f3.readline()
line4 = f4.readline()


Splitted_line1 = line1.split(' ')
Splitted_line2 = line2.split(' ')
Splitted_line3 = line3.split(' ')
Splitted_line4 = line4.split(' ')

print("The mats have {} {} {} nonzero components ".format(len(Splitted_line1),len(Splitted_line2),len(Splitted_line3)))

x1 = []
y1 = []
x2 = []
y2 = []
x3 = []
y3 = []
x4 = []
y4 = []






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


for i in range(len(Splitted_line4)):
    Pair = Splitted_line4[i].split(',')
    if(len(Pair)<2):
        continue
    x4.append(float(Pair[0]))
    y4.append(float(Pair[1]))



# plt.scatter(x3,y3)
plt.scatter(x2,y2,s=0.2)
# plt.scatter(x1,y1)
plt.axis('equal')
# plt.show()


plt.hist(x1,bins= int(np.max(x1)-np.min(x1)))
plt.hist(x2,bins= int(np.max(x2)-np.min(x2)))
plt.hist(x3,bins= int(np.max(x3)-np.min(x3)))

# plt.show()

plt.scatter(x4,y4,s=0.2)
plt.axis('equal')   
# plt.show()
plt.clf()





x5 = []
y5 = []
x6 = []
y6 = []


for line in f5:
    Pair = line.split(',')
    if(len(Pair)<2):
        continue
    x5.append(float(Pair[0]))
    y5.append(float(Pair[1]))

for line in f6:
    Pair = line.split(',')
    if(len(Pair)<2):
        continue
    x6.append(float(Pair[0]))
    y6.append(float(Pair[1]))

# plt.axis('equal')
plt.scatter(x6,y6,s=0.4,label="Hessian debug",color='red')
plt.axis('equal')
plt.scatter(x5,y5,s=0.8,label="Hessian structure",color='black')
print("The hessian structure has {} nonzero components and the hessian debug has {} nonzero components".format(len(x5),len(x6)))
plt.legend()
plt.show()


# plt.plot(x5,'o',label="Hessian structure",color='black')
# plt.plot(x6,'o',label="Hessian debug",color='red')
plt.legend()
plt.show()