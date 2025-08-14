import numpy as np 
import matplotlib.pyplot as plt 




import csv

f = open("../../../input/L2_harmonic_mesh.obj", "w")
counter = 1
with open('../../../input/coordinates_config_1.csv', mode ='r')as file:
  csvFile = csv.reader(file)
  for lines in csvFile:
        # print(lines)
        f.write("v {} {} {}\n".format(lines[0],lines[1],lines[2]))

with open('../../../input/faces_config_1.csv', mode ='r')as file:
  csvFile = csv.reader(file)
  for lines in csvFile:
        # print(lines)
        f.write("f {} {} {}\n".format(int(lines[0])+1,int(lines[1])+1,int(lines[2])+1))




f.close()