import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import subprocess

set_cord = [[], []]
dataset = []
cmd = './output3.o'
process = subprocess.Popen(cmd, stdout=subprocess.PIPE)

for line in process.stdout:
    temp = (line.decode("ASCII"))[:-1]
    if temp == "End":
        dataset.append(set_cord)
        set_cord = [[], []]
        continue
    temp = temp.split()
    set_cord[0].append(float(temp[0])) #theta
    set_cord[1].append(float(temp[1])) #theta_momentum

fig, ax = plt.subplots()
ax.set_xlabel('θ')
ax.set_ylabel('pθ')

N = len(dataset)
for i in range(N):
    ax.scatter(dataset[i][0], dataset[i][1])

plt.show()
