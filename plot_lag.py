import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import subprocess

dataset = []
set_cord = [[], [], [], []]
cmd = './output.o'
process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
i = 0

for line in process.stdout:
    temp = (line.decode("ASCII"))[:-1]
    if temp == 'End':
        i+=1
        dataset.append(set_cord)
        set_cord = [[], [], [], []]
        continue
    temp = temp.split()
    set_cord[0].append(float(temp[1])) #theta
    set_cord[1].append(float(temp[2])) #theta dot
    set_cord[2].append(float(temp[3])) #phi
    set_cord[3].append(float(temp[4])) #phi dot

fig, ax = plt.subplots(1, 2)
N = len(dataset)
ax[0].plot(set_cord[0], set_cord[1])
ax[1].plot(set_cord[2], set_cord[3])
plt.show()
