import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import subprocess

dataset = []
set_cord = [[], [], [], []]
cmd = './output2.o'
process = subprocess.Popen(cmd, stdout=subprocess.PIPE)

for line in process.stdout:
    temp = (line.decode("ASCII"))[:-1]
    if temp == 'End':
        dataset.append(set_cord)
        set_cord = [[], [], [], []]
        continue
    temp = temp.split()
    set_cord[0].append(float(temp[1])) #theta
    set_cord[1].append(float(temp[2])) #momentum theta
    set_cord[2].append(float(temp[3])) #phi
    set_cord[3].append(float(temp[4])) #momentum phi

fig, ax = plt.subplots(1, 2)
ax[0].plot(set_cord[0], set_cord[1])
ax[1].plot(set_cord[2], set_cord[3])
ax[0].set_xlabel('θ')
ax[1].set_xlabel('φ')
ax[0].set_ylabel('pθ')
ax[1].set_ylabel('pφ')
plt.show()
