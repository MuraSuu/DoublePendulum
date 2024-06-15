import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import subprocess

set_cord = [[], []]
cmd = './daniel.out'
process = subprocess.Popen(cmd, stdout=subprocess.PIPE)

for line in process.stdout:
    temp = (line.decode("ASCII"))[:-1]
    temp = temp.split()
    set_cord[0].append(float(temp[0])) #angle
    set_cord[1].append(float(temp[1])) #action

fig, ax = plt.subplots()
ax.scatter(set_cord[0], set_cord[1])
ax.set_xlabel('θ')
ax.set_xlabel('J')
plt.show()