import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import subprocess

set_cord = [[], []]
cmd = './output3.o'
process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
i = 0

for line in process.stdout:
    temp = (line.decode("ASCII"))[:-1]
    temp = temp.split()
    set_cord[0].append(float(temp[0])) #phi
    set_cord[1].append(float(temp[1])) #phi_momentum

fig, ax = plt.subplots()
ax.scatter(set_cord[0], set_cord[1])
plt.show()
