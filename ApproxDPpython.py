import matplotlib.pyplot as plt
import subprocess

cmd = './a.out'
process = subprocess.Popen(cmd, stdout=subprocess.PIPE)

set1 = []
set2 = []

for line in process.stdout:
    temp = (line.decode("ASCII"))[:-1]
    temp = temp.split()
    theta = float(temp[1])
    phi = float(temp[2])
    ptheta = float(temp[3])
    pphi = float(temp[4])
    
    set1.append([theta, ptheta])
    set2.append([phi, pphi])

plt.plot([point[0] for point in set1], [point[1] for point in set1])
plt.plot([point[0] for point in set2], [point[1] for point in set2])
plt.show()

