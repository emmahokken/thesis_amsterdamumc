import numpy as np
import matplotlib.pyplot as plt 
import csv 
from collections import defaultdict

structure = defaultdict(lambda : {'gt':[], 'rim':[]}) 
begin = 3
stop = 16
step = 3
for i in range(begin, stop, step):
    with open(f'../../008_{i}_FWHM.csv','r') as f:
        reader = csv.reader(f)
        lc = 0
        for row in reader:
            if lc > 0:
                print(row)
                structure[i]['gt'].append(float(row[0]))
                structure[i]['rim'].append(float(row[1]))
            lc += 1

print(structure[15])
for i in range(begin, stop, step):
    acc_factor = [i for j in range(len(structure[i]['gt']))]
    plt.scatter(acc_factor, structure[i]['rim'], c='rebeccapurple', marker='o')
    plt.scatter(acc_factor, structure[i]['gt'], c='orange', marker='^')
plt.legend(['rim','gt'])
plt.xticks(range(3,16,3))
plt.xlabel('Acceleration factor')
plt.ylabel('Sharpness in FWHM')
plt.title('Sharpness as a function of acceleration factor for one subject')
plt.savefig('../../plots_saved/008_FWHM_all.pdf')
plt.show()