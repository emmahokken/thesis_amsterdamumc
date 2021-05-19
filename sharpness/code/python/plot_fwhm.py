import csv 
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt 

def plot_fwhm(subj, begin, stop):
    '''
    Generates a scatter plot of the calculated FWHM (sharpness). 

    Args:
        subj: subject to investigate
        begin: starting acceleration rate (factor of 3)
        end: ending acceleration rate (factor of 3, < 15)
    '''
    
    # define structure in which the data will be temporarily stored 
    structure = defaultdict(lambda : {'name':[], 'gt':[], 'rim':[]}) 
    step = 3
    
    for i in range(begin, stop + 1, step):
        with open(f'../../{subj:03}_{i}_FWHM.csv','r') as f:
            reader = csv.reader(f)
            lc = 0
            for row in reader:
                if lc > 0:
                    # print(row)
                    structure[i]['name'].append(row[0][4s:])
                    structure[i]['gt'].append(float(row[1]))
                    structure[i]['rim'].append(float(row[2]))
                lc += 1

    for i in range(begin, stop, step):
        acc_factor = [i for j in range(len(structure[i]['gt']))]
        plt.scatter(acc_factor, structure[i]['rim'], c='rebeccapurple', marker='o')
        plt.scatter(acc_factor, structure[i]['gt'], c='orange', marker='^')
        for j in range(len(structure[i]['gt'])):
            plt.annotate(structure[i]['name'][j], (acc_factor[j], structure[i]['rim'][j]))
            plt.annotate(structure[i]['name'][j], (acc_factor[j], structure[i]['gt'][j]))
    plt.legend(['rim','gt'])
    plt.xticks(range(3,16,3))
    plt.xlabel('Acceleration factor')
    plt.ylabel('Sharpness in FWHM')
    plt.title('Sharpness as a function of acceleration factor for one subject')
    plt.savefig(f'../../plots_saved/{subj:03}_FWHM_all.pdf')
    plt.show()

if __name__ == '__main__':

    subj = 25
    begin = 3
    stop = 12

    plot_fwhm(subj, begin, stop)