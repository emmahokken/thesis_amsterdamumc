import csv 
from collections import defaultdict
import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import scipy
from statsmodels.stats.anova import AnovaRM 

import plot_fwhm 

def main():
    # subjects = [5, 8, 18, 25, 31, 64, 77, 98, 105]
    subjects = [5,8]
    
    # gather data
    sharpness = read_data(subjects)
    print(sharpness.columns)
    # caclulate basic stats
    mean_gt = sharpness.sigma_gt.mean()
    mean_rim = sharpness.sigma_rim.mean()
    std_gt = sharpness.sigma_gt.std()
    std_rim = sharpness.sigma_rim.std()

    ventricles = sharpness.loc[(sharpness.fields == 'ventl') | (sharpness.fields =='ventr')]
    plot_fwhm.plot_fwhm_per_subject(sharpness, 5)
    plot_fwhm.plot_per_region(sharpness)

def remove_subj(field):
    return field[4:]


def read_data(subjects):

    step = 3
    begin = 3
    stop = 12
    df = pd.concat([pd.read_csv(f'../../results/{s:03}_{i}_FWHM.csv') for s in subjects for i in range(begin,stop+1,step)], ignore_index=True)
    df.fields = df.fields.apply(remove_subj)    
    
    return df



if __name__ == '__main__':
    main()