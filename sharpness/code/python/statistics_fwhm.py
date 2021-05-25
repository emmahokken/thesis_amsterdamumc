import csv 
from collections import defaultdict
import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import scipy
from statsmodels.stats.anova import AnovaRM 

import plot_fwhm 

def main():
    subjects = [5, 8, 18, 25, 31, 64, 77, 98, 105]
    # subjects = [5,8]
    
    # gather data
    df = read_data(subjects)
    print(df.columns)
    # caclulate basic stats
    mean_gt = df.sigma_gt.mean()
    mean_rim = df.sigma_rim.mean()
    std_gt = df.sigma_gt.std()
    std_rim = df.sigma_rim.std()
    
    ventricles = df.loc[(df.fields == 'ventl') | (df.fields =='ventr')]
    # plot_fwhm.plot_fwhm_per_subject(df, 8)
    # plot_fwhm.plot_per_region(df.loc[df.fields != 'vent4'], 'gt')
    # plot_fwhm.plot_per_region(df.loc[df.fields != 'vent4'], 'rim')
    plot_fwhm.boxplot(df)
    # for s in subjects:
    #     a = df.loc[(df.fields == 'strr') & (df.subj_id == s)]
    #     print(a)
    #     plt.plot(a.sigma_gt)
    #     plt.show()

    rmanova(df.dropna())


def rmanova(df):
    '''
    Perform a repeated measuers ANOVA

    Args:
        df: DataFrame object containing FWHM information
    '''
    df = df.drop(columns=['sigma_gt'])
    print(df)
    rmanova = AnovaRM(data=df, depvar='sigma_rim', subject='subj_id', within=['acc_factor','fields']).fit()
    print(rmanova)

def remove_subj(field):
    return field[4:]


def read_data(subjects):

    step = 3
    begin = 3
    stop = 12

    # read files for all subjects and all acceleration factors and concatinate them 
    df = pd.concat([pd.read_csv(f'../../results/{s:03}_{i}_FWHM.csv') for s in subjects for i in range(begin,stop+1,step)], ignore_index=True)
    df.fields = df.fields.apply(remove_subj)    
    
    return df



if __name__ == '__main__':
    main()