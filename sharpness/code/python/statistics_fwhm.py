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
    gp = df.loc[(df.fields == 'gpl') | (df.fields == 'gpr')]
    df = determine_outliers(df)
    # t_test(df)
    print(df[df.fields == 'thal'])
    # plot_fwhm.plot_fwhm_per_subject(df, 8)
    plot_fwhm.plot_per_region(df.loc[df.fields != 'vent4'], 'gt')
    # plot_fwhm.plot_per_region(df.loc[df.fields != 'vent4'], 'rim')
    plot_fwhm.boxplot(df)
    # for s in subjects:
    #     a = df.loc[(df.fields == 'strr') & (df.subj_id == s)]
    #     print(a)
    #     plt.plot(a.sigma_gt)
    #     plt.show()

    # plt.hist(df.loc[df.acc_factor == 6].sigma_rim, color='rebeccapurple')
    # plt.title('Histogram of FWHM for reconstructed image,s all acceleration rates')
    # plt.show()
    # print(scipy.stats.normaltest(df.loc[df.acc_factor == 6].sigma_gt))
    # rmanova(df.dropna())

def determine_outliers(df):
    fields = df.fields.unique()
    subj = df.subj_id.unique()
    
    for s in subj:
        for f in fields:
            sigma_gt = df.loc[(df.fields == f) & (df.subj_id == s)].sigma_gt
            if len(sigma_gt.unique()) > 1:
                for val in sigma_gt.unique():
                    if len(sigma_gt[sigma_gt == val]) < 3:
                        ind = sigma_gt[sigma_gt == val].index
                        df = df.drop(ind)
                        break 

    return df 
            

def rmanova(df):
    '''
    Perform a repeated measuers ANOVA

    Args:
        df: DataFrame object containing FWHM information
    '''
    df = df.drop(columns=['sigma_gt'])
    # print(df)
    rmanova = AnovaRM(data=df, depvar='sigma_rim', subject='subj_id', within=['acc_factor','fields']).fit()
    print(rmanova)

def t_test(df):
    ''' Perform a t-test on the data. 

    Args:
        df: DataFrame object containing FWHM information
    '''

    globus = df[(df.fields == 'gpl') | (df.fields == 'gpr')]
    globus_3 = globus[globus.acc_factor == 3].sigma_rim
    globus_12 = globus[globus.acc_factor == 12].sigma_rim
    print(globus[globus.acc_factor == 12])
    levene = scipy.stats.levene(globus_3, globus_12)
    print(levene)
    ttest = scipy.stats.ttest_ind(globus_3, globus_12)
    print(ttest)
    linres = scipy.stats.linregress(globus_3, globus_12)
    print(linres)
    # plt.plot(globus_3.sigma_rim, label='3')
    # plt.plot(globus_12.sigma_rim, label='12')
    # plt.legend()
    # plt.show()

def remove_subj(field):
    return field[4:]


def read_data(subjects):

    step = 3
    begin = 3
    stop = 12

    # read files for all subjects and all acceleration factors and concatinate them 
    df = pd.concat([pd.read_csv(f'../../results/{s:03}_{i}_FWHM.csv') for s in subjects for i in range(begin,stop+1,step)], ignore_index=True)
    df.fields = df.fields.apply(remove_subj)    
    
    return df.dropna()



if __name__ == '__main__':
    main()