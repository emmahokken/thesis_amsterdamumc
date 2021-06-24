import csv 
from collections import defaultdict
import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import scipy
from statsmodels.stats.anova import AnovaRM 
import seaborn as sns 
from sklearn.model_selection import train_test_split
from sklearn import linear_model
import statsmodels.api as sm 
import statsmodels.formula.api as smf 

from plot_fwhm import *

def main():
    subjects = [5, 8, 18, 25, 31, 64, 77, 98, 105]
    # subjects = [5,8]
    
    # gather data
    df = read_data(subjects)

    # caclulate basic stats
    mean_gt = df.sigma_gt.mean()
    mean_rim = df.sigma_rim.mean()
    std_gt = df.sigma_gt.std()
    std_rim = df.sigma_rim.std()
    
    ventricles = df.loc[(df.fields == 'ventl') | (df.fields =='ventr')]
    gp = df.loc[(df.fields == 'gpl') | (df.fields == 'gpr')]
    df = determine_outliers(df)

    # df = merge_gt(df)
    t_test(df)
    # linreg(df)


    # plot_fwhm_per_subject(df, 8)
    # plot_per_region(df.loc[df.fields != 'vent4'], 'gt')
    # plot_per_region(df.loc[df.fields != 'vent4'], 'rim')
    boxplot(df)
    # scatter_fwhm(df)
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
                        # ind = sigma_gt[sigma_gt == val].index
                        df = df.drop(sigma_gt.index)
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
    print('\n\nStatistics\n')
    # df = df[df.subj_id != 105]
    df = df[df.subj_id != 18]
    globus_l = df[df.fields == 'gpl']
    globus_r = df[df.fields == 'gpr']
    print(globus_l.shape, globus_r.shape)
    gp = pd.DataFrame({'gpl':np.array(globus_l.sigma_rim), 'gpr': np.array(globus_r.sigma_rim), 'gp':np.mean(np.array([globus_l.sigma_rim, globus_r.sigma_rim]), axis=0), 'acc_factor':np.array(globus_l.acc_factor), 'subj_id':np.array(globus_l.subj_id)})
    tha_l = df[(df.fields == 'thal') & (df.subj_id != 105)]
    tha_r = df[df.fields == 'thar']
    tha = pd.DataFrame({'thal':np.array(tha_l.sigma_rim), 'thar': np.array(tha_r.sigma_rim), 'tha':np.mean(np.array([tha_l.sigma_rim, tha_r.sigma_rim]), axis=0), 'acc_factor':np.array(tha_l.acc_factor), 'subj_id':np.array(tha_l.subj_id)})
    vent_l = df[(df.fields == 'ventl') & (df.subj_id != 105) & (df.subj_id != 25) & (df.subj_id != 5)]
    vent_r = df[(df.fields == 'ventr') & (df.subj_id != 8)]
    print(vent_l, vent_r)
    vent = pd.DataFrame({'ventl':np.array(vent_l.sigma_rim), 'ventr': np.array(vent_r.sigma_rim), 'vent':np.mean(np.array([vent_l.sigma_rim, vent_r.sigma_rim]), axis=0), 'acc_factor':np.array(vent_l.acc_factor), 'subj_id':np.array(vent_l.subj_id)})
    str_l = df[(df.fields == 'strl') & (df.subj_id != 105)]
    str_r = df[df.fields == 'strr']
    print(str_l, str_r)
    stri = pd.DataFrame({'strl':np.array(str_l.sigma_rim), 'strr': np.array(str_r.sigma_rim), 'str':np.mean(np.array([str_l.sigma_rim, str_r.sigma_rim]), axis=0), 'acc_factor':np.array(str_l.acc_factor), 'subj_id':np.array(str_l.subj_id)})
    gp = gp.merge(tha, on=['subj_id','acc_factor'], how='outer')
    gp = gp.merge(vent, on=['subj_id','acc_factor'], how='outer')
    globus = gp.merge(stri, on=['subj_id','acc_factor'], how='outer')

    
    sigmas_3 = np.array([globus[globus.acc_factor == 3][f'{region}'].dropna().to_list() for region in ['gp', 'vent', 'tha', 'str']])
    sigmas_12 = np.array([globus[globus.acc_factor == 12][f'{region}'].dropna().to_list() for region in ['gp', 'vent', 'tha', 'str']])
    globus_3 = globus[globus.acc_factor == 3].gp
    globus_12 = globus[globus.acc_factor == 12].gp

    sigmas_3 = [item for sublist in sigmas_3 for item in sublist]
    sigmas_12 = [item for sublist in sigmas_12 for item in sublist]

    # exit()
    normal = scipy.stats.normaltest(sigmas_3)
    print(normal)
    normal = scipy.stats.normaltest(sigmas_12)
    print(normal)
    levene = scipy.stats.levene(sigmas_3, sigmas_12)
    print(levene)
    ttest = scipy.stats.ttest_rel(sigmas_3, sigmas_12)
    print('t_test')
    print(ttest)
    slope = []
    # linres = scipy.stats.linregress(acc_factors, sigma)
    # print(linres)
    for subj in globus.subj_id.unique():
        for region in ['gp', 'vent', 'tha', 'str']:
            globus_subj = globus[globus.subj_id == subj]
            globus_subj = globus_subj[(globus_subj.acc_factor == 3) | (globus_subj.acc_factor == 12)]

            fwhm = globus_subj[f'{region}'].dropna()
            if len(fwhm) > 0:
                acc_factors = globus_subj.acc_factor[:len(fwhm)]
                b, m = np.polynomial.polynomial.polyfit(acc_factors, fwhm, 1)
                slope.append(m)

        # plt.plot(acc_factors,gp, 'o')
        # plt.plot(acc_factors, acc_factors*m + b)
        # plt.show()
        
        # t_test_subj = scipy.stats.ttest_rel(globus_subj[globus_subj.acc_factor == 3].gp, globus_subj[globus_subj.acc_factor == 12].gp)
        # print(t_test_subj)
    slope_mean = np.mean(slope)
    # slope = slope +  slope
    print(slope_mean)
    one_sided = scipy.stats.ttest_1samp(slope, 0.0)
    print(one_sided)
    x = np.array([3, 6, 9, 12])
    print(x*slope_mean)
    plt.plot(x, x*slope_mean)
    plt.show()    

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

def merge_gt(df):
    '''
    Merges the ground truth values into the RIM values DataFrame
    
    Args:
        df: DataFrame object containing FWHM information
    '''

    def return_one(x):
        return 1

    sigma_rim = df.drop(columns=['sigma_gt'])
    sigma_gt = df[df.acc_factor == 3].drop(columns=['sigma_rim'])
    sigma_gt.acc_factor = sigma_gt.acc_factor.apply(return_one)
    sigma_gt = sigma_gt.rename(columns={'sigma_gt': 'sigma_rim'})
    df = pd.concat([sigma_gt, sigma_rim])
    print(df)
    return df 

def linreg(df):
    '''
    Perform fit a linear regression model to the data. 

    Sources: https://medium.com/analytics-vidhya/implementing-linear-regression-using-sklearn-76264a3c073c 
             https://medium.com/swlh/interpreting-linear-regression-through-statsmodels-summary-4796d359035a 
    '''

    df.sigma_rim.hist(color='rebeccapurple')
    plt.show()

    sigma = df.sigma_rim
    acc_factors = df.acc_factor
    # acc_factors = df.acc_factor.apply(lambda x: np.log(x))
    df.acc_factor = df.acc_factor.apply(lambda x: str(x))
    x = pd.get_dummies(data=df['acc_factor'])
    y = df.sigma_rim 

    x_train, x_test, y_train, y_test = train_test_split(x,y,test_size=0.2, random_state=42)
    # print(x_train.shape, y_train.shape, x_test.shape, y_test.shape)
    linreg = linear_model.LinearRegression()
    linreg.fit(x_train,y_train)
    pred = linreg.predict(x_test)

    sns.regplot(y_test,pred, color='rebeccapurple')
    plt.show()
    
    x = x.rename(columns={'1': 'one', '3': 'three', '6':'six', '9':'nine', '12':'twelve'})
    x['sigma_rim'] = df.sigma_rim
    lsf = smf.ols(formula='sigma_rim ~ one + twelve + three + six + nine', data=x).fit()
    print(lsf.summary())
    params = lsf.params
    exit()
    # normal = scipy.stats.normaltest(sigma)
    # print(normal)
    # levene = scipy.stats.levene(sigma, acc_factors)
    # print(levene)
    # linres = scipy.stats.linregress(acc_factors, sigma)
    # print(linres)
    
    plt.plot(x, 'o', label='Original data', color='rebeccapurple')
    plt.plot(y, linres.intercept + linres.slope*y, label='Fitted line', color='orange')
    plt.legend()
    plt.xticks(acc_factors.unique())
    plt.ylabel('Sharpness in FWHM')
    plt.xlabel('Acceleration factor')
    plt.savefig(f'../../plots_saved/linreg.pdf')
    plt.show()


if __name__ == '__main__':
    main()