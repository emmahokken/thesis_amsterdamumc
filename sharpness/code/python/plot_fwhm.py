import csv 
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt 

def plot_fwhm_per_subject(df, subj):
    '''
    Generates a scatter plot of the calculated FWHM (sharpness). 

    Args:
        df: DataFrame object containing FWHM information
        subj: subject to investigate
    '''
    
    fwhm = df.loc[(df.subj_id == subj) & (df.fields != 'ventl') &  (df.fields != 'ventr')]
    print(fwhm)
  
    plt.scatter(fwhm.acc_factor, fwhm.sigma_gt, c='rebeccapurple', marker='o')
    plt.scatter(fwhm.acc_factor, fwhm.sigma_rim, c='orange', marker='^')
    for i in range(len(fwhm.sigma_gt)):
        # plt.scatter(fwhm.acc_factor.iloc[i], fwhm.sigma_gt.iloc[i], c='rebeccapurple', marker='o')
        # plt.scatter(fwhm.acc_factor.iloc[i], fwhm.sigma_rim.iloc[i], c='orange', marker='^')
        plt.annotate(text=fwhm.fields.iloc[i],xy=(fwhm.acc_factor.iloc[i],fwhm.sigma_gt.iloc[i]))
        plt.annotate(text=fwhm.fields.iloc[i],xy=(fwhm.acc_factor.iloc[i],fwhm.sigma_rim.iloc[i]))

    plt.legend(['rim','gt'])
    plt.xticks(range(3,16,3))
    plt.xlabel('Acceleration factor')
    plt.ylabel('Sharpness in FWHM')
    plt.title('Sharpness as a function of acceleration factor for one subject')
    plt.savefig(f'../../plots_saved/{subj:03}_FWHM_all.pdf')
    plt.show()


def plot_per_region(df):
    '''
    Generates a line plot of the calculated FWHM (sharpness) per region for all axxeleration rates. 

    Args:
        df: DataFrame object containing FWHM information
    '''

    regions = df.fields.unique()
    acc_factors = df.acc_factor.unique()
    acc_factors = np.insert(acc_factors,0,0)
    cmap = plt.cm.jet(np.linspace(0,1,len(regions)))
    fig, ax = plt.subplots()
    for i, r in enumerate(regions):
        values = df.loc[df.fields == r]
        fwhm = []
        for a in acc_factors:
            if a == 0:
                acc = values.loc[values.acc_factor == 3]
                fwhm.append(acc.sigma_gt.mean())
                continue
            acc = values.loc[values.acc_factor == a]
            fwhm.append(acc.sigma_gt.mean())
        ax.plot(acc_factors,fwhm, color=cmap[i], label=r)
        ax.annotate(text=r,xy=(acc_factors[-1],fwhm[-1]),color=cmap[i],xytext=(5,0),textcoords='offset points',va='center')

    # plt.legend()
    plt.xticks(acc_factors)
    plt.xlabel('Acceleration rate')
    plt.ylabel('Sharpness in FWHM')
    plt.title('Sharpness of different ROIs over acceleration rates')
    plt.savefig(f'../../plots_saved/FWHM_per_region.pdf')
    plt.show()


if __name__ == '__main__':

    subj = 25
    begin = 3
    stop = 12

    plot_fwhm(subj, begin, stop)