import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd 

def plot_fwhm_per_subject(df, subj):
    '''
    Generates a scatter plot of the calculated FWHM (sharpness). 

    Args:
        df: DataFrame object containing FWHM information
        subj: subject to investigate
    '''
    
    fwhm = df.loc[(df.subj_id == subj) & (df.fields != 'ventl') &  (df.fields != 'ventr')]
  
    plt.scatter(fwhm.acc_factor, fwhm.sigma_gt, c='rebeccapurple', marker='o', label='GT')
    plt.scatter(fwhm.acc_factor, fwhm.sigma_rim, c='orange', marker='^', label='RIM')
    for i in range(len(fwhm.sigma_gt)):
        # plt.scatter(fwhm.acc_factor.iloc[i], fwhm.sigma_gt.iloc[i], c='rebeccapurple', marker='o')
        # plt.scatter(fwhm.acc_factor.iloc[i], fwhm.sigma_rim.iloc[i], c='orange', marker='^')
        plt.annotate(text=fwhm.fields.iloc[i],xy=(fwhm.acc_factor.iloc[i],fwhm.sigma_gt.iloc[i]))
        plt.annotate(text=fwhm.fields.iloc[i],xy=(fwhm.acc_factor.iloc[i],fwhm.sigma_rim.iloc[i]))

    plt.legend()
    plt.xticks(range(3,16,3))
    plt.xlabel('Acceleration factor')
    plt.ylabel('Sharpness in FWHM')
    plt.title('Sharpness as a function of acceleration factor for one subject')
    plt.savefig(f'../../plots_saved/{subj:03}_FWHM_all.pdf')
    plt.show()


def plot_per_region(df, version):
    '''
    Generates a line plot of the calculated FWHM (sharpness) per region for all axxeleration rates. 

    Args:
        df: DataFrame object containing FWHM information
        version: whether to only plot the ground truth or the reconstructed results
    '''

    regions = df.fields.unique()
    acc_factors = df.acc_factor.unique()
    acc_factors = np.insert(acc_factors,0,0)
    cmap = plt.cm.jet(np.linspace(0,1,len(regions)))
    
    for i, r in enumerate(regions):
        values = df.loc[df.fields == r]
        fwhm = []
        for a in acc_factors:
            if a == 0:
                acc = values.loc[values.acc_factor == 3]
                if version == 'gt':
                    fwhm.append(acc.sigma_gt.mean())
                else:
                    fwhm.append(acc.sigma_rim.mean())
                continue
            acc = values.loc[values.acc_factor == a]
            if version == 'gt':
                fwhm.append(acc.sigma_gt.mean())
            else:
                fwhm.append(acc.sigma_rim.mean())
        plt.plot(acc_factors,fwhm, color=cmap[i], label=r)
        # plt.annotate(text=r,xy=(acc_factors[-1],fwhm[-1]),color=cmap[i],xytext=(5,0),textcoords='offset points',va='center')

    plt.legend(bbox_to_anchor=(1.05,1),loc='upper left')
    plt.xticks(acc_factors)
    plt.xlabel('Acceleration rate')
    plt.ylabel('Sharpness in FWHM')
    plt.title(f'Sharpness of different ROIs over acceleration rates for the {version.upper()}')
    plt.tight_layout()
    plt.savefig(f'../../plots_saved/FWHM_per_region_{version}.pdf')
    plt.show()


def boxplot(df):
    '''
    Generate a boxplot for the FWHM data

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

    df.boxplot(column=['sigma_rim'],by='acc_factor', grid=True, boxprops=dict(color='orange'),whiskerprops=dict(color='orange'),medianprops=dict(color='rebeccapurple'))
    plt.title('Distribution of sharpness in FWHM per acceleration factor.')
    plt.suptitle('')
    plt.xlabel('Acceleration factor')
    plt.ylabel('Sharpness in FWHM')
    plt.savefig(f'../../plots_saved/FWHM_boxplot.pdf')
    plt.show()
