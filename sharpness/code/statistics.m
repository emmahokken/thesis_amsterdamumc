% Load in files 
df = readtable('../results/FWHM.csv')
ran = ranova(df, 'WithinModel', 'acc_factor')