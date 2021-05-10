import os 
from shutil import copyfile

dirs = [d[1] for d in os.walk('../../data/recon')][0]

print(dirs)

for d in dirs:

    acc = [a for a in d.split('_') if a.isdigit()][0]
    print(acc)
    gt = [f[2] for f in os.walk(f'../../data/recon/{d}/R2star_map_gt')][0]
    rim = [f[2] for f in os.walk(f'../../data/recon/{d}/R2star_map_gt')][0]
    # print(rim)
    
    # for r in rim:
    #     subj = [z for z in r.split('_') if z.isdigit()][0][1:]
    #     copyfile(f'../../data/recon/{d}/R2star_map_rim/{r}',  f'data/{subj}_r2star_rim_{acc}.nii')
        
    for g in gt:
        subj = [z for z in g.split('_') if z.isdigit()][0][1:]
        copyfile(f'../../data/recon/{d}/R2star_map_gt/{g}',  f'data/{subj}_r2star_gt.nii')

    exit()

    # '005_r2star_rim'
    # '005+r2star_gt'