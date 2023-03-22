import os
import glob
from astropy.io import fits
from hst_flystar.microlens import microlens_reduce
from hst_flystar import starlists
import pdb

# xym2mat/bar script for data in 
# /g/lu/data/microlens/hst/OB110462/v2022_11_11.
# It assumes everything's been reduced by hst1pass,
# all it is doing is the collation step. 

# Move to the top level, right above all the dates.

data_dir = './../00.DATA/' 
os.chdir(data_dir)

# Check whether single or multi exposure times
fits_files = glob.glob('*flt.fits')
exptimes = []
for ffile in fits_files:
    hdul = fits.open(ffile)
    exptime = hdul[0].header['EXPTIME']
    exptimes.append(exptime)

if len(exptimes) > 2:
    nn = len(exptimes) - 1
elif len(exptimes) == 2:
    nn = 2
elif len(exptimes) == 1:
    print('Only one frame, cannot reduce.')
    pdb.set_trace()
else:
    pdb.set_trace()

os.chdir('../')
    
# Other parameters.
out_format = 'uvwXYMpq'
work_dir = '01.XYM'
lowmag_lim = -6
minMagGood = -13.5
maxMagGood = -10
plotMagRange = [-22, 0]

microlens_reduce.microlens_reduce('F814W', nn, lowmag_lim, file_suffix=out_format,
                                  work_dir=work_dir, use_distortion=False)
os.chdir('./' + work_dir)
starlists.plot_matchup_err('MATCHUP_XYMEEE_final',
                           minMagGood=minMagGood, maxMagGood=maxMagGood,
                           plotMagRange=plotMagRange)
