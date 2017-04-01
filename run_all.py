from Bio import SeqIO
import cgr as cg
import ssim
import numpy as np
import compute_mds as mds
import pandas as pd
from ggplot import *
import time
import os, glob


TEST_RUN = False
NUM_TEST_ITEMS = 20
VERTEBRATA_ONLY = True


def clear_tmpdir():
    filelist = glob.glob("tmp/*.*")
    for f in filelist:
        os.remove(f)


def run_sequences(data_file):
    seq_recs_med = SeqIO.parse(data_file, "gb")
    cgr = cg.CGR()
    idx = 0
    for seq in seq_recs_med:
        # break if a test run to only analyze 10 sequences
        if TEST_RUN and idx > NUM_TEST_ITEMS:
            break
        if len(seq.annotations['taxonomy']) < 5 or (VERTEBRATA_ONLY and seq.annotations['taxonomy'][4] != 'Vertebrata'):
            print "SKIPPING"
            print seq.annotations['taxonomy']
            continue
        idx += 1
        cgr.run(seq, idx)
    print "total sequnces: %d" % idx


def compute_ssim():
    DSSIM_matrix = ssim.run()
    np.save('DSSIM', DSSIM_matrix)
    print(DSSIM_matrix.shape)


def compute_mds():
    DSSIM_matrix = np.load('DSSIM.npy')
    coords, annotations = mds.run(DSSIM_matrix)
    np.save('coords', coords)
    np.save('annotations', annotations)


def report():
    coords = np.load('coords.npy')
    annotations = np.load('annotations.npy')
    df = pd.DataFrame()
    tax_level = 6  # i think that level 6 is reported in fig 2 of the paper
    df['Species'] = [a['taxonomy'][tax_level] for a in annotations]
    df['x'] = [x[0] for x in coords]
    df['y'] = [x[1] for x in coords]
    print(df.head())
    p = ggplot(df, aes('x', 'y', color='Species')) + geom_point() + theme_bw()
    #sns.set_context("notebook", font_scale=1.1)
    #sns.set_style("ticks")
    #sns.lmplot('x', 'y', data=df, hue='tax', scatter_kws={'marker': 'D', 's': 100}, fit_reg=False)
    #plt.title('Molecular Distance Map using Multidimensional Scaling Embedding')
    #plt.savefig('MDM_MDS_euc.png')
    p.save('MDM_MDS_euc.png')


T0 = 0
def tic():
    global T0
    T0 = time.time()


def toc():
    t1 = time.time()
    print "(%.5f seconds)" % (t1-T0)


def main():
    print 'Clearing tmpdir...'
    clear_tmpdir()
    print 'Running sequences...'
    tic()
    run_sequences('data/mitochondrion.1.genomic.gbff')
    toc()
    print 'Computing SSIM...'
    tic()
    compute_ssim()
    toc()
    print 'Computing MDS...'
    tic()
    compute_mds()
    toc()
    print 'Reporting...'
    tic()
    report()
    toc()


if __name__ == '__main__':
    main()
