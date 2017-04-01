## SSIM related stuff

from skimage.measure import compare_ssim as ssim
from skimage import data, img_as_float, io
import numpy as np
import glob


def run():
    all_image_fnames = glob.glob('tmp/*.png')

    # sort to ensure that our annotations line up with our image comparisons
    sorted(all_image_fnames)

    n = len(all_image_fnames)

    DSSIM = np.zeros((n, n))
    for i in range(n):
        im1 = img_as_float(io.imread(all_image_fnames[i], as_grey=True))
        for j in range(n):
            im2 = img_as_float(io.imread(all_image_fnames[j], as_grey=True))
            DSSIM[i, j] = 1-ssim(im1, im2)
            #print "(%d,%d)" % (i,j)
    return DSSIM

if __name__ == '__main__':
    A = run()
    print(A)