import numpy as np
import astropy.io.fits as pyfits
import subprocess, sys

def calc_metric(file_name, mask_file):

    # block method
    trim = 200
    patch_size = 20

    # load files
    im = pyfits.open(file_name)[0].data
    mask = pyfits.open(mask_file)[0].data

    # trim the border
    im = im[trim:-1*trim,trim:-1*trim]
    mask = mask[trim:-1*trim,trim:-1*trim]

    # compute the patch means and RMSs
    means = []
    rms = []
    for i in range(int(im.shape[0]/patch_size)-1):
        for j in range(int(im.shape[1]/patch_size)-1):

            # cut the patches
            patch = im[i*patch_size:(i+1)*patch_size, j*patch_size:(j+1)*patch_size]
            pmask = mask[i*patch_size:(i+1)*patch_size, j*patch_size:(j+1)*patch_size]

            # test to see if there is >50% free pixels
            # mask is non-zero where there are sources
            if np.size(pmask[pmask==0]) >= 0.5*patch_size*patch_size:
                means.append(np.mean(patch[pmask==0]))
                rms.append(np.std(patch[pmask==0]))

    means = np.array(means)
    rms = np.array(rms)
    print(np.std(means), np.average(rms))


# main prog. --- add in mask
if __name__ == '__main__':

    # get filename arg. - sort this out a bit better.
    if len(sys.argv) < 2:
        print('Error: Missing input filenames.')
        print('Call as: python quick_metric.py <fits image> <object mask>')
        sys.exit()
    else:
        input_file = sys.argv[1]
        mask = sys.argv[2]

    calc_metric(input_file, mask)
