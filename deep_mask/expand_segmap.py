import astropy.io.fits as pyfits
import numpy as np
import sys
import subprocess
from scipy import signal

class Seg_Map(object):
    
    def __init__(self, image_name, exp_size=8):

        self.exp_size = exp_size
        self.construct_mask(image_name)


    def expand_seg(self, mask):
        mask = signal.convolve2d(mask, np.ones((self.exp_size,self.exp_size)), mode='same')
        mask = np.floor(mask/self.exp_size**2)
        return mask

        
    def construct_mask(self, image_name):
        # read in segmap
        segmap = pyfits.open(image_name)
        mask = segmap[0].data
        self.head = segmap[0].header
        segmap.close()
        
        # change values such that source pixels are 0 and background are 1
        # is this in the right sense??
        mask[mask==1] = 2
        mask[mask==0] = 1
        mask[mask>1] = 0

        # expand sources in segmap by pre-defined amount
        self.mask = self.expand_seg(mask)



def run_sex(fname, weight_name, conf_file, param_name):
    # run sextractor
    command = 'sex -c {0} {1} -WEIGHT_IMAGE {2} -PARAMETER_FILE {3} -CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME seg.fits'.format(conf_file, fname, weight_name, param_name)
    subprocess.call(['/bin/bash', '-i', '-c', command]) # this form because I use an alias


def save_mask(mask, fname='mask.fits'):
    # save the mask
    hdu = pyfits.PrimaryHDU(mask.mask, header=mask.head)
    hdu.writeto(fname, overwrite=True)


# main prog.
if __name__ == '__main__':

    # get filename arg. - sort this out a bit better.
    if len(sys.argv) < 4:
        print('Error: No filename supplied.')
        print('Call as: python clean_im.py <image filename> <weight filename> <config filename>')
        sys.exit()
    else:
        fname = sys.argv[1]
        wname = sys.argv[2]
        conf_name = sys.argv[3]
        try:
            param_name = sys.argv[4]
        except:
            param_name = "../conf/default.param"

            
    run_sex(fname, wname, conf_name, param_name)
    mask = Seg_Map("seg.fits")
    save_mask(mask)
    
