import numpy as np
from reproject import reproject_interp
import astropy.io.fits as pyfits
import subprocess


def save_mask(mask, header):
    pyfits.writeto('mask.fits', mask, header, overwrite=True)

    
def remap_mask(mask_hdu, target_header):
    # use the astropy function to shift wcs and resample into the desired frame
    new_mask, _ = reproject_interp(mask_hdu, target_header)
    return new_mask


def run_sextr(fname):
    # construct the output file name
    out_name = fname.rsplit(".", 1)[0]+"_bksub.fits"

    # SExtractor command #-CHECKIMAGE_TYPE -BACKGROUND
    # Find a better way to do the config path
    sextr_cmd = "sex -c ../conf/backsub.conf {0} -CHECKIMAGE_NAME {1}".format(fname, out_name)
    subprocess.call(sextr_cmd, shell=True)


def subtract_image_backgrounds(deep_mask, image_list):
    # loop over each item in the list file
    for im in image_list:
        # read the target image header
        im_head = pyfits.open(im)[0].header
        
        # move the deep_mask into the right frame
        mask = remap_mask(deep_mask, im_head)

        # write to disk so SExtractor can access it
        save_mask(mask, im_head)

        # run SExtractor to get background subtracted image
        run_sextr(fname)



# main prog.
if __name__ == '__main__':

    # get filename arg. - sort this out a bit better.
    if len(sys.argv) < 3:
        print('Error: Missing input filenames.')
        print('Call as: python subtract_background.py <mask filename> <list file of images>')
        sys.exit()
    else:
        mask_name = sys.argv[1]
        list_file = sys.argv[2]

    deep_mask = pyfits.open(mask_name)[0]
    image_list = np.genfromtxt(list_file)

    subtract_image_backgrounds(deep_mask, image_list)

