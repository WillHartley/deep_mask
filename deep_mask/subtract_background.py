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


def invert_mask(mask):
    return -1*mask + 1
    
def combine_masks(obj_mask, pipe_mask):
    # take in the object mask, and the mask already present in the input file and combine
    # 1=masked, 0=ok.
    # The SE mask is the opposite (for ease of applying the filter), so need to convert
    pipe_mask = invert_mask(pipe_mask)
    obj_mask = obj_mask*pipe_mask
    obj_mask = invert_mask(obj_mask)
    return obj_mask
    
    

def subtract_image_backgrounds(deep_mask, image_list):
    # loop over each item in the list file
    for im_file in image_list:
        # read the target image header
        im = pyfits.open(im_file)
        
        # move the deep_mask into the right frame
        mask = remap_mask(deep_mask, im[0].header)

        # combine the new mask with the one from the single-epoch file
        mask = combine_masks(mask, im[1].data)

        # write to disk so SExtractor can access it
        save_mask(mask, im[0].header)

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

