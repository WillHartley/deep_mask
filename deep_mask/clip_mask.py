import numpy as np
import astropy.io.fits as pyfits
import sys, copy
from collections import namedtuple
from reproject import reproject_interp
from skimage.measure import EllipseModel
from scipy.spatial import ConvexHull

import matplotlib.pyplot as plt
plt.rc('font',**{'family':'serif','size':18})
plt.rc('text', usetex=True)


def parse_args(config, args):

    return config


# https://stackoverflow.com/questions/62980280/finding-neighboring-pixels-python
def candidate_neighbors(node):
    return ((node[0]-1, node[1]-1), (node[0]-1, node[1]), (node[0]-1, node[1]+1), (node[0], node[1]-1), 
            (node[0], node[1]+1), (node[0]+1, node[1]-1), (node[0]+1, node[1]), (node[0]+1, node[1]+1))

def neighboring_groups(nodes):
    remain = set(nodes)
    while len(remain) > 0:
        visit = [remain.pop()]
        group = []
        while len(visit) > 0:
            node = visit.pop()
            group.append(node)
            for nb in candidate_neighbors(node):
                if nb in remain:
                    remain.remove(nb)
                    visit.append(nb)
        yield tuple(group)


def save_mask(mask_data, name_base):
    # construct output name
    out_name = name_base.rsplit('.',1)[0]+"_clipmask.fits"
    # write out
    pyfits.writeto(out_name, mask_data, pyfits.Header(), overwrite=True)

    
def shift_mask(mask, new_head, orig_head):
    # use the astropy function to shift wcs and resample into the desired frame
    # build the mask hdu from the array and original header
    mask_hdu = pyfits.ImageHDU(mask, header=orig_head)
    new_mask, _ = reproject_interp(mask_hdu, new_head)
    return new_mask

def expand_pix(pix, exp, dims):
    # grow the mask region
    out_pix = np.flip(np.array(pix, dtype=int), 1)
    if exp==0:
        return out_pix
    else:
        # https://carstenschelp.github.io/2018/09/14/Plot_Confidence_Ellipse_001.html

        # get the convex hull
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.ConvexHull.html
        hull = out_pix[ConvexHull(out_pix).vertices,:]
        print(hull)
        
        # fit an ellipse to the hull
        # (https://stackoverflow.com/questions/39693869/fitting-an-ellipse-to-a-set-of-data-points-in-python)
        ell = EllipseModel()
        ell.estimate(hull)
        if ell.estimate(hull) is True:

            # is theta in rads or deg?
            xc, yc, a, b, theta = ell.params

            print("center = ",  (xc, yc))
            print("angle of rotation = ",  theta)
            print("axes = ", (a,b))

            # expand the axes by the desired amount.
            a *= (1+exp)
            b *= (1+exp)

            # define the set of enclosed pixels. - do I need a 1px offset here????
            ####### CHECK ############
            indx_arr = np.meshgrid(np.arange(dims[0]),np.arange(dims[1]))
            ix = indx_arr[0].flatten()
            iy = indx_arr[1].flatten()
            #xt = xc + a*cos(theta)*cos(t) - b*sin(theta)*sin(t)
            #yt = yc + a*sin(theta)*cos(t) + b*cos(theta)*sin(t)
            xell = (ix-xc)*np.cos(theta) + (iy-yc)*np.sin(theta)
            yell = -(ix-xc)*np.sin(theta) + (iy-yc)*np.cos(theta)
            mask_px = np.where(((xell/a)**2+(yell/b)**2 < 1.))[0]
            xy = np.array([(ix[i], iy[i]) for i in mask_px])
            
            return np.array(xy, dtype=int)

        else:
            print("Warning: ellipse estimate failed for group, ",pix)
            return out_pix



class clip_mask():

    def __init__(self, config):

        # I/O and set up
        self.config = config
        self.coadd_head = pyfits.open(config['coadd_file'])[0].header
        self.clip = np.genfromtxt(config['clip_file'])
        self.sefs = np.genfromtxt(config['se_list'], dtype=str)

        # remove clipped pixels that coincide with known bright stars
        if self.config['star_list'] is not None:
            self.delete_collisions()

        # gather clipped pixels into groups
        self.groups = self.build_groups()

        # create the coadd-level products
        self.coadd_mask = self.build_mask()
        save_mask(self.coadd_mask, self.config['coadd_file'])
        if self.config['output_regions'] is True:
            self.output_reg(self.groups, name_base=self.config['coadd_file'])
        if 'plot_name' in self.config.keys():
            self.save_plot(plot_name=self.config['plot_name'])

        # Single-epoch level products
        self.generate_masks()
        


    def generate_masks(self):
        # main function to produce the individual frame masks and translate them via wcs
        # get number of frames
        N_frames = len(self.sefs)

        # loop over sefs
        for iframe in range(N_frames):
            mask = self.build_mask(frame_num=iframe)
            sef_head = pyfits.open(self.sefs[iframe])[0].header
            mask = shift_mask(mask, sef_head, self.coadd_head)
            mask[mask>self.config['threshold']*self.config['mask_value']] = self.config['mask_value']
            mask[mask<self.config['mask_value']] = 0
            save_mask(mask, self.sefs[iframe])
            if self.config['output_regions'] is True:
                output_reg(self.groups, name_base=self.sefs[iframe], frame_num=iframe)

    def build_mask(self, frame_num=None):
        # make a blank array to insert mask bits into
        mask = np.zeros((self.coadd_head['NAXIS1'], self.coadd_head['NAXIS2'])).T
        indx_arr = np.meshgrid(np.arange(self.coadd_head['NAXIS1']),np.arange(self.coadd_head['NAXIS2']))
        # temporary: change all pixels in a group to 16
        for g in self.groups:
            # test to see whether we want to use this pixel group in the mask.
            if frame_num is None or frame_num==g.frame:
                #mask[(indx_arr[0]-g.x)**2+(indx_arr[1]-g.y)**2 < (g.size*np.pi)] = self.config['mask_value']
                pix_arr = expand_pix(g.pixels, g.expand, mask.shape)
                # the -1 is to convert from SExtractor unity indexing to python 0 indexing.
                mask[pix_arr[:,0]-1,pix_arr[:,1]-1] = self.config['mask_value']
        return mask
        
    def build_groups(self):
        # build a named tuple to contain group info
        group = namedtuple('group', ['x', 'y', 'size', 'frame', 'pixels', 'expand'])
        groups = []

        for frame in range(int(np.max(self.clip[:,0]))):
            for f in config['filters']:
                subset = np.where((self.clip[:,0]==frame)&(np.abs(self.clip[:,3])>f[0]))[0]
                pos = tuple([tuple(self.clip[s,1:3]) for s in subset])
                out = tuple(neighboring_groups(pos))
    
                # try a min number of pixels, 30
                for g in out:
                    if len(g) > f[1]:
                        groups.append(group(x=np.mean([g[i][0] for i in range(len(g))]),
                                            y=np.mean([g[i][1] for i in range(len(g))]),
                                            size=len(g),
                                            frame=frame,
                                            pixels=g,
                                            expand=f[2]))  
        return groups

    def delete_collisions(self):
        # remove clipped pixels that coincide with stars from being possibly masked.
        stars = np.genfromtxt(self.config['star_list'])
        # ........

        
    def output_reg(groups, name_base=None, frame_num=None):
        # need to add the header, make it obvious it is for the coadd - can make for SEFs too.
        # re-write to take named tuples.
        # try except the name, else use: 'clipped_mask.reg'
        # change to RA, Dec so that we don't need to translate coords.
        [print("circle({},{},{})".format(x_bar[i],y_bar[i],np.sqrt(group_size[i])*np.pi)) for i in range(len(x_bar))]
        # the regions were a bit small when dividing by pi

    
    def save_plot(self, plot_name='clipped_mask.png'):
        # update for the new group structure.
        fig = plt.figure(figsize=(9,6))
        ax = fig.add_subplot(111)
        ax.set_xlabel(r"\textbf{x}", fontsize=24)
        ax.set_ylabel(r'\textbf{y}', fontsize=24)
        ax.scatter(np.array(x_bar), np.array(y_bar),s=np.array(group_size))
        plt.tight_layout()
        plt.savefig(plot_name) 
        plt.show()


# main prog.
if __name__ == '__main__':

    # default config,
    # filters: a tuple of tuple containing pairs of (significance, min number of pix, expansion factor)
    # [here, expansion factor refers to....]
    # min number of pixels mustbe >6 if an exansion factor is used.
    # mask_ext: the fits extension of the mask, -1 will create new files, one for each input frame
    # mask_value: the value used in the output maskfiles - to be combined with existing values, bitwise.
    # star_list: a list of star positions - clipped pixels covering these positions will be deleted. List should contain, RA, Dec, radius for each object.
    # output_regions: whether to save a regions file of circles representing the size of each group.
    # threshold: the reprojection causes fractional values in the single-epoch masks. This value determines whether to keep a fraction or not (real: 0 to 1).
    config = {'filters': ((5, 400, 1.), (5, 30, 0.5), (15, 8, 0.2), (25, 1, 0.)),
              'mask_ext': -1,
              'mask_value': 1,
              'star_list': None,
              'output_regions': False,
              'threshold': 0.4}
        

    # get filename arg. - sort this out a bit better.
    if len(sys.argv) < 3:
        print('Error: Missing input filenames.')
        print('Call as: python clip_mask.py <clipped pixels filename> <coadded fits image> <list file of single-epoch images>')
        sys.exit()
    else:
        clip_file = sys.argv[1]
        coadd_file = sys.argv[2]
        se_list = sys.argv[3]

    config['clip_file'] = clip_file
    config['coadd_file'] = coadd_file
    config['se_list'] = se_list

    # parse remaining args
    if len(sys.argv) > 3:
        config = parse_args(config, sys.argv[3:])

    clip_mask(config)
