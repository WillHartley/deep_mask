#### WIP - doesn't run yet #####


import numpy as np
import astropy.io.fits as pyfits

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



class clip_mask():

    def __init__(self, config):

        self.config = config
        self.coadd_head = pyfits.open(config['coadd_file'])[0].header
        self.clip = np.genfromtxt(config['clip_file'])
        self.sefs = np.genfromtxt(config['se_list'], dtype=str)

        #self.delete_collisions()
        
        self.groups = self.build_groups()
        print(self.groups)

        
    def build_groups(self):
        x_bar = []
        y_bar = []
        group_size = []
        frame_id = []

        for frame in range(int(np.max(clip[:,0]))):
            subset = np.where((clip[:,0]==frame)&(np.abs(clip[:,3])>5))[0]
            pos = tuple([tuple(clip[s,1:3]) for s in subset])
            out = tuple(neighboring_groups(pos))
    
            # try a min number of pixels, 30
            for group in out:
                if len(group) > 29:
                    group_size.append(len(group))
                    x_bar.append(np.mean([group[i][0] for i in range(len(group))]))
                    y_bar.append(np.mean([group[i][1] for i in range(len(group))]))
                    frame_id.append(frame)
                    
        return [x_bar, y_bar, group_size, frame_id]

    def delete_collisions(self):
        # remove clipped pixels that coincide with stars from being possibly masked.
        

    def output_reg(self, reg_name='clipped_mask.reg'):
        # need to add the header, make it obvious it is for the coadd.
        [print("circle({},{},{})".format(x_bar[i],y_bar[i],np.sqrt(group_size[i])*np.pi)) for i in range(len(x_bar))]
        # the regions were a bit small when dividing by pi

    
    def save_plot(self, plot_name='clipped_mask.png'):
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
    # mask_ext: the fits extension of the mask, -1 will create new files, one for each input frame
    # mask_value: the value used in the output maskfiles - to be combined with existing values, bitwise.
    # star_list: a list of star positions - clipped pixels covering these positions will be deleted. List should contain, RA, Dec, radius for each object.
    config = {'filters': ((5, 30, 0.2)),
              'mask_ext': -1,
              'mask_value': 16,
              'star_list': None}
        

    # get filename arg. - sort this out a bit better.
    if len(sys.argv) < 3:
        print('Error: Missing input filenames.')
        print('Call as: python clip_mask.py <clipped pixels filename> <coadded fits image> <list file of single-epoch images>')
        sys.exit()
    else:
        clip_file = sys.argv[1]
        coadd_file = sys.argv[2]
        se_list = sys.argv[2]

    config['clip_file'] = clip_file
    config['coadd_file'] = coadd_file
    config['se_list'] = se_list

    # parse remaining args
    if len(sys.argv) > 3:
        config = parse_args(config, sys.argv[3:])

    clip_mask(config)
