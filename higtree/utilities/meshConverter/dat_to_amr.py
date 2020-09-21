#!/usr/bin/env python

import sys
import struct
import numpy as np

def fortran_reader(fd):
    while True:
        header = fd.read(4)
        if not header:
            break
        header, = struct.unpack('I', header)
        yield fd.read(header+4)[:header]

def calc_size(level):
    # Calculate the size of the domain
    d = np.array(level[:3])
    n = np.array(level[-1][0][3:6])
    return n * d

def main():
    if len(sys.argv) < 3:
        print 'Usage:\n {} [--normalize] in_file.dat out_file.amr\n'.format(sys.argv[0])
        return

    if sys.argv[1] == '--normalize':
        normalize = True
        in_filename = sys.argv[2]
        out_filename = sys.argv[3]
    else:
        normalize = False
        in_filename = sys.argv[1]
        out_filename = sys.argv[2]

    with open(in_filename, 'rb') as dat:
        dat = fortran_reader(dat)
        levels = []
        try:
            while True:
                buf = next(dat) # expecting dx, dy, dz, n_patches
                l_info = struct.unpack('dddi', buf)
                n_patches = l_info[-1]
                patches = []
                for i in xrange(n_patches):
                    # expecting 13 ints: ix, iy, iz, nx, ny, nz, plus 7 ints
                    buf = next(dat)
                    patches.append(struct.unpack('i'*13, buf)[:6])
                levels.append(list(l_info[:3]) + [patches])
        except StopIteration:
            pass

    # Calculate the size of the domain
    dsize = calc_size(levels[-1])

    # Find first level whose patch size is different from the domain size
    for i, level in enumerate(reversed(levels)):
        if (dsize != calc_size(level)).any() or len(level[-1]) != 1:
            break

    # Filter out virtual levels
    i = 1 - i
    if i:
        levels = levels[:i]

    if normalize:
        scale = 1/dsize
    else:
        scale = np.ones(3)

    with open(out_filename, 'w') as amr:
        # Write header
        amr.write('0 {} 0 {} 0 {}\n'.format(*(dsize * scale)))
        amr.write('{}\n'.format(len(levels)))

        # Write level
        for level in reversed(levels):
            amr.write("{} {} {}\n".format(*(level[:3] * scale)))
            amr.write("{}\n".format(len(level[-1])))
            # write each patch
            for patch in level[-1]:
                amr.write(("{}" + " {}" * 5 + "\n").format(*patch))

if __name__ == '__main__':
    main()
