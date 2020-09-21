#!/usr/bin/env python

import struct
from sys import stdout as out

def unit_reader(fname, fmt):
    size = struct.calcsize(fmt)
    with open(fname, 'rb') as fd:
        while True:
            b = fd.read(size)
            if len(b) < size:
                return
            yield struct.unpack(fmt, b)[0]

rowr = unit_reader('matrix-row-refs.dat', 'I')
cols = unit_reader('matrix-col-ids.dat', 'I')
elems = unit_reader('matrix-elems.dat', 'd')
rhs = unit_reader('matrix-rhs.dat', 'd')

ini = next(rowr)
for i, end in enumerate(rowr):
    out.write('{}:'.format(i))
    while ini != end:
        out.write(' {} ({:e}),'.format(next(cols), next(elems)))
        ini += 1
    out.write(' | {:e}\n'.format(next(rhs)))
