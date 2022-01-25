#!/usr/bin/env python

import sys
import re

try:
    fname = sys.argv[1]
except IndexError:
    print "Log file not specified!"
    sys.exit(1)

fd = open(fname, 'r')

p = re.compile(r'^==(\d+)== .*$')

outs = {}

for line in fd:
    m = p.match(line)
    if m:
        pid = m.group(1)
        try:
            out = outs[pid]
        except KeyError:
            out_fname = '.'.join([fname, pid])
            out = open(out_fname, 'w')
            print 'Writing file "{}".'.format(out_fname)
            outs[pid] = out
        out.write(line)
