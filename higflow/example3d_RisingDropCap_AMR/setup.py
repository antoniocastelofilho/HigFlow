#!/usr/bin/env python3
import re, sys

def pos_in_args(key):
    for i, a in enumerate(sys.argv[1:]):
        if a.startswith(key):
            return i
    return -1

def get_arg(key, default=None):
    p = pos_in_args(key)
    if p != -1:
        return sys.argv[p+1].split('=')[1]
    return default

yaml_path = 'input/example-3d.load.par.contr.yaml'
with open(yaml_path, 'r') as f:
    text = f.read()

for p in ['Re', 'Ca', 'Fr', 'dt', 'numsteps', 'dts', 'dtp']:
    val = get_arg(p + '=')
    if val is not None:
        text = re.sub(r'(^\s*' + p + r': )\S+', r'\g<1>' + val, text, flags=re.MULTILINE)

with open(yaml_path, 'w') as f:
    f.write(text)

Re = get_arg('Re=', '35')
Ca = get_arg('Ca=', '3.57')
out = f'RisingDrop3D-Re{Re}-Ca{Ca}'
print(out, end='')
