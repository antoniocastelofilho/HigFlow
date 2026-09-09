#!/usr/bin/env python3
import ruamel.yaml
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

yaml_path = 'input/load.par.contr.yaml'
with open(yaml_path, 'r') as f:
    text = f.read()

for p in ['Re', 'Ca', 'Fr', 'dt', 'numsteps', 'dts', 'dtp']:
    val = get_arg(p + '=')
    if val is not None:
        text = re.sub(r'(?<=^' + p + r': )\S+', val, text, flags=re.MULTILINE)

grav = get_arg('gravity=')
if grav is not None:
    text = re.sub(r'(?<=^  gravity: )\w+', grav, text, flags=re.MULTILINE)

with open(yaml_path, 'w') as f:
    f.write(text)

Re = get_arg('Re=', '50000')
Ca = get_arg('Ca=', '0.045')
dt = get_arg('dt=', '0.0001')
numsteps = get_arg('numsteps=', '80000')
dts = get_arg('dts=', '0.05')
dtp = get_arg('dtp=', '0.005')
out = f'DamBreak-Re{Re}-Ca{Ca}'
print(out, end='')
