#!/usr/bin/env python3
import ruamel.yaml
import re, os, sys

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

# Read and modify par.contr.yaml
yaml_path = 'input/load.par.contr.yaml'
with open(yaml_path, 'r') as f:
    text = f.read()

# Parse parameters
for p in ['Re', 'Ca', 'Fr', 'dt', 'numsteps', 'dts', 'dtp']:
    val = get_arg(p + '=')
    if val is not None:
        text = re.sub(r'(?<=^' + p + r': )\S+', val, text, flags=re.MULTILINE)

# Gravity
grav = get_arg('gravity=')
if grav is not None:
    text = re.sub(r'(?<=^  gravity: )\w+', grav, text, flags=re.MULTILINE)

with open(yaml_path, 'w') as f:
    f.write(text)

# Update ADAPT_ENABLED / ADAPT_FREQ in header
adapt = get_arg('adapt=')
if adapt is not None:
    enabled = 1 if adapt.lower() in ('true', '1') else 0
    with open('ns-example-2d.h', 'r') as f:
        htext = f.read()
    htext = re.sub(r'(#define\s+ADAPT_ENABLED\s+)[01]', f'\\g<1>{enabled}', htext)
    with open('ns-example-2d.h', 'w') as f:
        f.write(htext)

freq = get_arg('adapt_freq=')
if freq is not None:
    with open('ns-example-2d.h', 'r') as f:
        htext = f.read()
    htext = re.sub(r'(#define\s+ADAPT_FREQ\s+)\d+', f'\\g<1>{freq}', htext)
    with open('ns-example-2d.h', 'w') as f:
        f.write(htext)

# Output name
Re = get_arg('Re=', '3130000')
Ca = get_arg('Ca=', '0.045')
dt = get_arg('dt=', '0.0005')
numsteps = get_arg('numsteps=', '16000')
dts = get_arg('dts=', '0.05')
dtp = get_arg('dtp=', '0.005')
out = f'DamBreak-Re{Re}-Ca{Ca}-AMR'
print(out, end='')
