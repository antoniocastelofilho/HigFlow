#!/usr/bin/env python3
"""Apply IN=... overrides to the YAML and print the output-name stem."""
import re
import sys


def pos_in_args(key):
    """Return the index (in argv[1:]) of the first token starting with key."""
    for i, a in enumerate(sys.argv[1:]):
        if a.startswith(key):
            return i
    return -1


def get_arg(key, default=None):
    """Return the value of a key=value token, or default if absent."""
    p = pos_in_args(key)
    if p != -1:
        return sys.argv[p + 1].split('=')[1]
    return default


def main():
    """Rewrite the control YAML from CLI overrides and emit the name stem."""
    yaml_path = 'input/load.par.contr.yaml'
    with open(yaml_path, 'r', encoding='utf-8') as f:
        text = f.read()

    for p in ['Re', 'Ca', 'Fr', 'dt', 'numsteps', 'dts', 'dtp']:
        val = get_arg(p + '=')
        if val is not None:
            # Keys are indented under their YAML section, so allow (and
            # preserve) leading whitespace before the key name.
            pattern = r'(?m)^(\s*' + re.escape(p) + r': )\S+'
            text = re.sub(pattern, lambda m, v=val: m.group(1) + v, text)

    grav = get_arg('gravity=')
    if grav is not None:
        text = re.sub(r'(?m)^(\s*gravity: )\w+',
                      lambda m, v=grav: m.group(1) + v, text)

    with open(yaml_path, 'w', encoding='utf-8') as f:
        f.write(text)

    re_num = get_arg('Re=', '35')
    ca_num = get_arg('Ca=', '3.57')
    print(f'RisingDrop2D-Re{re_num}-Ca{ca_num}', end='')


if __name__ == '__main__':
    main()
