#!/usr/bin/env python3
"""Fix VTK files affected by the CELL_TYPES local_str_size bug.

The bug in higflow_print_vtk3D_parallel_single used local_str_size=2
for CELL_TYPES, but VTK_CELL_TYPE=12 needs 3 bytes ("12 ").
This produced merged digits like "121212..." instead of "12 12 12...".
"""
import re, sys, os

def fix_vtk_file(filepath):
    with open(filepath, 'r') as f:
        content = f.read()

    pattern = r'(CELL_TYPES (\d+)\n)(\d+)'
    match = re.search(pattern, content)
    if not match:
        return False

    header = match.group(1)
    count = int(match.group(2))
    data = match.group(3)

    if len(data) == count * 2:
        fixed = ' '.join(data[i:i+2] for i in range(0, len(data), 2))
        content = content.replace(match.group(0), header + fixed)
        with open(filepath, 'w') as f:
            f.write(content)
        return True
    return False

def main():
    if len(sys.argv) < 2:
        print("Usage: fix_vtk_celltypes.py <file.vtk> [file2.vtk ...]")
        print("       fix_vtk_celltypes.py <directory/>")
        return

    paths = []
    for arg in sys.argv[1:]:
        if os.path.isdir(arg):
            for root, _, files in os.walk(arg):
                for f in files:
                    if f.endswith('.vtk'):
                        paths.append(os.path.join(root, f))
        else:
            paths.append(arg)

    fixed_count = 0
    for p in paths:
        if fix_vtk_file(p):
            print(f"Fixed: {p}")
            fixed_count += 1
        else:
            print(f"Skipped (no CELL_TYPES or already correct): {p}")

    print(f"\n{fixed_count} files fixed.")

if __name__ == "__main__":
    main()
