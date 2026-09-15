#!/usr/bin/env python3

from pathlib import Path
from sys import argv
from enum import Enum
import math

width_factor_dict = {
    "square": 1,
    "short" : 2,
    "medium" : 4,
    "long"  : 8,
    "vlong"  : 16
}
ini_Ny_dict = {
    "std" : 40,
    "coarse" : 20,
    "fine": 80,
    "vcoarse": 10,
    "vfine": 160
}
posy_type_dict = {
    "sym" : [-1.0, 1.0],
    "br"  : [0.0, 1.0]
}
# enum of special names
class SpecialNames(Enum):
    NONE = 0
    DROP = 1
# special names variable
spec_name = SpecialNames.NONE

# get program args
args = argv

if(len(args)==1):
    print("Usage: python3 gen_mesh_channel.py num_ref width_name(opt) ini_Ny_name(opt) posy_type(opt) dx/dy(opt) mesh_name(opt)")
    print("num_ref: number of refinements")
    print("width_name: width factor of the mesh: square (1x1), short (1x2), medium (1x4), long (1x8), vlong (1x16)")
    print("ini_Ny_name: initial Ny: std (40), coarse (20), fine (80), vcoarse (10), vfine (160) or custom number (if given a positive integer)")
    print("posy_type: type of posy: sym ([-1, 1]), br ([0, 1])")
    print("dx/dy: ratio of dx to dy (default is 1)")
    print("mesh_name: extra name of the mesh at the beginning - special extra names include:")
    print("---> 'drop<radius>_<num_ref_drop>' for refinement in the middle of the channel for a drop of radius <radius>")
    print("-----> and <num_ref_drop> middle refinements (default is <num_ref>)")
    exit(1)

num_ref = int(args[1])
if(num_ref < 1):
    print("num_ref must be a positive integer")
    exit(1)
    
if(len(args)>=3):
    width_name = args[2]
    if(width_name not in width_factor_dict): # explain the current format
        print("width_name must be one of the following:")
        print(width_factor_dict.keys())
        exit(1)
else:
    width_name = "medium"
width_factor = width_factor_dict[width_name]


if(len(args)>=4):
    ini_Ny_name = args[3]
    # check if ini_Ny_name is a positive integer
    if(ini_Ny_name.isdigit()):
        ini_Ny = int(ini_Ny_name)
        if(ini_Ny < 1):
            print("ini_Ny must be a positive integer")
            exit(1)
    else:
        if(ini_Ny_name not in ini_Ny_dict): # explain the current format
            print("ini_Ny_name must be one of the following:")
            print(ini_Ny_dict.keys())
            exit(1)
        ini_Ny = ini_Ny_dict[ini_Ny_name]
        #ini_Ny_name = str(ini_Ny)
else:
    ini_Ny_name = "std"
    ini_Ny = ini_Ny_dict[ini_Ny_name]

if(len(args)>=5):
    posy_type = args[4]
    if(posy_type not in posy_type_dict): # explain the current format
        print("posy_type must be one of the following:")
        print(posy_type_dict.keys())
        exit(1)
else:
    posy_type = "sym"
posy = posy_type_dict[posy_type]

if(len(args)>=6):
    try:
        h2v_ratio = float(args[5])
    except ValueError:
        print("dx/dy must be a number")
        exit(1)
else:
    h2v_ratio = 1

mesh_name = ""
if(len(args) == 7):
    mesh_name = args[6]
    if('drop' in mesh_name):
        if(len(mesh_name.split('_'))>2):
            print("invalid mesh name")
            exit(1)
        if('_' not in mesh_name):
            num_ref_drop = num_ref
            if(num_ref_drop > num_ref):
                print("num_ref_drop must be less than or equal to num_ref")
                exit(1)
            mesh_name0 = mesh_name
        else:
            mesh_name0 = mesh_name.split('_')[0]
            mesh_name1 = mesh_name.split('_')[1]
            # check if string is postive integer
            if(not mesh_name1.isdigit()):
                print("invalid mesh name")
                exit(1)
            num_ref_drop = int(mesh_name1)
            if(num_ref_drop < 1):
                print("num_ref must be a positive integer")
                exit(1)

        if(mesh_name0[:4] != "drop"):
            print("invalid mesh name")
            exit(1)
        # check if string is valid float
        try:
            radius = float(mesh_name0[4:])
        except ValueError:
            print("radius must be a number")
            exit(1)
        if(radius <= 0):
            print("radius must be a positive number")
            exit(1)
        elif(radius >= 1):
            print("radius must be less than 1")
            exit(1)
        spec_name = SpecialNames.DROP
        drop_radius = radius


        

name = mesh_name
if(num_ref>1):
    name = name + f"_ref{num_ref}"
if(width_name != "medium"):
    name = name + f"_{width_name}"
if(ini_Ny_name != "std"):
    name = name + f"_{ini_Ny_name}"
if(posy_type != "sym"):
    name = name + f"_{posy_type}"
if(h2v_ratio != 1):
    name = name + f"_hv{h2v_ratio}"

if(name == ""):
    name = "std"

# trim underline if first character
if(name[0] == "_"):
    name = name[1:]

ref_Ny = [16]*(num_ref-1)
ref_Ny.append(16)
ref_Nx = ref_Ny

# create folder with mesh name
directory = Path(name)
directory.mkdir(parents=True, exist_ok=True)
directory = Path(name + "/domain")
directory.mkdir(parents=True, exist_ok=True)
directory = Path(name + "/bc")
directory.mkdir(parents=True, exist_ok=True)

# channel positions and dimensions
Ly = posy[1] - posy[0]
widthx = width_factor * Ly
posx = [0.0, 0.0 + widthx]
Lx = posx[1] - posx[0]

ini_Nx = int(ini_Ny * Lx / Ly / h2v_ratio)
ini_dx = Lx / ini_Nx
ini_dy = Ly / ini_Ny

# create file for domain
filename_domain = name + "/domain/" + name + "-d.amr"
file_path = Path(filename_domain)
file_path.touch()

with open(filename_domain, "w") as f:
    f.write(f"{posx[0]} {posx[1]} {posy[0]} {posy[1]}\n")

    num_levels = num_ref
    if(spec_name == SpecialNames.DROP and num_ref_drop < num_ref):
        num_levels = num_levels + 1

    f.write(f"{num_levels}\n")
    f.write(f"{ini_dx} {ini_dy} 1\n")
    f.write(f"1 1 {ini_Nx} {ini_Ny}\n")
    curr_dx, curr_dy = ini_dx, ini_dy
    curr_Nx, curr_Ny = ini_Nx, ini_Ny
    for i in range(1,num_ref):
        curr_dx, curr_dy = curr_dx/2, curr_dy/2
        curr_Nx, curr_Ny = curr_Nx*2, curr_Ny*2

        num_patches = 4
        if(spec_name == SpecialNames.DROP and i < num_ref_drop):
            num_patches = num_patches + 1

        f.write(f"{curr_dx} {curr_dy} {num_patches}\n")
        # horizontal patches
        f.write(f"1 1 {curr_Nx} {ref_Ny[i]}\n")
        f.write(f"1 {1+curr_Ny-ref_Ny[i]} {curr_Nx} {ref_Ny[i]}\n")
        # vertical patches
        f.write(f"1 {ref_Ny[i]+1} {ref_Nx[i]} {curr_Ny - 2*ref_Ny[i]}\n")
        f.write(f"{1+curr_Nx-ref_Nx[i]} {ref_Ny[i]+1} {ref_Nx[i]} {curr_Ny - 2*ref_Ny[i]}\n")
        if(spec_name == SpecialNames.DROP): # middle patches
            if(i < num_ref_drop):
                mid_Ny = (curr_Ny+1)/2
                if(posy_type == "sym"):
                    rad_inc = int(drop_radius * curr_Ny/2)  
                elif(posy_type == "br"):
                    rad_inc = int(drop_radius * curr_Ny)
                rad_inc = max(rad_inc,1) + 4
                upper_Ny = math.ceil(mid_Ny) + rad_inc
                lower_Ny = math.floor(mid_Ny) - rad_inc
                for j in range(i+1,num_ref_drop):
                    sub_ref_inc = max(int(ref_Ny[-1]/2**(j-i)),1)
                    upper_Ny = upper_Ny + sub_ref_inc
                    lower_Ny = lower_Ny - sub_ref_inc
                f.write(f"{ref_Nx[i]+1} {lower_Ny} {curr_Nx - 2*ref_Nx[i]} {upper_Ny - lower_Ny + 1}\n")

    if(spec_name == SpecialNames.DROP and num_ref_drop < num_ref): # multiphase mesh thickness
        div_factor = 2**(num_ref_drop-1)
        f.write(f"{ini_dx/div_factor} {ini_dy/div_factor} 0\n") 


# create file for boundary conditions
filename_bc0 = name + "/bc/" + name + "-bc-0.amr"
file_path = Path(filename_bc0)
file_path.touch()

with open(filename_bc0, "w") as f:
    f.write(f"{posx[0]} {posx[0]} {posy[0]} {posy[1]}\n")
    f.write(f"{num_ref}\n")
    curr_dy, curr_Ny = ini_dy, ini_Ny
    for i in range(num_ref):
        f.write(f"0.0 {curr_dy} 1\n")
        f.write(f"1 1 1 {curr_Ny}\n")
        curr_dy, curr_Ny = curr_dy/2, curr_Ny*2

filename_bc2 = name + "/bc/" + name + "-bc-2.amr"
file_path = Path(filename_bc2)
file_path.touch()

with open(filename_bc2, "w") as f:
    f.write(f"{posx[1]} {posx[1]} {posy[0]} {posy[1]}\n")
    f.write(f"{num_ref}\n")
    curr_dy, curr_Ny = ini_dy, ini_Ny
    for i in range(num_ref):
        f.write(f"0.0 {curr_dy} 1\n")
        f.write(f"1 1 1 {curr_Ny}\n")
        curr_dy, curr_Ny = curr_dy/2, curr_Ny*2

filename_bc1 = name + "/bc/" + name + "-bc-1.amr"
file_path = Path(filename_bc1)
file_path.touch()

with open(filename_bc1, "w") as f:
    f.write(f"{posx[0]} {posx[1]} {posy[1]} {posy[1]}\n")
    f.write(f"{num_ref}\n")
    curr_dx, curr_Nx = ini_dx, ini_Nx
    for i in range(num_ref):
        f.write(f"{curr_dx} 0.0 1\n")
        f.write(f"1 1 {curr_Nx} 1\n")
        curr_dx, curr_Nx = curr_dx/2, curr_Nx*2

filename_bc3 = name + "/bc/" + name + "-bc-3.amr"
file_path = Path(filename_bc3)
file_path.touch()

with open(filename_bc3, "w") as f:
    f.write(f"{posx[0]} {posx[1]} {posy[0]} {posy[0]}\n")
    f.write(f"{num_ref}\n")
    curr_dx, curr_Nx = ini_dx, ini_Nx
    for i in range(num_ref):
        f.write(f"{curr_dx} 0.0 1\n")
        f.write(f"1 1 {curr_Nx} 1\n")
        curr_dx, curr_Nx = curr_dx/2, curr_Nx*2