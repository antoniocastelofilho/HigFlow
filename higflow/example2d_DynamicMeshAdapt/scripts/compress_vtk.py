from glob import glob
from re import search as re_search
import os
from argparse import ArgumentParser
import vtk

def get_files_with_pattern(directory, pattern, regex):
    # find all vtk files that match the pattern
    vtk_files = sorted(glob(os.path.join(directory, pattern)))

    files_dict = {}

    for vtk_file in vtk_files:
        # extract the iteration number from the file name
        match = re_search(regex, vtk_file)
        if match:
            iteration = int(match.group(1))
            files_dict[iteration] = vtk_file

    return files_dict

def read_vtk(filename):
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.ReadAllTensorsOn()
    reader.Update()
    return reader.GetOutput()

def write_vtu(data, filename):
    writer = vtk.vtkXMLDataSetWriter()
    writer.SetFileName(filename)
    writer.SetInputData(data)
    writer.Write()


def read_vtk_files_write_vtu(files_dict):
    for (iteration, vtk_file) in files_dict.items():
        print(f"Processing {vtk_file}...")
        
        # Read VTK file
        data = read_vtk(vtk_file)
        
        # Generate output filename
        base_name = os.path.basename(vtk_file)
        name_without_ext = os.path.splitext(base_name)[0]
        output_file = os.path.join(directory, f"{name_without_ext}.vtu")
        
        # Write VTU file
        write_vtu(data, output_file)
        
        print(f"Converted to {output_file}")


# create an argument parser to receive the directory name as an argument
parser = ArgumentParser()
parser.add_argument("directory", help="directory containing vtk files")
args = parser.parse_args()

# get the directory name from the command-line argument
directory = args.directory

# specify the patterns for the vtk files
vtk_pattern = "*.print_*.vtk"
regex = r"print_(\d+)\.vtk"

files_dict = get_files_with_pattern(directory, vtk_pattern, regex)
read_vtk_files_write_vtu(files_dict)

vtk_pattern = "*.print_mult_*.vtk"
regex = r"print_mult_(\d+)\.vtk"

files_dict = get_files_with_pattern(directory, vtk_pattern, regex)
read_vtk_files_write_vtu(files_dict)

