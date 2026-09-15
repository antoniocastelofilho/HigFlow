import vtk
from glob import glob
from re import search as re_search
import os
from argparse import ArgumentParser

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

# Function to read a single timestep
def read_timestep(xdmf_file):
    reader = vtk.vtkXdmfReader()
    reader.SetFileName(xdmf_file)
    reader.Update()
    return reader.GetOutputDataObject(0)

# Function to merge datasets from multiple processes
def merge_datasets(datasets: vtk.vtkMultiBlockDataSet):
    appender = vtk.vtkAppendDataSets()
    numblocks = datasets.GetNumberOfBlocks()
    for i in range(numblocks):
        dataset = datasets.GetBlock(i)
        appender.AddInputData(dataset)
    appender.Update()
    return appender.GetOutput()


# create an argument parser to receive the directory name as an argument
parser = ArgumentParser()
parser.add_argument("directory", help="directory containing vtk files")
args = parser.parse_args()

# get the directory name from the command-line argument
directory = args.directory

# specify the patterns for the vtk files
vtk_pattern = "*.print.t*.xdmf"
regex = r"print.t(\d+)\.xdmf"

xdmf_files = get_files_with_pattern(directory, vtk_pattern, regex)

# Read and merge data for each timestep
merged_timesteps = []
for step in xdmf_files:
    dataset = read_timestep(xdmf_files[step])
    merged_dataset = merge_datasets(dataset)  # If the XDMF already combines multiple processes, this step might be unnecessary
    merged_timesteps.append(merged_dataset)

# Write the merged dataset
for i, dataset in enumerate(merged_timesteps):
    #print with padded zeros
    filename = os.path.join(directory, f"ns.print.t{str(i).zfill(6)}.vtu")
    writer = vtk.vtkXMLDataSetWriter()
    writer.SetFileName(filename)
    writer.SetInputData(dataset)
    writer.Write()

# remove the xdmf files
for step in xdmf_files:
    os.remove(xdmf_files[step])
    print(f"Removed {xdmf_files[step]}")
# remove the h5 files
h5_files = glob(os.path.join(directory, "*.h5"))
for h5_file in h5_files:
    os.remove(h5_file)
    print(f"Removed {h5_file}")
