# DOES NOT WORK WHEN PRINTING MULTIPLE VECTOR OR TENSOR FIELDS

from glob import glob
from re import search as re_search
from os.path import join as os_path_join
from os import remove as os_remove
from argparse import ArgumentParser
from vtk import vtkMultiPieceDataSet, vtkMultiBlockDataSet, vtkUnstructuredGridReader, vtkCompositeDataSet
from vtk import vtkXdmfWriter, vtkXMLMultiBlockDataWriter, vtkXMLDataSetWriter

# create an argument parser to receive the directory name as an argument
parser = ArgumentParser()
parser.add_argument("directory", help="directory containing vtk files")
args = parser.parse_args()

# get the directory name from the command-line argument
directory = args.directory

# specify the pattern for the vtk files
vtk_pattern = "*.print_*-*.vtk"

# find all vtk files that match the pattern
vtk_files = glob(os_path_join(directory, vtk_pattern))

# create a dictionary to hold the grouped vtk files
files_dict = {}

# iterate over the vtk files
for vtk_file in vtk_files:
    # extract the iteration and process number from the file name
    match = re_search(r"print_(\d+)-(\d+)\.vtk", vtk_file)
    if match:
        process = int(match.group(1))
        iteration = int(match.group(2))
        # add the vtk file to the corresponding group
        if iteration not in files_dict:
            files_dict[iteration] = {}
        files_dict[iteration][process] = vtk_file
        beg_name = vtk_file.split(".print_", 1)[0]

print(f"merging files at {beg_name}")

writer = vtkXMLDataSetWriter()
#writer = vtkXMLMultiBlockDataWriter()
#writer.SetNumberOfTimeSteps(len(files_dict))
# iterate over the vtk dictionary
for iteration in files_dict:
    # create a new vtk file with the desired naming convention
    new_filename = beg_name + f".print_grouped_{iteration}.xmf"
    # multiblock is deprecated
    #multipiece = vtkMultiPieceDataSet()
    multiblock = vtkMultiBlockDataSet()
    for process, filename in files_dict[iteration].items():
        reader = vtkUnstructuredGridReader()
        reader.ReadAllFieldsOn()
        reader.ReadAllVectorsOn()
        reader.ReadAllTensorsOn()
        reader.ReadAllScalarsOn()
        reader.SetFileName(filename)
        reader.Update()

        # Get the point data
        point_data = reader.GetOutput().GetPointData()
        cell_data = reader.GetOutput().GetCellData()

        # Set Tensor (has to be unique)
        for i in range(point_data.GetNumberOfArrays()):
            array = point_data.GetArray(i)
            if array.GetNumberOfComponents() == 9:
                point_data.SetTensors(array)

        #multipiece.SetPiece(process, reader.GetOutput())
        multiblock.SetBlock(process, reader.GetOutput())
        multiblock.GetMetaData(process).Set(vtkCompositeDataSet.NAME(), f'Block {process}')

    # write the merged vtk file to disk
    writer.SetFileName(new_filename)
    #writer.SetInputData(multipiece)
    writer.SetInputData(multiblock)
    writer.Write()

# delete vtk files after writing
# for iteration in files_dict:
#     for process, filename in files_dict[iteration].items():
#         os_remove(filename)


