from matplotlib import pyplot as plt
from os import listdir
from numpy import array
from sys import argv

def create_file_dict(filename):
    data_dict = {}
    with open(filename, 'r') as file:
        # Read the first line to get the column labels
        columns = file.readline().split()
        # Initialize dictionary values as empty lists
        for column in columns:
            key = column
            data_dict[key] = []
        # Read the rest of the lines
        for line in file:
            values = line.split()
            step_val = values[0]
            data_dict["step"].append(int(step_val))
            mem_val = values[-1]
            data_dict["mem"].append(int(mem_val))
            # concatenates the rest for the description
            desc_val = values[1]
            for val in values[2:-1]:
                desc_val = desc_val + " " + val 
            data_dict["desc"].append(desc_val)
            
    # Convert lists to NumPy arrays
    for key in data_dict:
        if key != "desc":
            data_dict[key] = array(data_dict[key])
    return data_dict

# folder = "newt-ref2_coarse/save/"
# get folder as argument from the user
folder = argv[1]
filename = folder + "/save/mem.txt"
data = create_file_dict(filename)

fig, ax = plt.subplots(figsize=(10, 5))
title = "memory usage"
ax.set_title(title)
ax.set_xlabel("step")
ax.set_ylabel("memory usage (MB)")
# y-axis starting at zero
ax.set_ylim(bottom=0, top=max(data["mem"]/1000)*1.1)

ax.plot(data["step"], data["mem"]/1000, linestyle='-', color='b', label='memory usage (MB)')
ax.legend(loc='lower right')

# create high resolution png of the plot
plt.savefig(f"{folder}/memory.png", dpi=300)