#!/opt/anaconda3/bin/python3
# -*- coding: utf-8 -*-

title = """

//  ================================================================================
//  ||                                                                            ||
//  ||              backup_science                                                ||
//  ||              ------------------------------------------------------        ||
//  ||                           B A C K U P   S C I E N C E                      ||
//  ||              ------------------------------------------------------        ||
//  ||                                                                            ||
//  ||                __author__      = Ethan Ayari                               ||
//  ||                IMPACT/LASP, CU Boulder                                     ||
//  ||                                                                            ||
//  ||                For: IDEX Flight Model Integration and Test, L0-L1A         ||
//  ||                                                                            ||
//  ||                2023                                                        ||
//  ||                                                                            ||
//  ||                                                                            ||
//  ||                Works with Python 3.10.4                                    ||
//  ||                                                                            ||
//  ================================================================================


backup_science: A backup script for the science tool.

"""

# %%
# || Python libraries
import subprocess
import shutil
import os

print(title)
# %%
# || OS-independent most recent file/folder selector
def get_most_recent_file_or_folder(directory):
    # Get a list of all files and folders in the directory
    items = os.listdir(directory)

    # Sort the items based on their creation time (mtime)
    items.sort(key=lambda item: os.path.getmtime(os.path.join(directory, item)))

    # Get the most recently created item (last item in the sorted list)
    most_recent_item = items[-1] if items else None

    return most_recent_item

# %%
# || OS-independent tar command
def compress_and_move_folder(targetFolder):
    # Compress the folder using tar
    tar_command = ["tar", "-cvzf", f"{targetFolder}.tar.gz", targetFolder]
    subprocess.run(tar_command, check=True)

    # Move the compressed file to the "Outgoing" directory
    mv_command = ["mv", f"{targetFolder}.tar.gz", "Outgoing"]
    subprocess.run(mv_command, check=True)

# || Executable code
if __name__ == "__main__":
    # Get the ois output filename
    fname = get_most_recent_file_or_folder("output_from_ois")

    if(os.path.exists(os.path.join(os.getcwd(),"Outgoing"))):
       pass
    else:
        os.makedirs(os.path.join(os.getcwd(),"Outgoing"))

    # Compress and move the plot folder
    PlotFolder = os.path.join(os.getcwd(), f"Plots/{fname}")
    compress_and_move_folder(PlotFolder)
    # Copy over the HDF5 file
    shutil.copy(os.path.join(os.getcwd(), f"HDF5/{fname}.h5"), "Outgoing")



