# IDEXScienceTool

A Python library dedicated to the decompressing, decommutating, storing and plotting of IDEX science packets.


Author: Ethan Ayari, Institute for Modeling Plasma, Atmosphere and Cosmic Dust

## Getting started

To get started, please import my pip environment:

### pip install -r idex.txt

Also, if the "output_from_ois" directory is not present, please create it in the same directory as the repository:

### mkdir output_from_ois

## Script usage

The science tool is a combination of the read_from_ois+idex_packet scripts. It is an all-in-one tool to retrieve, store, and plot IDEX packets from the OASIS server.

To test the science tool, please issue the command:  

### ./science_tool.py

To use the quicklook tools, the following command will prompt you to select an HDF5 file:

### ./IDEX-quicklook.py

To decommutate, store and plot individual packet files, issue the command (replace Data/ois_output_05182023_204405 with your output file):

### ./idex_packet.py --file Data/ois_output_05182023_204405

or equivalently,

### ./idex_packet.py -f Data/ois_output_05182023_204405

This will populate the "Plots" and "HDF5" folders with the corresponding plots and data.

In order to decompress data, issue the following command prior to running idex_packet.py on it:

### ./rice_decode.py --sourcefile sciData_2023_102_14_24_55 --targetfile uncompressed_sciData_2023_102_14_24_55   

If the science tool has failed for some reason, use "read_from_ois.py" to generate a packet file, and "idex_packet.py" to decommutate it and plot it. 
