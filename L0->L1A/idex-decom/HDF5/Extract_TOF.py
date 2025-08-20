import h5py
import numpy as np

# Path to your HDF5 file
hdf5_file_path = "ois_output_12182023_184030.h5"
output_file_path = "tof_h_values.txt"  # Change to .csv if needed

# Open the HDF5 file
with h5py.File(hdf5_file_path, "r") as f:
    # Access the TOF H dataset (adjust the group/dataset name as needed)
    tof_h_data = f["14/TOF H"][:]
    
# Save the data to a text file
np.savetxt(output_file_path, tof_h_data, fmt="%d")

print(f"TOF H values saved to {output_file_path}")
