#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import h5py


def main():
    # Use a filename provided as command line argument or default to "data.txt"
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "idex_test_data_lon_lat_degree.txt"
    
    # Read the file.
    # Assumes the file has three columns: event_number, latitude, longitude.
    data = pd.read_csv(filename, delim_whitespace=True, header=None, names=["event", "lat", "lon"])
    
    # For debugging, print first few rows.
    print("Data preview:")
    print(data.head())

    # Create an xarray.Dataset.
    # Here, we use the CSV row index as the 'healpix_index' coordinate.
    ds = xr.Dataset(
        data_vars={
            "event": (("healpix_index",), data["event"].values),
            "latitude": (("healpix_index",), data["lat"].values),
            "longitude": (("healpix_index",), data["lon"].values)
        },
        coords={
            "healpix_index": data.index.values
        }
    )
    
    # Show a summary of the created Dataset
    print("\nCreated xarray.Dataset:")
    print(ds)

    hdf5_filename = "imap_l2c_example.h5"

    with h5py.File(hdf5_filename, "w") as f:
        for idx in ds.healpix_index.values:
            # Create a group named after the current healpix index
            grp = f.create_group(str(idx))
            
            # Extract scalar values for each variable at this index
            event_value = ds.event.sel(healpix_index=idx).values.item()
            latitude_value = ds.latitude.sel(healpix_index=idx).values.item()
            longitude_value = ds.longitude.sel(healpix_index=idx).values.item()
            
            # Create datasets within the group
            grp.create_dataset("event", data=event_value)
            grp.create_dataset("latitude", data=latitude_value)
            grp.create_dataset("longitude", data=longitude_value)
            
    print(f"Dataset written to {hdf5_filename}")
    
    
    # Create a 2D histogram: X axis -> longitude, Y axis -> latitude.
    # You can adjust the number of bins (here, 50 in each dimension)
    heatmap, xedges, yedges = np.histogram2d(data["lat"], data["lon"], bins=20)
    
    # Plot the heat map
    plt.figure(figsize=(8, 6))
    plt.imshow(heatmap.T, origin='lower',
               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
               aspect='auto', cmap='hot')
    plt.colorbar(label="Binned Count")
    # We change the fontsize of minor ticks label 
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.tick_params(axis='both', which='minor', labelsize=20)
    plt.xlabel("Longitude [°]", fontsize=40, font="Helvetica")
    plt.ylabel("Latitude [°]", fontsize=40, font="Helvetica")
    plt.title("Heat Map of Binned Counts")
    plt.show()

if __name__ == "__main__":
    main()
