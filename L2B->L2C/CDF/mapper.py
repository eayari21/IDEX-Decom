import numpy as np
import matplotlib.pyplot as plt
from spacepy import pycdf

# Load the CDF
cdf = pycdf.CDF("imap_idex_l2c_rectangular-map-1week_20231218_v999.cdf")

# Pull out counts — shape is (1, N_lat, N_lon)
counts = np.array(cdf['counts'][0])  # shape (N_lat, N_lon)
n_lat, n_lon = counts.shape

# Synthetic bin edges (can be replaced by cdf['rectangular_lat_pixel_label'] etc.)
lat_min, lat_max = -90, 90
lon_min, lon_max = -180, 180
lat_edges = np.linspace(lat_min, lat_max, n_lat + 1)
lon_edges = np.linspace(lon_min, lon_max, n_lon + 1)

lon_grid, lat_grid = np.meshgrid(lon_edges, lat_edges)

# === Plot ===
fig, ax = plt.subplots(figsize=(12, 10), facecolor='black')
fig.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.15)

pcm = ax.pcolormesh(lon_grid, lat_grid, counts, shading='auto', cmap='rainbow', rasterized=True)

cb = plt.colorbar(pcm, ax=ax, orientation='horizontal', pad=0.15)
cb.set_label("Counts", font="Helvetica", fontsize=14, color='white')
cb.ax.tick_params(colors='white')

# Style like IBEX map
ax.set_title("IMAP-IDEX Map", fontsize=16, color='white', pad=20)
ax.set_xlabel("Longitude [deg]", fontsize=16, color='white')
ax.set_ylabel("Latitude [deg]", fontsize=16, color='white')
ax.set_xlim([-180, 180])
ax.set_ylim([-90, 90])
ax.set_aspect('equal')

ax.tick_params(colors='white')
ax.grid(color='green', linestyle='--', linewidth=0.5, alpha=0.5)

plt.show()
