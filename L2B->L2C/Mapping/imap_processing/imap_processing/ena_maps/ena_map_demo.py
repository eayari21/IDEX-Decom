#!/usr/bin/env python
# coding: utf-8

# # Demo of ENA Map Creation tools 
# as it might be used for L2 map creation
# 
# Date: 2025-03-25
# 
# <!-- Selected subset of code to show/demo -- 2025-03-25 -->

# In[1]:


import healpy as hp_vis
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

from imap_processing.ena_maps import ena_maps
from imap_processing.tests.ultra.test_data.mock_data import (
    mock_l1c_pset_product_healpix,
    mock_l1c_pset_product_rectangular,
)
from imap_processing.spice import geometry


# ### This is a rough sketch of some code which makes some synthetic PSETs
# 
# I'd mostly like to ignore the stuff going into the creation of these. There are 2 groups of PSETs:
# 
# 1. Rectangular (lon/lat) gridded PointingSets, covering the entire sky with spacing = 0.5˚.
# 2. Healpix tessellated PointingSets, covering the entire sky with nside = 128 .
# 
# Both contain integer counts pulled from some binomial distribution. To distinguish them, there is a moving bright stripe centered on a latitude.

# ## For portability, manually furnish the necessary SPICE kernels:
# 
# Paths assume that this Notebook is running from `imap_processing/imap_processing/ena_maps/`. You can check that below:

# In[2]:


# get_ipython().system('pwd')


# In[3]:


import spiceypy as spice
spice.furnsh(
    (
        "../tests/spice/test_data/imap_sclk_0000.tsc",
        "../tests/spice/test_data/naif0012.tls",
        "../tests/spice/test_data/imap_science_0001.tf",
        "../tests/spice/test_data/sim_1yr_imap_pointing_frame.bc",
    )
)


# In[ ]:


# Make fake L1c products - both the new healpix style of Ultra PSETs,
# and the old rectangular style

# Spatial Parameters for the fake L1c products
l1c_nside = 128
l1c_spatial_bin_spacing_deg = 0.5

fake_l1c_products_ultra = [
            mock_l1c_pset_product_healpix(
                nside=l1c_nside,
                stripe_center_lat=mid_latitude,
                width_scale=5,
                counts_scaling_params=(50, 0.5),
                peak_exposure=1000,
                timestr=f"2025-05-{4 * i + 1:02d}T12:00:00",
                head=("90"),
            )
            for i, mid_latitude in enumerate(
                np.arange(
                    -90,
                    90,
                    22.5,
                )
            )
        ]

fake_l1c_products_rect = [
            mock_l1c_pset_product_rectangular(
                spacing_deg=l1c_spatial_bin_spacing_deg,
                stripe_center_lat=mid_latitude,
                head="90",
                timestr=f"2025-05-{4 * i + 1:02d}T12:00:00",
            )
            for i, mid_latitude in enumerate(
                np.arange(
                    -90,
                    90,
                    22.5,
                )
            )
        ]


# ### Create the `PointingSet` objects from the synthetic PSETs

# In[ ]:


ultra_psets = [
    ena_maps.UltraPointingSet(
        spice_reference_frame=geometry.SpiceFrame.IMAP_DPS, l1c_dataset=l1c_product,
    )
    for l1c_product in fake_l1c_products_ultra
]

rect_psets = [
    ena_maps.RectangularPointingSet(
        spice_reference_frame=geometry.SpiceFrame.IMAP_DPS, l1c_dataset=l1c_product,
    )
    for l1c_product in fake_l1c_products_rect
]


# #### And take a look at one of them:

# In[ ]:


print(ultra_psets[0].data)


# Plot the middle PSET's counts

# In[ ]:


rect_psets[len(rect_psets)//2].data['counts'].mean(
    dim=("epoch", "energy")).plot(x="longitude", y="latitude")
plt.title("Middle Rectangular PSET Counts")
hp_vis.mollview(
    ultra_psets[len(ultra_psets)//2].data['counts'].mean(
    dim=("epoch", "energy")),
    title="Middle Healpix PSET Counts"
)


# And plot a very rough view of all the counts, naively summed across all the PSETs

# In[ ]:


plt.imshow(np.array([rps.data['counts'].mean(dim=("epoch", "energy")).values for rps in rect_psets]).sum(axis=0).T)
plt.title("Naive sum of all Rectangular PSETS")

hp_vis.mollview(np.array([ups.data['counts'].mean(dim=("epoch", "energy")).values for ups in ultra_psets]).sum(axis=0).T, title="Naive sum of all Healpix PSETS")


# ## Make Maps and Project the Pointing Sets 
# Make 4 total maps: 2 healpix and 2 rectangular - and project each type of PSET to each types of map

# In[ ]:


hp_map_from_hp_pset = ena_maps.HealpixSkyMap(
    nside=64, spice_frame=geometry.SpiceFrame.ECLIPJ2000)
hp_map_from_rect_pset = ena_maps.HealpixSkyMap(
    nside=64, spice_frame=geometry.SpiceFrame.ECLIPJ2000)
rect_map_from_hp_pset = ena_maps.RectangularSkyMap(
    spacing_deg=2, spice_frame=geometry.SpiceFrame.ECLIPJ2000)
rect_map_from_rect_pset = ena_maps.RectangularSkyMap(
    spacing_deg=2, spice_frame=geometry.SpiceFrame.ECLIPJ2000)

# Project the Healpix tessellated Ultra PSETs to each type of SkyMap
for ultra_pset in ultra_psets:
    for skymap in [hp_map_from_hp_pset, rect_map_from_hp_pset]:
        skymap.project_pset_values_to_map(
            ultra_pset, ["counts", "exposure_time", "sensitivity"],
            index_match_method=ena_maps.IndexMatchMethod.PULL)

# Project the Rectangularly tiled PSETS to each type of SkyMap
for rect_pset in rect_psets:
    for skymap in [hp_map_from_rect_pset, rect_map_from_rect_pset]:
        skymap.project_pset_values_to_map(
            rect_pset, ["counts", "exposure_time", "sensitivity"],
            index_match_method=ena_maps.IndexMatchMethod.PULL)


# ### Produce an output dataset and take a look

# In[3]:


hp_map_from_hp_pset_ds = hp_map_from_hp_pset.to_dataset()
pset_epochs = np.array([ps.epoch for ps in ultra_psets])

# The first pset's epoch should be the map's epoch
np.testing.assert_equal(pset_epochs[0], hp_map_from_hp_pset_ds.epoch)

# Let's take a look at the output dataset
hp_map_from_hp_pset_ds


# Plot the `HealpixSkyMap`'s counts

# In[ ]:


hp_vis.mollview(
    hp_map_from_hp_pset_ds["counts"].mean(dim=("epoch", "energy")),
    title="HealpixSkyMap Counts from Healpix UltraPointingSets"
    )
hp_vis.mollview(
    hp_map_from_rect_pset.to_dataset()["counts"].mean(dim=("epoch", "energy")),
    title="HealpixSkyMap Counts from RectangularPointingSets"
)


# Plot the `RectangularSkyMap`'s counts

# In[2]:


var_to_plot = "counts"
for (name, rect_map) in zip(
    ["RectangularSkyMap from Healpix UltraPointingSets",
     "RectangularSkyMap from RectangularPointingSets"],
    [rect_map_from_hp_pset, rect_map_from_rect_pset]
):
    rect_map.to_dataset()[var_to_plot].mean(dim=("epoch", "energy")).plot(x="longitude", y="latitude")
    plt.title(name)
    plt.show()


# In[ ]:




