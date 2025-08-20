import os
import h5py
import numpy as np
import pandas as pd
from datetime import datetime, timedelta

# Path to the directory containing your HDF5 files
data_dir = "/Users/etay8828/Desktop/IDEX_Cal/idex-decom/HDF5/Flight_Packets_Analysis_3"

# Bin edges
mass_edges = np.array([
    6.31e-17, 1.00e-16, 1.58e-16, 2.51e-16, 3.98e-16,
    6.31e-16, 1.00e-15, 1.58e-15, 2.51e-15, 3.98e-15, 1.00e-14
])
charge_edges = np.array([
    1.00e-01, 3.16e-01, 1.00e+00, 3.16e+00, 1.00e+01,
    3.16e+01, 1.00e+02, 3.16e+02, 1.00e+03, 3.16e+03, 1.00e+04
])

# Bin label helpers
def format_mass_bin(i):
    if i == 0:
        return "Mass Bin <1"
    elif i == 10 or i >= 11:
        return "Mass Bin 10*"
    else:
        return f"Mass Bin {i}*"

def format_charge_bin(i):
    if i == 0:
        return "Charge Bin <1"
    elif i == 10 or i >= 11:
        return "Charge Bin 10*"
    else:
        return f"Charge Bin {i}*"

# Prepare filenames and assign dates
all_files = sorted([f for f in os.listdir(data_dir) if f.endswith(".h5")])
start_date = datetime(2025, 4, 6)
num_days = 15
date_bins = [start_date + timedelta(days=i) for i in range(num_days)]
files_per_day = np.array_split(all_files, num_days)

# Create rows for final DataFrame
final_rows = []

for day, files in zip(date_bins, files_per_day):
    mass_counts = {f"Mass Bin {i}*" if i != 1 and i != 10 else f"Mass Bin {i}*": 0 for i in range(1, 11)}
    charge_counts = {f"Charge Bin {i}*" if i != 1 and i != 10 else f"Charge Bin {i}*": 0 for i in range(1, 11)}

    for fname in files:
        fpath = os.path.join(data_dir, fname)
        try:
            with h5py.File(fpath, "r") as h5file:
                for top_level_key in h5file.keys():
                    group = h5file[top_level_key]
                    if "Analysis" in group:
                        analysis = group["Analysis"]
                        try:
                            mass_ds = analysis["Target LMassEstimate"]
                            charge_ds = analysis["Target LImpactCharge"]

                            mass_vals = [mass_ds[()]] if np.isscalar(mass_ds[()]) else mass_ds[:]
                            charge_vals = [charge_ds[()]] if np.isscalar(charge_ds[()]) else charge_ds[:]

                            # Convert units
                            mass_vals = 1e-15 * np.array(mass_vals)  # fg to kg
                            charge_vals = np.array(charge_vals)      # already in pC

                            # Bin
                            mass_bins = np.digitize(mass_vals, bins=mass_edges)
                            charge_bins = np.digitize(charge_vals, bins=charge_edges)

                            for m_bin in mass_bins:
                                label = format_mass_bin(m_bin)
                                if label in mass_counts:
                                    mass_counts[label] += 1

                            for c_bin in charge_bins:
                                label = format_charge_bin(c_bin)
                                if label in charge_counts:
                                    charge_counts[label] += 1
                        except KeyError:
                            continue
        except Exception:
            continue

    row = {"Date": day.strftime("%Y-%m-%d")}
    row.update(mass_counts)
    row.update(charge_counts)
    final_rows.append(row)

# Create final DataFrame
df = pd.DataFrame(final_rows)
df["Date"] = pd.to_datetime(df["Date"]).dt.dayofyear
df.rename(columns={"Date": "DOY"}, inplace=True)

# Display or export
print(df)
df.to_csv("mass_charge_daily_counts.csv", index=False)


# Create a separate copy for rates
df_rates = df.copy()

# Load uptime percentages
on_df = pd.read_csv("on_percentage_by_day.csv")

# Convert Date to datetime and extract DOY
on_df["Date"] = pd.to_datetime(on_df["Date"])
on_df["DOY"] = on_df["Date"].dt.dayofyear

# Merge with df_rates on DOY
df_rates = pd.merge(df_rates, on_df[["DOY", "Percent On"]], on="DOY", how="left")

# Compute rates using the formula
count_columns = [col for col in df.columns if col != "DOY"]
for col in count_columns:
    df_rates[col] = df_rates[col] / (0.01 * df_rates["Percent On"] * 86400)

# Set rates to zero where Percent On is zero
df_rates.loc[df_rates["Percent On"] == 0, count_columns] = -1

# Optional: remove Percent On column if not needed
# df_rates.drop(columns=["Percent On"], inplace=True)
df_rates.to_csv("mass_charge_daily_rates.csv", index=False)

# ✅ `df` has the original counts
# ✅ `df_rates` has the calculated rates

# === INPUTS ===
df = pd.read_csv("mass_charge_daily_counts.csv")      # includes "DOY", "Mass Bin 1*", ..., "Charge Bin 10*"
df_rates = pd.read_csv("mass_charge_daily_rates.csv") # same structure, plus "Percent On"
output_path = "IDEX_L2B_Example.h5"

# === DATE & WOY TAGGING ===
def doy_to_date(doy, year=2025):
    return datetime(year, 1, 1) + timedelta(days=doy - 1)

df["Date"] = df["DOY"].apply(lambda x: doy_to_date(x))
df["WOY"] = df["Date"].apply(lambda x: x.isocalendar().week)

df_rates["Date"] = df["DOY"].apply(lambda x: doy_to_date(x))
df_rates["WOY"] = df_rates["Date"].apply(lambda x: x.isocalendar().week)

# === FORMAT CLEANING ===
def clean_label(colname):
    return colname.lower().replace(" ", "_").replace("*", "")

count_cols = [col for col in df.columns if col.startswith("Mass Bin") or col.startswith("Charge Bin")]
rate_cols = [col for col in df_rates.columns if col.startswith("Mass Bin") or col.startswith("Charge Bin")]

# === HDF5 WRITING ===
with h5py.File(output_path, "w") as h5file:
    for woy in sorted(df["WOY"].unique()):
        woy_group = h5file.create_group(f"WOY_{woy}")

        # Filter by WOY
        df_woy = df[df["WOY"] == woy]
        df_rates_woy = df_rates[df_rates["WOY"] == woy]

        # Add Epoch = first DOY in the week
        first_doy = int(df_woy["DOY"].min())
        woy_group.create_dataset("Epoch", data=first_doy)

        # Loop through each DOY in this week
        for _, row in df_woy.iterrows():
            doy = int(row["DOY"])
            doy_group = woy_group.create_group(f"DOY_{doy}")

            # Write count data
            for col in count_cols:
                dataset_name = f"counts_{clean_label(col)}"
                doy_group.create_dataset(dataset_name, data=[row[col]])

            # Write rate data
            row_rates = df_rates_woy[df_rates_woy["DOY"] == doy].iloc[0]
            for col in rate_cols:
                dataset_name = f"rate_{clean_label(col)}"
                doy_group.create_dataset(dataset_name, data=[row_rates[col]])

            # Compute total counts and assign randomly to spin phase bins
            total_count = sum(row[col] for col in count_cols)
            if total_count > 0:
                # Generate 4 random weights that sum to 1
                weights = np.random.dirichlet(np.ones(4))
                spin_counts = np.random.multinomial(total_count, weights)
            else:
                spin_counts = [0, 0, 0, 0]

            # Store spin phase bin counts
            for i in range(4):
                doy_group.create_dataset(f"spin_phase_bin{i+1}", data=[spin_counts[i]])

print(f"✅ HDF5 file written with per-DOY groups under WOY: {output_path}")