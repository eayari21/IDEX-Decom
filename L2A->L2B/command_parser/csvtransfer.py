import pandas as pd
import csv
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, DateFormatter
from datetime import datetime, timedelta
from collections import defaultdict


# Load CSV
df = pd.read_csv("oasis_events.csv")

# Create empty list for results
results = []

# Iterate over each row and apply the filter logic
for _, row in df.iterrows():
    message = str(row["message"])
    date_str = row["timeGps (yyyy-D-'T'HH:mm:ss)"]

    # Parse "2025-96-T00:02:28" into "2025-096-00:02:28"
    try:
        parts = date_str.split("-T")
        year, doy = parts[0].split("-")
        time = parts[1]
        padded_date = f"{year}-{int(doy):03d}-{time}"
        date_obj = datetime.strptime(padded_date, "%Y-%j-%H:%M:%S")
    except Exception as e:
        print(f"⚠️ Skipping unparsable date: {date_str} ({e})")
        continue

    if "ACQSETUP) ==&#62; sciState16Dictionary(ACQ)" in message:
        results.append({"date": date_obj, "command": "on"})
    elif "ACQ) ==&#62; sciState16Dictionary(CHILL)" in message:
        results.append({"date": date_obj, "command": "off"})

# Convert results to DataFrame
parsed_df = pd.DataFrame(results)
parsed_df.to_csv("IDEX_DOY_94-116_2025.csv")

# Show preview
print(f"Parsed_df has {len(parsed_df)} elements.")
print(parsed_df.head())

# Skip if empty
if parsed_df.empty:
    print("🚫 No transitions found.")
    exit()

# Sort by datetime
parsed_df = parsed_df.sort_values("date").reset_index(drop=True)

# Build segments for plotting
segments = []
for i in range(len(parsed_df) - 1):
    start = parsed_df.loc[i, "date"]
    end = parsed_df.loc[i + 1, "date"]
    state = parsed_df.loc[i, "command"]
    color = "lightblue" if state == "on" else "lightgray"
    segments.append((date2num(start), date2num(end), color))

# Plotting
fig, ax = plt.subplots(figsize=(10, 1.5))
for start, end, color in segments:
    ax.axvspan(start, end, ymin=0, ymax=1, color=color)

ax.set_yticks([])
ax.set_xlabel("Date")
ax.set_title("ACQ State Timeline (Light Blue = ON, Light Gray = OFF)")
ax.set_xlim(date2num(parsed_df["date"].min()), date2num(parsed_df["date"].max()))

# ✅ Fix weird date labels
ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d %H:%M'))
fig.autofmt_xdate()

plt.tight_layout()
plt.show()

# Track total and 'on' durations per day
daily_totals = defaultdict(timedelta)
daily_on = defaultdict(timedelta)

for i in range(len(parsed_df) - 1):
    start = parsed_df.loc[i, "date"]
    end = parsed_df.loc[i + 1, "date"]
    state = parsed_df.loc[i, "command"]

    # Split time span by day boundaries
    current = start
    while current < end:
        next_day = (current + timedelta(days=1)).replace(hour=0, minute=0, second=0, microsecond=0)
        segment_end = min(end, next_day)

        duration = segment_end - current
        day_str = current.date().isoformat()

        daily_totals[day_str] += duration
        if state == "on":
            daily_on[day_str] += duration

        current = segment_end

# Sort and print results
print("\n📊 Percentage of time in 'on' state per day:\n")
for day in sorted(daily_totals.keys()):
    total = daily_totals[day].total_seconds()
    on_time = daily_on[day].total_seconds()
    pct_on = (on_time / total) * 100 if total > 0 else 0
    print(f"{day}: {pct_on:.2f}%")


# Save to CSV
with open("on_percentage_by_day.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Date", "Percent On"])
    for day in sorted(daily_totals.keys()):
        total = daily_totals[day].total_seconds()
        on_time = daily_on[day].total_seconds()
        pct_on = (on_time / total) * 100 if total > 0 else 0
        writer.writerow([day, f"{pct_on:.2f}"])