import pandas as pd
from glob import glob
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sb
import matplotlib
import numpy as np

sb.set_style("whitegrid")
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.sans-serif'] = "Droid Sans"

reds = ["#ffbaba", "#a70000"] # communication colours: light red (TRGTOL) and dark red (TRLTOM)
blues = ["#bac2ff", "#00059f"] # computation colors: light blue (FTDIR) and dark blue (LTDIR)

# Load files for all tasks
df_list = []
for i, file in enumerate(sorted(glob("../build/fort.*"))):
    df = pd.read_csv(file, delim_whitespace=True, names=["event", "step", "batch", "stage", "time"])
    df["task"] = i + 1 # adjust for Fortran counting (starts at 1)
    df_list.append(df)

if len(df_list) == 0:
    raise SystemExit("No files found")

# Compute min and max SYSTEM_CLOCK times, for adjusting x-axis
min_time = min([df["time"].min() for df in df_list])
max_time = max([df["time"].max() for df in df_list])

# Get number of tasks and batches
num_tasks = len(df_list)
num_batches = df_list[0]["batch"].nunique()

# Set up axes
fig, ax = plt.subplots(figsize=(10, 1.6*num_tasks))

for df_task in df_list:
    for event, colors in zip(["COMM", "COMP"], [reds, blues]):
        # Only loop over the opening events (step = 1)
        for index, row in df_task[(df_task["event"] == event) & (df_task["step"] == 1)].iterrows():
            batch = row["batch"]
            stage = row["stage"]
            start_time = 1000.0*(row["time"] - min_time) # Compute start time of this event in ms
            task = row["task"]

            # Find matching closing event (step = 3)
            end_row_bool = (df_task["event"] == event) & (df_task["step"] == 3) & \
                           (df_task["batch"] == batch) & (df_task["stage"] == stage)

            # Compute duration of this event in ms
            duration = 1000.0*(df_task[index:][end_row_bool].iloc[0]["time"] - min_time) - start_time

            # Add rectangle to plot for this event
            rect = patches.Rectangle((start_time, (task - 1) * num_batches + batch), duration, 1.0,
                                    facecolor=colors[stage-1])
            ax.add_patch(rect)

# Plot task labels
task_label_x = 1.13 * 1000.0 * (max_time - min_time)
for task in range(num_tasks):
    ax.plot([task_label_x, task_label_x], [task * num_batches + 1, (task + 1) * num_batches],
            marker='o', clip_on=False, color="silver")
    ax.text(task_label_x*1.018, task * num_batches + 1 + num_batches/2, f"Task {task+1}",
            horizontalalignment="center", verticalalignment="center", rotation="vertical")

# Plot phony data to allow legend
ax.plot(1000 * task_label_x, 0, label="TRGTOL", color=reds[0])
ax.plot(1000 * task_label_x, 0, label="TRLTOM", color=reds[1])
ax.plot(1000 * task_label_x, 0, label="FTDIR", color=blues[0])
ax.plot(1000 * task_label_x, 0, label="LTDIR", color=blues[1])
ax.legend(ncol=4, loc="center", bbox_to_anchor=(0.5, -0.01), frameon=False)

plt.xlim(0.0, 1.1 * 1000.0 * (max_time - min_time))
plt.ylim(1.0, num_tasks * num_batches + 1)

ax.set_yticks(np.arange(1.5, num_tasks * num_batches + 1))
labels = [item.get_text() for item in ax.get_yticklabels()]
labels = list(range(1, num_batches+1)) * num_tasks
ax.set_yticklabels(labels)
ax.set_ylabel("Batch")

ax.set_xlabel("Time (ms)")

plt.savefig("trace.png", bbox_inches="tight")
plt.savefig("trace.pdf", bbox_inches="tight")
