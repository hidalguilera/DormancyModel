import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Fig size for two-column article
fontsize = 10
fig_width = 3.375
fig_height = fig_width * (4/6)

plt.rcParams.update({
    "figure.figsize": (fig_width, fig_height),
    "font.size": fontsize,
    "axes.labelsize": fontsize,
    "xtick.labelsize": fontsize * 0.8,
    "ytick.labelsize": fontsize * 0.8,
    "legend.fontsize": fontsize * 0.8,
    "text.usetex": False,  # Use mathtext instead of LaTeX
    "mathtext.fontset": "cm",  # Computer Modern math fonts
    "font.family": "serif",
    "font.serif": ["Computer Modern"],
})

# Read and combine data
data1 = pd.read_csv("parallel_time_correlation_individual_replicates_2025-11-14_17-12-16.csv")
data2 = pd.read_csv("parallel_time_correlation_individual_replicates_2025-11-18_15-41-45.csv")
data = pd.concat([data1, data2], ignore_index=True)

# Group and aggregate
comb = data.groupby(['sigma', 'init_a']).agg(
    mean_final_awake_rate=('final_awake_rate', lambda x: (np.mean(np.log10(x)))),
    std_final_awake_rate=('final_awake_rate', lambda x: (np.std(np.log10(x))))
).reset_index()

# Filter sigma > 89
filtered = comb[comb['sigma'] > 89]


# Create scatter plot
fig, ax = plt.subplots()
ax.errorbar(filtered['init_a'], filtered['final_awake_rate'], 
            fmt='o', ms=6, capsize=2.5, alpha=0.6)

ax.set_xlabel(r'Initial awake rate')
ax.set_ylabel(r'Mean final awake rate')
ax.set_xscale('log')
ax.set_yscale('log')
fig.set_figwidth(3)
fig.set_figheight(3)

plt.tight_layout()
plt.savefig('scatter_plot_all2.pdf', dpi=300, bbox_inches='tight')