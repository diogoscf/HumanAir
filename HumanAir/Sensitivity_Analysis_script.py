import matplotlib.pyplot as plt

# Given data without Parameter 3
percentages = [
    -15,
    -14,
    -13,
    -12,
    -11,
    -10,
    -9,
    -8,
    -7,
    -6,
    -5,
    -4,
    -3,
    -2,
    -1,
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
]
data_arrays = [
    [
        0.34,
        0.34,
        0.33,
        0.32,
        0.32,
        0.31,
        0.3,
        0.3,
        0.29,
        0.28,
        0.28,
        0.27,
        0.27,
        0.26,
        0.25,
        0.25,
        0.24,
        0.23,
        0.23,
        0.22,
        0.22,
        0.21,
        0.2,
        0.2,
        0.19,
        0.18,
        0.18,
        0.17,
        0.17,
        0.16,
    ],
    [
        0.63,
        0.63,
        0.63,
        0.62,
        0.62,
        0.62,
        0.62,
        0.61,
        0.61,
        0.61,
        0.6,
        0.6,
        0.6,
        0.6,
        0.59,
        0.59,
        0.59,
        0.58,
        0.58,
        0.58,
        0.58,
        0.57,
        0.57,
        0.57,
        0.56,
        0.56,
        0.56,
        0.56,
        0.55,
        0.55,
    ],
    [
        4.712478248031918,
        4.712478248031918,
        4.712478248031918,
        4.712478248031918,
        4.712478248031918,
        4.712478248031918,
        4.712478248031918,
        4.712478248031918,
        4.712482827732589,
        4.712482827732589,
        4.712482827732589,
        4.712482827732589,
        4.712482827732589,
        4.712482827732589,
        4.712482827732589,
        4.712482827732589,
        4.95982513711783,
        5.9597399500498565,
        5.9597399500498565,
        5.9597399500498565,
        5.9597399500498565,
        5.9597399500498565,
        5.9597399500498565,
        5.9597399500498565,
        5.9597399500498565,
        5.9597399500498565,
        5.836689271448258,
        5.9597399500498565,
        5.836689271448258,
        5.9597399500498565,
    ],
    [
        202.0289258979847,
        202.57908172144064,
        203.12923754489663,
        203.67939336835258,
        204.22954919180853,
        204.7797050152645,
        205.3298608387205,
        205.88001666217644,
        206.43017248563243,
        206.98032830908838,
        207.53048413254436,
        208.0806399560003,
        208.6307957794563,
        209.18095160291224,
        209.73110742636823,
        210.28126324982418,
        210.83141907328016,
        211.38157489673614,
        211.9317307201921,
        212.48188654364805,
        213.03204236710403,
        213.58219819055998,
        214.13235401401596,
        214.6825098374719,
        215.2326656609279,
        215.78282148438385,
        216.33297730783983,
        216.88313313129578,
        217.43328895475173,
        217.9834447782077,
    ],
]

# Find values at 0 percent
zero_percent_values = [arr[percentages.index(0)] for arr in data_arrays]

# Calculate relative differences in percentages
relative_differences = [
    [(value - zero_percent_values[i]) / zero_percent_values[i] * 100 for value in data_arrays[i]]
    for i in range(len(data_arrays))
]

# Plotting
plt.figure(figsize=(20, 6))

# Plot each parameter's relative difference with specific labels
labels = ["Current Technology CO2 Reduction", "Future Technology CO2 Reduction", "Sh", "Operation Cost"]
for i, rel_diff in enumerate(relative_differences):
    plt.plot(percentages, rel_diff, label=labels[i])

plt.axhline(0, color="black", linestyle="--", linewidth=0.5)  # Add a horizontal line at y=0 for reference

plt.xlabel("Percentage Variation in MTOW", fontsize=20)
plt.ylabel("Relative Difference (%)", fontsize=20)

plt.legend(fontsize=20)
plt.grid(True)
plt.tight_layout()
plt.savefig("Sensitivity_Analysis_general.pdf")
plt.show()
