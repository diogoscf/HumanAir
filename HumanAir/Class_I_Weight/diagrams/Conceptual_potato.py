import numpy as np
import matplotlib.pyplot as plt

cases = ["OEW", "OEW + PL", "OEW + PL + FL", "OEW + FL"]
case_numbers = [1, 2, 3, 4]

# conventional aircraft
xcg_conv = np.array([3.06, 4.18, 4.13, 3.10])  # m from nose
mass_conv = np.array([0, 675, 822.6, 147.6])  # kg

# canard
xcg_canard = np.array([6.15, 5.82, 5.95, 6.32])
mass_canard = np.array([0, 675, 800.9, 128.9])  # kg

# flying wing
xcg_fwing = np.array([1.92, 1.75, 1.74, 1.9])
mass_fwing = np.array([0, 675, 888.3, 213.3])  # kg

# Append the first points to the end to close the loop
xcg_conv_closed = np.append(xcg_conv, xcg_conv[0])
mass_conv_closed = np.append(mass_conv, mass_conv[0])

xcg_canard_closed = np.append(xcg_canard, xcg_canard[0])
mass_canard_closed = np.append(mass_canard, mass_canard[0])

xcg_fwing_closed = np.append(xcg_fwing, xcg_fwing[0])
mass_fwing_closed = np.append(mass_fwing, mass_fwing[0])

# Create the subfigures
fig, ax = plt.subplots(1, 3, figsize=(15, 5))
# Plotting for conventional aircraft
ax[0].plot(xcg_conv_closed, mass_conv_closed, "-o", label="_nolegend_")
for i, num in enumerate(case_numbers):
    color = "red" if xcg_conv[i] == min(xcg_conv) or xcg_conv[i] == max(xcg_conv) else "black"
    ax[0].plot(xcg_conv[i], mass_conv[i], "o", color=color)
    ax[0].annotate(num, (xcg_conv[i], mass_conv[i]), textcoords="offset points", xytext=(5, 5), ha="center")
ax[0].grid()
ax[0].set_title("Conventional")

# Custom legend for conventional aircraft
min_xcg_conv = np.min(xcg_conv)
max_xcg_conv = np.max(xcg_conv)
legend_elements_conv = [
    plt.Line2D(
        [0],
        [0],
        marker="",
        color="w",
        label=f"Fwd $X_{{cg}}$: {min_xcg_conv} ({cases[int(np.nonzero(xcg_conv == np.min(xcg_conv))[0])]})",
        markerfacecolor="red",
        markersize=10,
    ),
    plt.Line2D(
        [0],
        [0],
        marker="",
        color="w",
        label=f"Aft $X_{{cg}}$: {max_xcg_conv} ({cases[int(np.nonzero(xcg_conv == np.max(xcg_conv))[0])]})",
        markerfacecolor="red",
        markersize=10,
    ),
]
ax[0].legend(handles=legend_elements_conv)

# Plotting for canard
ax[1].plot(xcg_canard_closed, mass_canard_closed, "-o", label="_nolegend_")
for i, num in enumerate(case_numbers):
    color = "red" if xcg_canard[i] == min(xcg_canard) or xcg_canard[i] == max(xcg_canard) else "black"
    ax[1].plot(xcg_canard[i], mass_canard[i], "o", color=color)
    ax[1].annotate(num, (xcg_canard[i], mass_canard[i]), textcoords="offset points", xytext=(5, 5), ha="center")
ax[1].grid()
ax[1].set_title("Canard")

# Custom legend for canard
min_xcg_canard = np.min(xcg_canard)
max_xcg_canard = np.max(xcg_canard)
legend_elements_canard = [
    plt.Line2D(
        [0],
        [0],
        marker=" ",
        color="w",
        label=f"Fwd $X_{{cg}}$: {min_xcg_canard} ({cases[int(np.nonzero(xcg_canard == np.min(xcg_canard))[0])]})",
        markerfacecolor="red",
        markersize=10,
    ),
    plt.Line2D(
        [0],
        [0],
        marker=" ",
        color="w",
        label=f"Aft $X_{{cg}}$: {max_xcg_canard} ({cases[int(np.nonzero(xcg_canard == np.max(xcg_canard))[0])]})",
        markerfacecolor="red",
        markersize=10,
    ),
]
ax[1].legend(handles=legend_elements_canard)

# Plotting for flying wing
ax[2].plot(xcg_fwing_closed, mass_fwing_closed, "-o", label="_nolegend_")
for i, num in enumerate(case_numbers):
    color = "red" if xcg_fwing[i] == min(xcg_fwing) or xcg_fwing[i] == max(xcg_fwing) else "black"
    ax[2].plot(xcg_fwing[i], mass_fwing[i], "o", color=color)
    ax[2].annotate(num, (xcg_fwing[i], mass_fwing[i]), textcoords="offset points", xytext=(5, 5), ha="center")
ax[2].grid()
ax[2].set_title("Flying Wing")

# Custom legend for flying wing
min_xcg_fwing = np.min(xcg_fwing)
max_xcg_fwing = np.max(xcg_fwing)
legend_elements_fwing = [
    plt.Line2D(
        [0],
        [0],
        marker=" ",
        color="w",
        label=f"Fwd $X_{{cg}}$: {min_xcg_fwing} ({cases[int(np.nonzero(xcg_fwing == np.min(xcg_fwing))[0])]})",
        markerfacecolor="red",
        markersize=10,
    ),
    plt.Line2D(
        [0],
        [0],
        marker=" ",
        color="w",
        label=f"Aft $X_{{cg}}$: {max_xcg_fwing} ({cases[int(np.nonzero(xcg_fwing == np.min(xcg_fwing))[0])]})",
        markerfacecolor="red",
        markersize=10,
    ),
]
ax[2].legend(handles=legend_elements_fwing)

# Set common labels
fig.text(0.5, 0.02, "$X_{cg} [m]$", ha="center", fontsize=12)
fig.text(0.002, 0.5, "Mass [kg]", va="center", rotation="vertical", fontsize=12)

# Legend
handles = [
    plt.Line2D([0], [0], marker="", color="w", label=f"{num}: {case}", markerfacecolor="black", markersize=10)
    for num, case in zip(case_numbers, cases)
]
fig.legend(handles=handles, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1))

plt.tight_layout(rect=(0, 0.03, 1, 0.95))
plt.savefig("xcg_range.png")
plt.show()
