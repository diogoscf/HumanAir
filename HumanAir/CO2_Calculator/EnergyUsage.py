import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle, Patch
import os

phases = ["TO", "CR", "LD"]
# Define the frequency distribution for mission distances (in nautical miles)
distance_bins = np.linspace(0, 500, 11)

# Define frequency for each bin
distance_frequencies = np.ones(10) * 2  # Frequency for each bin

# Power requirements (assumed constant for takeoff and landing, varying for cruise)
P_takeoff = 298  # kW
P_landing = 0.25 * P_takeoff  # kW
P_cruise_per_mile = 0.476  # kW per nautical mile

# Define durations for takeoff and landing
takeoff_duration = 5 / 60  # 5 minutes in hours
landing_duration = 10 / 60  # 10 minutes in hours

# Total energy available
E_battery = 256  # in kWh
E_fuel = 3524  # in kWh


# Helper function to generate missions based on frequency distribution
def generate_missions(distance_bins, distance_frequencies):
    missions = []
    for i in range(len(distance_bins) - 1):
        bin_start = distance_bins[i]
        bin_end = distance_bins[i + 1]
        frequency = distance_frequencies[i]
        interval = np.linspace(bin_start, bin_end, int(frequency) + 1)
        for j in range(int(frequency)):
            distance = (interval[j] + interval[j + 1]) / 2
            cruise_duration = distance / 116.631  # Simplified duration calculation
            cruise_power = 238  # kW (constant power for cruise)
            missions.append(
                [
                    ("TO", takeoff_duration, P_takeoff, 0.31 + distance),
                    ("CR", cruise_duration, cruise_power, distance),
                    ("LD", landing_duration, P_landing, 0.31 + distance),
                ]
            )
    return missions


# Generate missions
missions = generate_missions(distance_bins, distance_frequencies)

# Allocate power sources and track usage for each distance bin
battery_usage_bins = np.zeros(len(distance_bins) - 1)
fuel_usage_bins = np.zeros(len(distance_bins) - 1)

for mission in missions:
    remaining_battery_energy = E_battery
    remaining_fuel_energy = E_fuel
    for segment in mission:
        name, duration, power, distance = segment
        energy_required = duration * power

        # Ensure bin index is within valid range
        bin_index = np.digitize(distance, distance_bins) - 1
        bin_index = min(bin_index, len(distance_bins) - 2)  # Correct the bin index if it goes out of bounds

        # Check if battery energy is sufficient for this segment
        if remaining_battery_energy >= energy_required:
            battery_usage_bins[bin_index] += energy_required
            remaining_battery_energy -= energy_required
        else:
            # If battery energy is insufficient, use remaining battery energy and cover the rest with fuel
            if remaining_battery_energy > 0:
                partial_duration = remaining_battery_energy / power
                battery_usage_bins[bin_index] += remaining_battery_energy
                energy_required -= remaining_battery_energy
                remaining_battery_energy = 0

            # Use fuel for the remaining energy required
            fuel_usage_bins[bin_index] += energy_required
            remaining_fuel_energy -= energy_required

# Calculate total energy used
total_energy_used = battery_usage_bins + fuel_usage_bins

# Calculate percentage of total energy used, handling division by zero
battery_percentage = np.where(total_energy_used != 0, (battery_usage_bins / total_energy_used) * 100, 0)
fuel_percentage = np.where(total_energy_used != 0, (fuel_usage_bins / total_energy_used) * 100, 0)

# Create the bar plot
bar_width = 0.45
x = np.arange(len(distance_bins) - 1)

fig, ax = plt.subplots(figsize=(20, 6))

# Plot battery and fuel usage as percentage of total energy available
bar1 = ax.bar(x - bar_width / 2, battery_percentage, bar_width, label="Battery", color="green")
bar2 = ax.bar(x + bar_width / 2, fuel_percentage, bar_width, label="Fuel", color="blue")

# Add labels and title with custom fonts
font_dict = {"fontsize": 20}
ax.set_xlabel("Mission Distance [NM]", fontdict=font_dict)
ax.set_ylabel("Energy Provided [%]", fontdict=font_dict)
ax.set_xticks(x)
ax.set_xticklabels([f"{int(distance_bins[i])}-{int(distance_bins[i+1])}" for i in range(len(distance_bins) - 1)])

# Customize tick parameters
ax.tick_params(axis="both", which="major", labelsize=15)
ax.tick_params(axis="both", which="minor", labelsize=15)


# Add annotations on top of each bar for energy usage phases
def get_phase_annotation(energy_type, battery_energy, fuel_energy, mission):
    if energy_type == "Battery":
        if fuel_energy > 0:
            battery_phases = [
                phase for phase, (seg_phase, _, _, _) in zip(phases, mission) if seg_phase in ("CR", "LD")
            ]
        else:
            battery_phases = [phase for phase, (seg_phase, _, _, _) in zip(phases, mission)]
        return ", ".join(battery_phases)

    elif energy_type == "Fuel" and fuel_energy > 0:
        fuel_phases = ["TO"]  # Always use fuel for takeoff if needed
        if battery_energy <= 7680:
            fuel_phases.append("CR")  # Use fuel for cruise if battery is depleted
        return ", ".join(fuel_phases)

    return ""


for bars, energy_type in zip([bar1, bar2], ["Battery", "Fuel"]):
    for i, bar in enumerate(bars):
        height = bar.get_height()
        battery_energy = battery_usage_bins[i]
        fuel_energy = fuel_usage_bins[i]

        # Get the phase annotation based on the energy type and usage
        phase_annotation = get_phase_annotation(energy_type, battery_energy, fuel_energy, missions[i])

        if phase_annotation:
            ax.annotate(
                phase_annotation,  # Annotation for energy usage during different phases
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),  # 3 points vertical offset
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=11,  # Font size for annotations
                fontweight="bold",  # Font weight for annotations
            )

# Create custom legend entries
legend_elements = [
    Patch(facecolor="green", edgecolor="black", label="Battery"),  # Example for 'Battery'
    Patch(facecolor="blue", edgecolor="black", label="Fuel"),  # Example for 'Fuel'
    Rectangle(
        (0, 0), 1, 1, fc="w", fill=False, edgecolor="none", linewidth=0, label="TO: Take-Off"
    ),  # Example for 'TO'
    Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor="none", linewidth=0, label="LD: Landing"),  # Example for 'LD'
    Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor="none", linewidth=0, label="CR: Cruise"),  # Example for 'CR'
]

# Add custom legend with the specified elements
ax.legend(handles=legend_elements, loc="upper center", prop={"size": 15})

# Save the figure
fig.savefig(
    os.path.join(os.path.dirname(__file__), "..", "..", "Figures", "energy-usage-plot.pdf"),
    bbox_inches="tight",
    dpi=200,
)

# Display the plot
plt.show()
