import numpy as np
import matplotlib.pyplot as plt

# from scipy.optimize import curve_fit


# generate flight profile with 3 phases. Ground, flight, and landing

factor = 4

# ground phase
ground_phase_start = 0
ground_phase_end = 20
ground_load = -5
x1 = np.linspace(ground_phase_start, ground_phase_end, (ground_phase_end - ground_phase_start) * factor)
y1 = np.ones_like(x1) * ground_load

# flight phase
flight_phase_start = 24
flight_phase_end = 44
flight_load = 10
x2 = np.linspace(flight_phase_start, flight_phase_end, (flight_phase_end - flight_phase_start) * factor)
y2 = np.ones_like(x2) * flight_load

# landing phase
landing_phase_start = 48
landing_load = -8
x3 = np.array([landing_phase_start])
y3 = np.ones_like(x3) * landing_load

# Second ground phase
ground_phase_start = 49
ground_phase_end = 68
ground_load = -5
x4 = np.linspace(ground_phase_start, ground_phase_end, (ground_phase_end - ground_phase_start) * factor)
y4 = np.ones_like(x4) * ground_load

# between ground and flight phase
transition_phase_start = 20
transition_phase_end = 24
transition_load = np.linspace(-5, 10, (transition_phase_end - transition_phase_start) * factor)
x5 = np.linspace(transition_phase_start, transition_phase_end, (transition_phase_end - transition_phase_start) * factor)
y5 = transition_load

# between flight and landing phase
transition_phase_start = 44
transition_phase_end = 48
transition_load = np.linspace(10, -8, (transition_phase_end - transition_phase_start) * factor)
x6 = np.linspace(transition_phase_start, transition_phase_end, (transition_phase_end - transition_phase_start) * factor)
y6 = transition_load

# between landing and second ground phase
transition_phase_start = 48
transition_phase_end = 49
transition_load = np.linspace(-8, -5, (transition_phase_end - transition_phase_start) * factor)
x7 = np.linspace(transition_phase_start, transition_phase_end, (transition_phase_end - transition_phase_start) * factor)
y7 = transition_load

# Apply noise to all phases deviating from the mean value with a normal distribution
noise = 0.3
y1 += np.random.normal(0, noise, y1.shape)  # type: ignore[arg-type]
y2 += np.random.normal(0, noise, y2.shape)  # type: ignore[arg-type]
y4 += np.random.normal(0, noise, y4.shape)  # type: ignore[arg-type]
y5 += np.random.normal(0, noise, y5.shape)  # type: ignore[arg-type]
y6 += np.random.normal(0, noise, y6.shape)  # type: ignore[arg-type]
y7 += np.random.normal(0, noise, y7.shape)  # type: ignore[arg-type]
xa = np.concatenate((x1, x5, x2, x6, x3, x7, x4))
ya = np.concatenate((y1, y5, y2, y6, y3, y7, y4))


# New curve:
# ground phase
ground_phase_start = 0
ground_phase_end = 20
ground_load = -10
x1 = np.linspace(ground_phase_start, ground_phase_end, (ground_phase_end - ground_phase_start) * factor)
y1 = np.ones_like(x1) * ground_load

# flight phase
flight_phase_start = 24
flight_phase_end = 44
flight_load = 10
x2 = np.linspace(flight_phase_start, flight_phase_end, (flight_phase_end - flight_phase_start) * factor)
y2 = np.ones_like(x2) * flight_load

# landing phase
landing_phase_start = 48
landing_load = -10
x3 = np.array([landing_phase_start])
y3 = np.ones_like(x3) * landing_load

# Second ground phase
ground_phase_start = 49
ground_phase_end = 68
ground_load = -10
x4 = np.linspace(ground_phase_start, ground_phase_end, (ground_phase_end - ground_phase_start) * factor)
y4 = np.ones_like(x4) * ground_load

# between ground and flight phase
transition_phase_start = 20
transition_phase_end = 24
transition_load = np.linspace(-10, 10, (transition_phase_end - transition_phase_start) * factor)
x5 = np.linspace(transition_phase_start, transition_phase_end, (transition_phase_end - transition_phase_start) * factor)
y5 = transition_load

# between flight and landing phase
transition_phase_start = 44
transition_phase_end = 48
transition_load = np.linspace(10, -10, (transition_phase_end - transition_phase_start) * factor)
x6 = np.linspace(transition_phase_start, transition_phase_end, (transition_phase_end - transition_phase_start) * factor)
y6 = transition_load

# between landing and second ground phase
transition_phase_start = 48
transition_phase_end = 49
transition_load = np.linspace(-10, -10, (transition_phase_end - transition_phase_start) * factor)
x7 = np.linspace(transition_phase_start, transition_phase_end, (transition_phase_end - transition_phase_start) * factor)
y7 = transition_load

# Apply noise to all phases deviating from the mean value with a normal distribution
noise = 0.3
y1 += np.random.normal(0, noise, y1.shape)  # type: ignore[arg-type]
y2 += np.random.normal(0, noise, y2.shape)  # type: ignore[arg-type]
y4 += np.random.normal(0, noise, y4.shape)  # type: ignore[arg-type]
y5 += np.random.normal(0, noise, y5.shape)  # type: ignore[arg-type]
y6 += np.random.normal(0, noise, y6.shape)  # type: ignore[arg-type]
y7 += np.random.normal(0, noise, y7.shape)  # type: ignore[arg-type]
xb = np.concatenate((x1, x5, x2, x6, x3, x7, x4))
yb = np.concatenate((y1, y5, y2, y6, y3, y7, y4))

# Sine curve
x = np.linspace(0, 68, 1000)
y = -10 * np.sin(2 * np.pi * x / 44)


# plot the flight profile
fig1 = plt.figure()
plt.plot(xa, ya, label="_nolegend_")
plt.plot(x, y, "r-", label="Cyclic Approximation")
plt.xlabel("Time", fontsize=14)
plt.ylabel("Load", fontsize=14)
plt.xticks([])
# Get current ticks
current_ticks = [-10, 0, 10]
# Define the tick labels
tick_labels = ["-S_max", "0", "S_max"]
# increase the size of the ticks
plt.yticks(fontsize=14)
# Set the new ticks and labels
plt.yticks(current_ticks, tick_labels)
# # Get rid of x and y ticks
# plt.yticks([-10, 0, 10])
# Plot only y = 0 line and indicate
plt.axhline(-10, color="black", lw=1, linestyle="--")
plt.axhline(0, color="black", lw=1)
plt.axhline(10, color="black", lw=1, linestyle="--")
plt.ylim(-12, 14.5)
# put legend center bottom
plt.legend(fontsize=14)
plt.savefig(
    "D:/a_Education/BSc - TUDelft/Aerospace Engineering/BSc year 3/DSE/HumanAir/Fatigue/flight_profile_approx.pdf",
    bbox_inches="tight",
    dpi=200,
)


# Plot a clean flight profile with names of the phases
fig2 = plt.figure()
plt.plot(xa, ya)
plt.xlabel("Time", fontsize=14)
plt.ylabel("Load", fontsize=14)
plt.xticks([])
# Get current ticks
current_ticks = [-8, 0, 10]
# Define the tick labels
tick_labels = ["S_min", "0", "S_max"]
# increase the size of the ticks
plt.yticks(fontsize=14)
# Set the new ticks and labels
plt.yticks(current_ticks, tick_labels)
# # Get rid of x and y ticks
# plt.yticks([-10, 0, 10])
# Plot only y = 0 line and indicate
plt.axhline(-8, color="black", lw=1, linestyle="--")
plt.axhline(0, color="black", lw=1)
plt.axhline(10, color="black", lw=1, linestyle="--")
plt.text(-1, -4.1, "Ground Phase", fontsize=14)
plt.text(25, 11, "Flight Phase", fontsize=14)  # plt.text(25, 8, 'Flight Phase', fontsize=14)
plt.text(35.5, -7.5, "Landing", fontsize=14)  # plt.text(49, -7.7, 'Landing')
plt.text(49, -4.1, "Ground Phase", fontsize=14)
plt.ylim(-12, 14.5)
plt.savefig(
    "D:/a_Education/BSc - TUDelft/Aerospace Engineering/BSc year 3/DSE/HumanAir/Fatigue/flight_profile.pdf",
    bbox_inches="tight",
    dpi=200,
)
plt.show()
