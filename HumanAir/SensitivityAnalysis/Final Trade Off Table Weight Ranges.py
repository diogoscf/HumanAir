import numpy as np
import itertools
import matplotlib.pyplot as plt

plt.rcParams["axes.labelsize"] = 24
plt.rcParams["axes.titlesize"] = 24
plt.rcParams["xtick.labelsize"] = 18
plt.rcParams["ytick.labelsize"] = 18
plt.rcParams["legend.fontsize"] = 24
# plt.rcParams["svg.fonttype"] = "none"
# plt.rcParams["figure.labelsize"] = 16
# plt.rcParams["font.family"] = "serif"
# plt.rcParams["font.serif"] = "CMU Serif"

# Define the trade off table
criteria = ["Ground Clearance", "Emission reduction", "Development Risk", "Cost", "Operability"]
initial_weights = np.array([0.1, 0.3, 0.2, 0.3, 0.1])
concepts = ["Flying Wing", "Canard", "Conventional"]

# Define importance
initial_importance = np.array([[3, 1, 1, 2, 1], [3, 4, 4, 4, 3], [4, 3, 5, 5, 5]])


# Function to calculate final scores
def calculate_scores(weights, importance):
    return np.sum(importance * weights, axis=1)


# Generate weight ranges within ±25% of the initial weights
weight_ranges_1 = [np.linspace(w * 0.75, w * 1.25, 10) for w in initial_weights]
weight_ranges_2 = [np.linspace(w * 0.5, w * 1.5, 10) for w in initial_weights]
weight_ranges_3 = [np.linspace(w * 0.25, w * 1.75, 10) for w in initial_weights]
weight_ranges_4 = [np.linspace(w * 0, w * 2, 10) for w in initial_weights]

# Generate all combinations of weights that sum to 1
weight_combinations_1 = []
weight_combinations_2 = []
weight_combinations_3 = []
weight_combinations_4 = []
for w in itertools.product(*weight_ranges_1):
    w = np.array(w)  # type: ignore[assignment]
    if np.isclose(np.sum(w), 1):
        weight_combinations_1.append(w)

for w in itertools.product(*weight_ranges_2):
    w = np.array(w)  # type: ignore[assignment]
    if np.isclose(np.sum(w), 1):
        weight_combinations_2.append(w)

for w in itertools.product(*weight_ranges_3):
    w = np.array(w)  # type: ignore[assignment]
    if np.isclose(np.sum(w), 1):
        weight_combinations_3.append(w)

for w in itertools.product(*weight_ranges_4):
    w = np.array(w)  # type: ignore[assignment]
    if np.isclose(np.sum(w), 1):
        weight_combinations_4.append(w)

# Prepare to store results
score_freq_1 = {concept: 0 for concept in concepts}
score_freq_2 = {concept: 0 for concept in concepts}
score_freq_3 = {concept: 0 for concept in concepts}
score_freq_4 = {concept: 0 for concept in concepts}

# Evaluate all combinations of weights and initial importance values for each range
for weights in weight_combinations_1:
    final_scores = calculate_scores(weights, initial_importance)
    winning_concept_index = np.argmax(final_scores)
    winning_concept = concepts[winning_concept_index]
    score_freq_1[winning_concept] += 1

for weights in weight_combinations_2:
    final_scores = calculate_scores(weights, initial_importance)
    winning_concept_index = np.argmax(final_scores)
    winning_concept = concepts[winning_concept_index]
    score_freq_2[winning_concept] += 1

for weights in weight_combinations_3:
    final_scores = calculate_scores(weights, initial_importance)
    winning_concept_index = np.argmax(final_scores)
    winning_concept = concepts[winning_concept_index]
    score_freq_3[winning_concept] += 1

for weights in weight_combinations_4:
    final_scores = calculate_scores(weights, initial_importance)
    winning_concept_index = np.argmax(final_scores)
    winning_concept = concepts[winning_concept_index]
    score_freq_4[winning_concept] += 1

# Plot the results
# Plot the results
plt.figure(figsize=(12, 6))
x_pos = np.arange(len(concepts))
bar_width = 0.2

total_cases = len(weight_combinations_1)  # Assuming weight_combinations_1 has the total number of cases

# Plot frequency for each range
plt.bar(
    x_pos - bar_width * 1.5,
    [score_freq_1[concept] * 100 / total_cases for concept in concepts],
    bar_width,
    label="±25%",
    color="blue",
)
plt.bar(
    x_pos - bar_width / 2,
    [score_freq_2[concept] * 100 / total_cases for concept in concepts],
    bar_width,
    label="±50%",
    color="orange",
)
plt.bar(
    x_pos + bar_width / 2,
    [score_freq_3[concept] * 100 / total_cases for concept in concepts],
    bar_width,
    label="±75%",
    color="red",
)
plt.bar(
    x_pos + bar_width * 1.5,
    [score_freq_4[concept] * 100 / total_cases for concept in concepts],
    bar_width,
    label="±100%",
    color="green",
)

# Add labels and title
# plt.xlabel('Concepts')
plt.ylabel("Win Frequency [%]")
plt.xticks(x_pos, concepts, fontsize=24)
plt.legend(loc="upper left")

# Add text labels for each bar
FONTSIZE_TXT = 16
for i in range(len(concepts)):
    plt.text(
        x_pos[i] - bar_width * 1.5,
        score_freq_1[concepts[i]] * 100 / total_cases,
        f"{(score_freq_1[concepts[i]] / total_cases) * 100:.0f}%",
        ha="center",
        va="bottom",
        fontsize=FONTSIZE_TXT,
    )
    plt.text(
        x_pos[i] - bar_width / 2,
        score_freq_2[concepts[i]] * 100 / total_cases,
        f"{(score_freq_2[concepts[i]] / total_cases) * 100:.0f}%",
        ha="center",
        va="bottom",
        fontsize=FONTSIZE_TXT,
    )
    plt.text(
        x_pos[i] + bar_width / 2,
        score_freq_3[concepts[i]] * 100 / total_cases,
        f"{(score_freq_3[concepts[i]] / total_cases) * 100:.0f}%",
        ha="center",
        va="bottom",
        fontsize=FONTSIZE_TXT,
    )
    plt.text(
        x_pos[i] + bar_width * 1.5,
        score_freq_4[concepts[i]] * 100 / total_cases,
        f"{(score_freq_4[concepts[i]] / total_cases) * 100:.0f}%",
        ha="center",
        va="bottom",
        fontsize=FONTSIZE_TXT,
    )

plt.tight_layout()
plt.savefig("Figures/sensitivity_analysis.pdf")
plt.show()
