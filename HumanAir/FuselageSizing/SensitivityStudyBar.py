import numpy as np
import itertools
import matplotlib.pyplot as plt

# Define the trade off table
criteria = ["Structure Weight", "Subsystem Integration", "Ease of Manufacturing", "Manufacturing Cost", "Aerodynamics"]
initial_weights = np.array([0.3, 0.1, 0.2, 0.2, 0.2])
concepts = ["Monocoque", "Semi-monocoque", "Truss", "Geodesic"]

# Define importance
initial_importance = np.array([[5, 3, 2, 3, 5], [4, 5, 3, 4, 5], [3, 1, 5, 5, 3], [3, 3, 1, 1, 4]])


# Function to calculate final scores
def calculate_scores(weights, importance):
    return np.sum(importance * weights, axis=1)


# Generate weight ranges within ±25% of the initial weights
weight_ranges_1 = [np.linspace(w * 0.75, w * 1.25, 15) for w in initial_weights]
weight_ranges_2 = [np.linspace(w * 0.5, w * 1.5, 15) for w in initial_weights]
weight_ranges_3 = [np.linspace(w * 0.25, w * 1.75, 15) for w in initial_weights]
weight_ranges_4 = [np.linspace(w * 0, w * 2, 15) for w in initial_weights]

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

# Plot frequency for each range
plt.bar(x_pos - bar_width * 1.5, [score_freq_1[concept] for concept in concepts], bar_width, label="±25%", color="blue")
plt.bar(x_pos - bar_width / 2, [score_freq_2[concept] for concept in concepts], bar_width, label="±50%", color="orange")
plt.bar(x_pos + bar_width / 2, [score_freq_3[concept] for concept in concepts], bar_width, label="±75%", color="red")
plt.bar(
    x_pos + bar_width * 1.5, [score_freq_4[concept] for concept in concepts], bar_width, label="±100%", color="green"
)

# Add labels and title
plt.xlabel("Concepts")
plt.ylabel("Frequency of Winning")
plt.xticks(x_pos, concepts)
plt.legend(loc="upper left")

# Add text labels for each bar
total_cases = len(weight_combinations_1)  # Assuming weight_combinations_1 has the total number of cases
for i in range(len(concepts)):
    plt.text(
        x_pos[i] - bar_width * 1.5,
        score_freq_1[concepts[i]],
        f" ({(score_freq_1[concepts[i]] / total_cases) * 100:.2f}%)",
        ha="center",
        va="bottom",
        fontsize=10,
    )
    plt.text(
        x_pos[i] - bar_width / 2,
        score_freq_2[concepts[i]],
        f" ({(score_freq_2[concepts[i]] / total_cases) * 100:.2f}%)",
        ha="center",
        va="bottom",
        fontsize=10,
    )
    plt.text(
        x_pos[i] + bar_width / 2,
        score_freq_3[concepts[i]],
        f" ({(score_freq_3[concepts[i]] / total_cases) * 100:.2f}%)",
        ha="center",
        va="bottom",
        fontsize=10,
    )
    plt.text(
        x_pos[i] + bar_width * 1.5,
        score_freq_4[concepts[i]],
        f" ({(score_freq_4[concepts[i]] / total_cases) * 100:.2f}%)",
        ha="center",
        va="bottom",
        fontsize=10,
    )

plt.tight_layout()
plt.savefig("sensitivity_analysis.pdf")
plt.show()
