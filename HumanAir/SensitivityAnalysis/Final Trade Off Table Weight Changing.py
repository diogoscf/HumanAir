import numpy as np
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


# Perform sensitivity analysis
weight_ranges = np.linspace(0, 1, 22)  # From 0 to 1 in steps of 0.05

# To store results for each criteria being varied
sensitivity_results = {concept: {criteria[i]: [] for i in range(len(criteria))} for concept in concepts}  # type: ignore[var-annotated]   # noqa: E501


# Iterate through each criteria weight
for i in range(len(initial_weights)):
    for w in weight_ranges:
        weights = initial_weights.copy()
        weights[i] = w
        weights /= np.sum(weights)  # Normalize to ensure weights sum to 1
        final_scores = np.sum(initial_importance * weights, axis=1)
        for j, concept in enumerate(concepts):
            sensitivity_results[concept][criteria[i]].append(final_scores[j])

# Plot the sensitivity results
fig, axes = plt.subplots(len(initial_weights), 1, figsize=(10, 15))
last = len(criteria) - 1
for i, criterion in enumerate(criteria):
    ax = axes[i]
    colours = "brg"
    for j, concept in enumerate(concepts):
        label = concept if i == last else None
        ax.plot(weight_ranges, sensitivity_results[concept][criterion], label=label, color=colours[j])
    if i == last:
        ax.set_xlabel("Weight Value")
        # ax.legend(loc = 'lower right', bbox_to_anchor=(1.0, 0.05))
    ax.set_ylabel("Final Score")
    ax.set_title(f"{criterion}")
    ax.grid()

    ax.set_yticks(np.arange(0, 10, 1))
    ax.set_yticks(np.arange(0, 10, 0.5), minor=True)
    ax.set_ylim(1.0, 5.0)

    ax.set_xticks(np.arange(0, 1.01, 0.1))
    ax.set_xticks(np.arange(-0.05, 1.06, 0.05), minor=True)
    # ax.set_xlim(1.0, 5.0)

    # ax.tick_params("both", length=10, width=1, which="major")
    # ax.tick_params("both", length=5, width=1, which="minor")


handles, labels = axes[-1].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=len(criteria), frameon=False)
plt.tight_layout()
plt.subplots_adjust(bottom=0.05, top=0.93, wspace=0, hspace=0.38)
plt.savefig("Figures/Sensitivity Analysis of Weight.pdf")
# plt.show()


# Perform sensitivity analysis for importance
importance_ranges = range(1, 6)  # From 1 to 5
sensitivity_results_importance = {concept: {criteria[i]: [] for i in range(len(criteria))} for concept in concepts}  # type: ignore[var-annotated]  # noqa: E501

for i in range(initial_importance.shape[1]):
    for imp in importance_ranges:
        importance = initial_importance.copy()
        importance[:, i] = imp
        final_scores = np.sum(importance * initial_weights, axis=1)
        for j, concept in enumerate(concepts):
            sensitivity_results_importance[concept][criteria[i]].append(final_scores[j])

# Plot the sensitivity results for importance
fig, axes = plt.subplots(len(criteria), 1, figsize=(10, 15))
last = len(criteria) - 1
for i, criterion in enumerate(criteria):
    ax = axes[i]
    colours = "brg"
    for j, concept in enumerate(concepts):
        label = concept if i == last else None
        ax.plot(importance_ranges, sensitivity_results_importance[concept][criterion], label=label, color=colours[j])
    ax.set_ylabel("Final Score")
    if i == last:
        ax.set_xlabel("Criterion Score")

    ax.set_title(f"{criterion}")
    # ax.legend(loc = 'lower right')
    ax.grid()

    ax.set_yticks(np.arange(0, 10, 1))
    ax.set_yticks(np.arange(0, 10, 0.5), minor=True)
    ax.set_ylim(1.0, 5.0)

    ax.set_xticks(np.arange(1, 5.01, 1))
    ax.set_xticks(np.arange(1, 5.01, 0.5), minor=True)
    # ax.set_xlim(1.0, 5.0)

    # ax.tick_params("both", length=10, width=1, which="major")
    # ax.tick_params("both", length=5, width=1, which="minor")

handles, labels = axes[-1].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=len(criteria), frameon=False)
plt.tight_layout()
plt.subplots_adjust(bottom=0.05, top=0.93, wspace=0, hspace=0.38)
plt.savefig("Figures/Sensitivity Analysis of Importance.pdf")
# plt.show()
