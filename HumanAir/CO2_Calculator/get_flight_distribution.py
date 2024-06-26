import numpy as np
import scipy.stats as stats  # type: ignore[import-untyped]
import matplotlib.pyplot as plt
import os
import pandas as pd

# import json
import pickle


def fit_log_normal(mean, median):
    mu = np.log(median)
    sigma = np.sqrt(2 * np.log(mean / median))

    return stats.lognorm(s=sigma, scale=np.exp(mu))


if __name__ == "__main__":  # pragma: no cover
    with open(os.path.join(os.path.dirname(__file__), "maf_flights_before_refuelling.csv"), "r") as f:
        flight_data = pd.read_csv(f)

    x = np.linspace(0, 800, 1000)
    colours = "brgcmy"
    combined = np.array([])
    plt.rcParams.update({"font.size": 30})
    plt.rcParams.update({"axes.labelsize": 40})
    fig, ax = plt.subplots()
    fig.tight_layout()
    ax.axhline(linewidth=2, color="k")
    ax.axvline(linewidth=2, color="k")

    for i, row in flight_data.iterrows():
        dist = fit_log_normal(row["distance_before_refuelling_avg"], row["distance_before_refuelling_median"])
        print(
            i, dist.mean(), dist.median(), row["distance_before_refuelling_min"], row["distance_before_refuelling_max"]
        )
        simulated_lengths = dist.rvs(size=10000)
        ax.plot(x, dist.pdf(x), colours[i] + "-", lw=2, label=row["aircraft"])  # type: ignore[index]
        combined = np.concatenate((combined, simulated_lengths))

    combined = combined[np.where(combined <= 600)]
    combined = combined[np.where(combined >= 0)]
    shape, loc, scale = stats.lognorm.fit(combined, floc=0)
    overall_dist = stats.lognorm(s=shape, loc=loc, scale=scale)

    ax.plot(x, overall_dist.pdf(x), "k--", lw=2, label="Overall")
    print("Overall", overall_dist.mean(), overall_dist.median())

    # Plot simulated data
    # plt.hist(simulated_lengths, bins=20, density=True, alpha=0.6, color='g', edgecolor='black')

    # x = np.linspace(0, 1, 1000)
    # plt.plot(x, dist.pdf(x), 'r-', lw=2)
    # ax.set_title("Simulated Flight Lengths Distribution")
    ax.set_xlabel("Flight Length [NM]")
    ax.set_ylabel("Probability Density")
    ax.grid()
    ax.legend()

    fig.set_size_inches(16, 9)
    fig.tight_layout()

    fig.savefig(
        os.path.join(os.path.dirname(__file__), "..", "..", "Figures", "flight-length-distribution.pdf"),
        bbox_inches="tight",
        dpi=200,
    )
    plt.show()

    if input("Pickle distribution? (y/N): ").lower() == "y":
        with open(os.path.join(os.path.dirname(__file__), "flightdist.pickle"), "wb+") as f:  # type: ignore[assignment]
            pickle.dump(overall_dist, f)  # type: ignore[arg-type]
