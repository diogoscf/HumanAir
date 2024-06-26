import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression  # type: ignore[import-untyped]

# from sklearn.preprocessing import PolynomialFeatures


"MTOW vs Payload"


def plot_Payload_MTOW_vs_years(combined=True, show=None, save=None, ax=None):
    # single propeller
    MTOW = np.array(
        [
            0.047467854,
            0.145433392,
            0.116166505,
            0.118872597,
            0.096581843,
            0.077345429,
            0.101114689,
            0.1,
            0.169380134,
            0.08,
            0.355,
            0.180995475,
        ]
    )  # fractions
    years = np.array([1946, 1944, 1947, 1970, 1983, 1947, 2011, 1935, 2003, 1944, 2018, 1966])  # years

    # plot data points
    ax.scatter(years, MTOW, label="Data points")
    x_data = np.linspace(min(years), max(years), 10000000)

    # show linear regression line
    regressor = LinearRegression()
    regressor.fit(years.reshape(-1, 1), MTOW)
    ax.plot(
        x_data,
        regressor.predict(x_data.reshape(-1, 1)),
        linestyle="dotted",
        markersize=1,
        color="orange",
        linewidth=3,
        alpha=0.9,
        label="$Linear; R^2: 0.378$",
    )

    # adding the grid
    ax.grid()

    # adding labels
    ax.set_ylabel("Payload/MTOW")
    ax.set_xlabel("Year")
    ax.legend(loc="upper left")

    if show:
        ax.show()

    # save the images based on the case
    if save:
        plt.savefig("MTOW_vs_Payload_FlyingWing.svg")


"MTOW vs Payload"


def plot_Payload_OEW_vs_years(combined=True, show=None, save=None, ax=None):
    # single propeller
    OEW = np.array(
        [
            0.109019551,
            0.206440958,
            0.201838977,
            0.253054393,
            0.171885509,
            0.169603909,
            0.156039755,
            0.117647059,
            0.260689655,
            0.103225806,
            0.774545455,
            0.332640333,
        ]
    )  # fractions
    years = np.array([1946, 1944, 1947, 1970, 1983, 1947, 2011, 1935, 2003, 1944, 2018, 1966])  # years

    # plot data points
    ax.scatter(years, OEW, label="Data points")
    x_data = np.linspace(min(years), max(years), 10000000)

    # show linear regression line
    regressor = LinearRegression()
    regressor.fit(years.reshape(-1, 1), OEW)
    ax.plot(
        x_data,
        regressor.predict(x_data.reshape(-1, 1)),
        linestyle="dotted",
        markersize=1,
        color="orange",
        linewidth=3,
        alpha=0.9,
        label="$Linear; R^2: 0.376$",
    )

    # adding the grid
    ax.grid()

    # adding labels
    ax.set_ylabel("Payload/OEW")
    ax.set_xlabel("Year")
    ax.legend(loc="upper left")

    if show:
        ax.show()

    # save the images based on the case
    if save:
        plt.savefig("MTOW_vs_Payload_FlyingWing.svg")


if __name__ == "__main__":
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))

    plot_Payload_MTOW_vs_years(show=False, combined=True, ax=axes[0])
    plot_Payload_OEW_vs_years(show=False, combined=True, ax=axes[1])

    plt.tight_layout()
    plt.savefig("FlyingWing_vs_Years.pdf")
    plt.show()
