import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression  # type: ignore[import-untyped]
from sklearn.preprocessing import PolynomialFeatures  # type: ignore[import-untyped]


"MTOW vs Payload"


def plot_MTOW_vs_Payload(combined=True, show=None, save=None, ax=None):
    # single propeller
    MTOW = np.array([6876, 15495, 23814, 20185, 1600, 2678, 2000, 600, 884]) * 9.81  # N
    Payload = np.array([1000, 1800, 2300, 2041, 160, 453.6, 160, 213, 160]) * 9.81  # N

    # plot data points
    ax.scatter(Payload, MTOW, label="Data points")
    x_data = np.linspace(min(Payload), max(Payload), 10000000)

    # show polynomial regression line
    poly_reg = PolynomialFeatures(degree=2)

    # fit_transform() takes the input and returns the transformed input
    payload_poly = poly_reg.fit_transform(Payload.reshape(-1, 1))
    x_data_poly = poly_reg.fit_transform(x_data.reshape(-1, 1))

    # fit the transformed input to the linear regression model
    regressor_poly = LinearRegression()
    regressor_poly.fit(payload_poly, MTOW)
    ax.plot(
        x_data,
        regressor_poly.predict(x_data_poly),
        color="green",
        markersize=1,
        linewidth=3,
        alpha=0.7,
        label="$Polynomial; R^2: 0.9967$",
    )

    # show linear regression line
    regressor = LinearRegression()
    regressor.fit(Payload.reshape(-1, 1), MTOW)
    ax.plot(
        x_data,
        regressor.predict(x_data.reshape(-1, 1)),
        linestyle="dotted",
        markersize=1,
        color="orange",
        linewidth=3,
        alpha=0.9,
        label="$Linear; R^2: 0.978$",
    )

    # show exponential regression line
    regressor_exp = LinearRegression()
    regressor_exp.fit(Payload.reshape(-1, 1), np.log(MTOW))
    ax.plot(
        x_data,
        np.exp(regressor_exp.predict(x_data.reshape(-1, 1))),
        color="blue",
        alpha=0.7,
        linewidth=2.5,
        label="$Exponential; R^2: 0.9637$",
    )

    # show logarithmic regression line
    regressor_log = LinearRegression()
    regressor_log.fit(np.log(Payload).reshape(-1, 1), MTOW)
    ax.plot(
        x_data,
        regressor_log.predict(np.log(x_data).reshape(-1, 1)),
        color="red",
        alpha=0.7,
        linewidth=2.5,
        label="$Logarithmic; R^2: 0.8502$",
    )

    # adding the grid
    ax.grid()

    # adding labels
    ax.set_ylabel("MTOW[N]")
    ax.set_xlabel("Payload[N]")
    legend = ax.legend(loc="upper left")
    # setting the legend's text color to red
    legend.get_texts()[4].set_color("red")

    if show:
        ax.show()

    # save the images based on the case
    if save:
        plt.savefig("MTOW_vs_Payload_FlyingWing.svg")


"MTOW vs OEW"


def plot_MTOW_vs_OEW(combined=True, show=None, save=None, ax=None):
    # single propeller
    MTOW = np.array([6876, 15495, 23814, 20185, 1600, 2678, 2000, 600, 884]) * 9.81  # N
    OEW = np.array([4844, 8918, 13381, 13080, 1360, 1740, 1550, 275, 481]) * 9.81  # N

    # plot data points
    ax.scatter(MTOW, OEW, label="Data points")
    x_data = np.linspace(min(MTOW), max(MTOW), 10000000)

    # show polynomial regression line
    poly_reg = PolynomialFeatures(degree=2)

    # fit_transform() takes the input and returns the transformed input
    mtow_poly = poly_reg.fit_transform(MTOW.reshape(-1, 1))
    x_data_poly = poly_reg.fit_transform(x_data.reshape(-1, 1))

    # fit the transformed input to the linear regression model
    regressor_poly = LinearRegression()
    regressor_poly.fit(mtow_poly, OEW)
    ax.plot(
        x_data,
        regressor_poly.predict(x_data_poly),
        color="green",
        markersize=1,
        linewidth=3,
        alpha=0.7,
        label="$Polynomial; R^2:0.9917$",
    )

    # show linear regression line
    regressor = LinearRegression()
    regressor.fit(MTOW.reshape(-1, 1), OEW)
    ax.plot(
        x_data,
        regressor.predict(x_data.reshape(-1, 1)),
        linestyle="dotted",
        markersize=1,
        color="orange",
        linewidth=3,
        alpha=0.9,
        label="$Linear; R^2: 0.9895$",
    )

    # show exponential regression line
    regressor_exp = LinearRegression()
    regressor_exp.fit(MTOW.reshape(-1, 1), np.log(OEW))
    ax.plot(
        x_data,
        np.exp(regressor_exp.predict(x_data.reshape(-1, 1))),
        color="blue",
        alpha=0.7,
        linewidth=2.5,
        label="$Exponential; R^2: 0.8524$",
    )

    # show logarithmic regression line
    regressor_log = LinearRegression()
    regressor_log.fit(np.log(MTOW).reshape(-1, 1), OEW)
    ax.plot(
        x_data,
        regressor_log.predict(np.log(x_data).reshape(-1, 1)),
        color="red",
        alpha=0.7,
        linewidth=2.5,
        label="$Logarithmic; R^2: 0.8969$",
    )

    # adding the grid
    ax.grid()

    # adding labels
    ax.set_ylabel("OEW'[N]")
    ax.set_xlabel("MTOW[N]")
    legend = ax.legend(loc="upper left")

    # setting the legend's text color to red
    legend.get_texts()[3].set_color("red")
    legend.get_texts()[4].set_color("red")

    if show:
        ax.show()

    # save the images based on the case
    if save:
        plt.savefig("MTOW_vs_OEW_FlyingWing.svg")


if __name__ == "__main__":
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))

    plot_MTOW_vs_Payload(show=False, combined=True, ax=axes[0])
    plot_MTOW_vs_OEW(show=False, combined=True, ax=axes[1])

    plt.tight_layout()
    plt.savefig("FlyingWing.pdf")
    plt.show()
