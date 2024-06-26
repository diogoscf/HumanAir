import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression  # type: ignore[import-untyped]
from sklearn.preprocessing import PolynomialFeatures  # type: ignore[import-untyped]


"MTOW vs Payload"


def plot_MTOW_vs_Payload(combined=True, show=None, save=None, ax=None):
    # single propeller
    MTOW_single_propeller = (
        np.array([680.388, 997.903, 1065.942, 1202.02, 1202.02, 952.544, 1315.418, 1315.418, 1339.912]) * 9.81
    )  # N
    Payload_single_propeller = (
        np.array([180.529, 318.421, 326.1320, 332.9360, 324.3180, 145.6030, 494.4150, 383.2850, 335.6580]) * 9.81
    )  # N

    # double propeller
    MTOW_double_propeller = (
        np.array(
            [
                3342.068582,
                4399.845989,
                3243.185446,
                1814.36948,
                2190.851147,
                3810.175908,
                3077.62423,
                4501.904272,
                4628.910136,
            ]
        )
        * 9.81
    )  # N
    Payload_double_propeller = (
        np.array(
            [
                583.7733802,
                668.5951534,
                594.6595971,
                332.4832072,
                537.9605508,
                861.825503,
                544.310844,
                725.747792,
                889.0410452,
            ]
        )
        * 9.81
    )  # N

    # combine data points for plotting
    if combined:
        MTOW_ = np.concatenate((MTOW_single_propeller, MTOW_double_propeller))
        Payload_ = np.concatenate((Payload_single_propeller, Payload_double_propeller))
        MTOW_ = [MTOW_]
        Payload_ = [Payload_]

    else:
        MTOW_ = [MTOW_single_propeller, MTOW_double_propeller]
        Payload_ = [Payload_single_propeller, Payload_double_propeller]

    # keep track of the combined parameter in order to know how to save the images
    step = 0

    for MTOW in MTOW_:
        for Payload in Payload_:
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
                label="$Polynomial; R^2: 0.8501$",
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
                label="$Linear; R^2: 0.8497$",
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
                label="$Exponential; R^2: 0.772$",
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
                label="$Logarithmic; R^2: 0.75$",
            )

            # adding the grid
            ax.grid()

            # adding labels
            ax.set_ylabel("MTOW[N]")
            ax.set_xlabel("Payload[N]")
            legend = ax.legend(loc="upper left")
            # setting the legend's text color to red
            legend.get_texts()[4].set_color("red")
            legend.get_texts()[3].set_color("red")

            if show:
                ax.show()

            # save the images based on the case
            if save:
                step += 1
                if step == 1 and combined:
                    plt.savefig("Combined_MTOW_vs_Payload.svg")

                elif step == 1 and not combined:
                    plt.savefig("Single_MTOW_vs_Payload.svg")

                elif step == 4 and not combined:
                    plt.savefig("Double_MTOW_vs_Payload.svg")


"MTOW vs OEW"


def plot_MTOW_vs_OEW(combined=True, show=None, save=None, ax=None):
    # single propeller
    MTOW_single_propeller = (
        np.array([680.388, 997.903, 1065.942, 1202.02, 1202.02, 952.544, 1315.418, 1315.418]) * 9.81
    )  # N
    OEW_single_propeller = (
        np.array(
            [210.0132673, 304.8140726, 306.6284421, 338.379908, 346.9981631, 335.2047614, 406.8723559, 389.6358458]
        )
        * 9.81
    )  # N

    # double propeller
    MTOW_double_propeller = (
        np.array(
            [
                3342.068582,
                4399.845989,
                3243.185446,
                1814.36948,
                2190.851147,
                3810.175908,
                3077.62423,
                4501.904272,
                4628.910136,
            ]
        )
        * 9.81
    )  # N
    OEW_double_propeller = (
        np.array(
            [
                1098.147128,
                1434.712666,
                994.7280674,
                573.3407557,
                602.8242597,
                1143.052772,
                1122.187523,
                1650.169042,
                1782.164422,
            ]
        )
        * 9.81
    )  # N

    # combine data points for plotting
    if combined:
        MTOW_ = np.concatenate((MTOW_single_propeller, MTOW_double_propeller))
        OEW_ = np.concatenate((OEW_single_propeller, OEW_double_propeller))
        MTOW_ = [MTOW_]
        OEW_ = [OEW_]

    else:
        MTOW_ = [MTOW_single_propeller, MTOW_double_propeller]
        OEW_ = [OEW_single_propeller, OEW_double_propeller]

    # keep track of the combined parameter in order to know how to save the images
    step = 0

    for MTOW in MTOW_:
        for OEW in OEW_:
            step += 1
            if step == 1 or step == 4:
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
                    label="$Polynomial; R^2:0.9814$",
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
                    label="$Linear; R^2: 0.9735$",
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
                    label="$Exponential; R^2: 0.9676$",
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
                    label="$Logarithmic; R^2: 0.9039$",
                )

                # adding the grid
                ax.grid()

                # adding labels
                ax.set_ylabel("OEW'[N]")
                ax.set_xlabel("MTOW[N]")
                legend = ax.legend(loc="upper left")

                # setting the legend's text color to red
                legend.get_texts()[4].set_color("red")

                if show:
                    ax.show()

                # save the images based on the case
                if save:
                    if step == 1 and combined:
                        plt.savefig("Combined_MTOW_vs_OEW.svg")

                    elif step == 1 and not combined:
                        plt.savefig("Single_MTOW_vs_OEW.svg")

                    elif step == 4 and not combined:
                        plt.savefig("Double_MTOW_vs_OEW.svg")


"MTOW vs WW"


def plot_MTOW_vs_WW(combined=True, show=None, save=None, ax=None):
    # single propeller
    MTOW_single_propeller = np.array([680.388, 997.903, 1065.942, 1202.02, 1202.02, 952.544, 1315.418]) * 9.81  # N
    WW_single_propeller = (
        np.array([97.97595192, 102.5118756, 102.965468, 106.594207, 106.594207, 107.9549841, 118.3876086]) * 9.81
    )  # N

    # double propeller
    MTOW_double_propeller = (
        np.array(
            [
                3342.068582,
                4399.845989,
                3243.185446,
                1814.36948,
                2190.851147,
                3810.175908,
                3077.62423,
                4501.904272,
                4628.910136,
            ]
        )
        * 9.81
    )  # N
    WW_double_propeller = (
        np.array(
            [
                303.9068879,
                396.4397314,
                297.5565947,
                207.7453055,
                205.4773436,
                390.0894382,
                289.3919321,
                395.986139,
                454.0459624,
            ]
        )
        * 9.81
    )  # N

    # combine data points for plotting
    if combined:
        MTOW_ = np.concatenate((MTOW_single_propeller, MTOW_double_propeller))
        WW_ = np.concatenate((WW_single_propeller, WW_double_propeller))
        MTOW_ = [MTOW_]
        WW_ = [WW_]

    else:
        MTOW_ = [MTOW_single_propeller, MTOW_double_propeller]
        WW_ = [WW_single_propeller, WW_double_propeller]

    # keep track of the combined parameter in order to know how to save the images
    step = 0

    for MTOW in MTOW_:
        for WW in WW_:
            step += 1
            if step == 1 or step == 4:
                # plot data points
                ax.scatter(MTOW, WW, label="Data points")
                x_data = np.linspace(min(MTOW), max(MTOW), 10000000)

                # show polynomial regression line
                poly_reg = PolynomialFeatures(degree=2)

                # fit_transform() takes the input and returns the transformed input
                mtow_poly = poly_reg.fit_transform(MTOW.reshape(-1, 1))
                x_data_poly = poly_reg.fit_transform(x_data.reshape(-1, 1))

                # fit the transformed input to the linear regression model
                regressor_poly = LinearRegression()
                regressor_poly.fit(mtow_poly, WW)
                ax.plot(
                    x_data,
                    regressor_poly.predict(x_data_poly),
                    color="green",
                    markersize=1,
                    linewidth=3,
                    alpha=0.7,
                    label="$Polynomial; R^2:0.981$",
                )

                # show linear regression line
                regressor = LinearRegression()
                regressor.fit(MTOW.reshape(-1, 1), WW)
                ax.plot(
                    x_data,
                    regressor.predict(x_data.reshape(-1, 1)),
                    linestyle="dotted",
                    markersize=1,
                    color="orange",
                    linewidth=3,
                    alpha=0.9,
                    label="$Linear; R^2: 0.9809$",
                )

                # show exponential regression line
                regressor_exp = LinearRegression()
                regressor_exp.fit(MTOW.reshape(-1, 1), np.log(WW))
                ax.plot(
                    x_data,
                    np.exp(regressor_exp.predict(x_data.reshape(-1, 1))),
                    color="blue",
                    alpha=0.7,
                    linewidth=2.5,
                    label="$Exponential; R^2: 0.9489$",
                )

                # show logarithmic regression line
                regressor_log = LinearRegression()
                regressor_log.fit(np.log(MTOW).reshape(-1, 1), WW)
                ax.plot(
                    x_data,
                    regressor_log.predict(np.log(x_data).reshape(-1, 1)),
                    color="red",
                    alpha=0.7,
                    linewidth=2.5,
                    label="$Logarithmic; R^2:0.9307$",
                )

                # adding the grid
                ax.grid()

                # adding labels
                ax.set_ylabel("WW[N]")
                ax.set_xlabel("MTOW[N]")
                ax.legend(loc="upper left")

                if show:
                    ax.show()

                # save the images based on the case
                if save:
                    if step == 1 and combined:
                        plt.savefig("Combined_MTOW_vs_WW.svg")

                    elif step == 1 and not combined:
                        plt.savefig("Single_MTOW_vs_WW.svg")

                    elif step == 4 and not combined:
                        plt.savefig("Double_MTOW_vs_WW.svg")


if __name__ == "__main__":
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    plot_MTOW_vs_Payload(show=False, combined=True, ax=axes[0])
    plot_MTOW_vs_OEW(show=False, combined=True, ax=axes[1])
    plot_MTOW_vs_WW(show=False, combined=True, ax=axes[2])

    plt.tight_layout()
    plt.savefig("Combined_MTOW_vs_Payload_OEW_WW.pdf")
    plt.show()
