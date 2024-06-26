import numpy as np
from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression  # type: ignore[import-untyped]
from sklearn.preprocessing import PolynomialFeatures  # type: ignore[import-untyped]

# set the plotting size and style

fig, ax = plt.subplots(1, 3, figsize=(15, 5))


"Conventional Battery Regression"

bat_percentage = np.array([5, 10, 15, 20, 25])  # percentages

P_max_conventional = np.array([148658.911, 185715.7703, 250282.1975, 305776.0373, 514703.2653])  # W


ax[0].scatter(bat_percentage, P_max_conventional, label="Data points")
x_data = np.linspace(min(bat_percentage), max(bat_percentage), 10000000)

# show polynomial regression line
poly_reg = PolynomialFeatures(degree=2)

# fit_transform() takes the input and returns the transformed input
battery_poly = poly_reg.fit_transform(bat_percentage.reshape(-1, 1))
x_data_poly = poly_reg.fit_transform(x_data.reshape(-1, 1))

# fit the transformed input to the linear regression model
regressor_poly = LinearRegression()
regressor_poly.fit(battery_poly, P_max_conventional)
ax[0].plot(
    x_data,
    regressor_poly.predict(x_data_poly),
    color="green",
    markersize=1,
    linewidth=3,
    alpha=0.7,
    label="$Polynomial$",
)
ax[0].legend(loc="upper left")

ax[0].grid()
ax[0].set_title("Conventional")


"Canard Battery Regression"
P_max_canard = np.array([140808.56, 173918.24, 229886.83, 305347.04, 433016.05])  # W


ax[1].scatter(bat_percentage, P_max_canard, label="Data points")
x_data = np.linspace(min(bat_percentage), max(bat_percentage), 10000000)

# show polynomial regression line
poly_reg = PolynomialFeatures(degree=2)

# fit_transform() takes the input and returns the transformed input
battery_poly = poly_reg.fit_transform(bat_percentage.reshape(-1, 1))
x_data_poly = poly_reg.fit_transform(x_data.reshape(-1, 1))

# fit the transformed input to the linear regression model
regressor_poly = LinearRegression()
regressor_poly.fit(battery_poly, P_max_canard)
ax[1].plot(
    x_data,
    regressor_poly.predict(x_data_poly),
    color="green",
    markersize=1,
    linewidth=3,
    alpha=0.7,
    label="$Polynomial$",
)
ax[1].legend(loc="upper left")

ax[1].grid()
ax[1].set_title("Canard")

"Flying Wing Battery Regression"

bat_percentage = np.array([2.5, 5, 7.5, 10, 12.5, 15])  # percentages
P_max_flying_wing = np.array([217187.6014, 253011.4544, 304494.2278, 372678.8043, 493252.7253, 712382.756])  # W


ax[2].scatter(bat_percentage, P_max_flying_wing, label="Data points")
x_data = np.linspace(min(bat_percentage), max(bat_percentage), 10000000)

# show polynomial regression line
poly_reg = PolynomialFeatures(degree=2)

# fit_transform() takes the input and returns the transformed input
battery_poly = poly_reg.fit_transform(bat_percentage.reshape(-1, 1))
x_data_poly = poly_reg.fit_transform(x_data.reshape(-1, 1))

# fit the transformed input to the linear regression model
regressor_poly = LinearRegression()
regressor_poly.fit(battery_poly, P_max_flying_wing)
ax[2].plot(
    x_data,
    regressor_poly.predict(x_data_poly),
    color="green",
    markersize=1,
    linewidth=3,
    alpha=0.7,
    label="$Polynomial$",
)
ax[2].legend(loc="upper left")

ax[2].grid()
ax[2].set_title("Flying Wing")


# set common labels for the whole big figure

fig.text(0.5, 0.02, "Battery Percentage", ha="center", fontsize=12)
fig.text(0.05, 0.5, "$P_{max} [W]$", va="center", rotation="vertical", fontsize=12)


plt.savefig("Battery_Regression.pdf")
plt.show()
