import numpy as np
from HumanAir.FlightPerformance import aircraft

import matplotlib.pyplot as plt


def plot_variation():

    fig = plt.figure(figsize=(8, 5))

    for i in [1, 2, 3, 4]:
        if i == 2:
            # ***** design *****
            elevation = 750
            temp_offset = 18
            slope = 0
            surface = "grass"
            label = "750m ISA+18"
        if i == 1:
            # **** sea level *****
            elevation = 0
            temp_offset = 0
            slope = 0
            surface = "grass"
            label = "0m ISA+0"
        if i == 4:
            # **** 1500m *****
            elevation = 1500
            temp_offset = 18
            slope = 0
            surface = "grass"
            label = "1500m ISA+18"
        if i == 3:
            # **** 1000m *****
            elevation = 1000
            temp_offset = 18
            slope = 0
            surface = "grass"
            label = "1000m ISA+18"
        if i == 5:
            # **** 1000m *****
            elevation = 2000
            temp_offset = 18
            slope = 0
            surface = "grass"
            label = "2000m ISA+18"

        acf = aircraft.Aircraft()
        W_list = np.linspace(acf.W_OE, acf.W_MTO, 25)

        acf = aircraft.Aircraft()

        TO_dist = []

        printed = False

        for W in W_list:
            dist = acf.takeoff_ground_run(W, elevation, temp_offset, slope, surface)[0]
            if i == 2 and dist > 500 and not printed:
                printed = True
                print(W_list[len(TO_dist) - 1])
                print(TO_dist[-1])
            TO_dist.append(dist)

        plt.plot(W_list / 9.80665, TO_dist, label=label)

    # plt.ylim(300, 900)
    plt.ylabel("Required runway length [m]")
    plt.xlabel("Gross takeoff weight [kg]")
    plt.xlim(acf.W_OE / 9.80665, acf.W_MTO / 9.80665)
    plt.legend()
    fig.axes[0].set_xticks([2800, 3000, 3200, 3400])
    plt.grid()
    plt.savefig("plots/takeoff.svg")
    plt.show()


if __name__ == "__main__":
    plot_variation()

# we can reduce wing area to 28m^2 if we
# use 50% payload and 50% fuel for 750m ISA+15 takeoff
# decrease MTOW max ROC to 4 m/s
# increase MTOW stall speed to 31 m/s (CS23)
