import numpy as np
from HumanAir.FlightPerformance import aircraft

import matplotlib.pyplot as plt


def plot_variation():

    for i in [1, 2, 3, 4, 5]:
        if i == 1:
            # ***** design *****
            TORA = 500
            elevation = 750
            temp_offset = 18
            slope = 0
            surface = "grass"
        elif i == 2:
            # ***** aibai *****
            TORA = 415
            elevation = 2060
            temp_offset = 25
            slope = -11
            surface = "grass"
        elif i == 3:
            # ***** nomane *****
            TORA = 485
            elevation = 1810
            temp_offset = 23
            slope = 0
            surface = "grass"
        if i == 4:
            # **** sea level *****
            TORA = 500
            elevation = 0
            temp_offset = 0
            slope = 0
            surface = "grass"
        if i == 5:
            # **** nairobi *****
            TORA = 500
            elevation = 1800
            temp_offset = 21
            slope = 0
            surface = "grass"

        acf = aircraft.Aircraft()
        W_list = np.linspace(acf.W_OE, acf.W_MTO, 25)

        plt.figure(figsize=(8, 5))

        for S in np.arange(20, 38, step=2):

            acf = aircraft.Aircraft()
            acf.S_clean = S

            TO_dist = []

            for W in W_list:
                TO_dist.append(acf.takeoff_ground_run(W, elevation, temp_offset, slope, surface)[0])

            plt.plot(W_list, TO_dist, label=f"S={S} m^2")

        plt.title(f"Takeoff at {elevation}m ISA+{temp_offset}, {slope}% slope, {surface} runway")
        # plt.ylim(300, 900)
        plt.ylabel("TO ground run [m]")
        plt.xlabel("Gross takeoff weight [N]")
        plt.axhline(TORA, color="black")
        plt.legend()
        plt.show()


def vary_area():

    # ***** design *****
    TORA = 500
    elevation = 750
    temp_offset = 18
    slope = 0
    surface = "grass"

    # ***** aibai *****
    # TORA        = 415
    # elevation   = 2060
    # temp_offset = 25
    # slope       = -11
    # surface     = "grass"

    # ***** nomane *****
    # TORA        = 485
    # elevation   = 1810
    # temp_offset = 23
    # slope       = 0
    # surface     = "grass"

    # **** sea level *****
    # TORA        = 500
    # elevation   = 0
    # temp_offset = 0
    # slope       = 0
    # surface     = "grass"

    for S in [26]:  # np.arange(20, 40, step=1):

        acf = aircraft.Aircraft()

        # reduced_weight = acf.W_OE + acf.W_pl_no_pilot + acf.W_MF

        acf.S_clean = S

        climbrate = acf.RC_max(acf.W_MTO, 0, 0)
        gradient_cl = acf.climb_slope_max(acf.W_MTO, 0, 0, gear="up", flaps="up")
        gradient_TO = acf.climb_slope_max(acf.W_MTO, 0, 0, gear="up", flaps="TO")
        gradient_ld = acf.climb_slope_max(acf.W_MTO, 0, 0, gear="down", flaps="land")
        V_S0 = acf.stall_speed(acf.W_MTO, 0, 0, flaps="land")
        V_S1 = acf.stall_speed(acf.W_MTO, 0, 0, flaps="up")
        TO_dist = acf.takeoff_ground_run(acf.W_MTO, elevation, temp_offset, slope, surface)
        land_dist = acf.landing_ground_distance(
            acf.W_MTO, elevation, temp_offset, -slope, surface, reversible_pitch=True
        )

        pass_cr = "pass" if climbrate[0] >= 5 else "FAIL"
        pass_gr_TO = "pass" if gradient_TO[0] >= 4 else "FAIL"
        pass_gr_ld = "pass" if gradient_ld[0] >= 3 else "FAIL"
        pass_vs = "pass" if V_S0 <= 25 else "FAIL"
        pass_to = "pass" if TO_dist[0] <= TORA else "FAIL"
        pass_ld = "pass" if land_dist <= TORA else "FAIL"

        print(
            f"""
              ****************** S = {S} ***********************
              Climb rate           = {climbrate[0]:>8.2f} m/s    {pass_cr}   at V={climbrate[1]:.2f} m/s
              Climb gradient clean = ({gradient_cl[0]:>8.2f}) % at V={gradient_cl[1]:.2f} m/s
              Climb gradient TO    = ({gradient_TO[0]:>8.2f}) %      {pass_gr_TO}   at V={gradient_TO[1]:.2f} m/s
              Climb gradient land  = ({gradient_ld[0]:>8.2f}) %      {pass_gr_ld}   at V={gradient_ld[1]:.2f} m/s
              V_S0                 = {V_S0:>8.2f} m/s    {pass_vs}
              V_S1                 = {V_S1:>8.2f} m/s
              Min. TORA take-off   = {TO_dist[0]:>8.0f} m      {pass_to}
              Min. TORA landing    = {land_dist:>8.0f} m      {pass_ld}
              Runway with {TORA}m TORA at {elevation}m ISA+{temp_offset}, {slope}% slope, {surface} runway
              *************************************************
              ! Note that the climb gradient calculations are
                inaccurate because of the drag approximation  !
              """
        )


def vary_weight():

    # ***** design *****
    TORA = 500
    elevation = 750
    temp_offset = 18
    slope = 0
    surface = "grass"

    # ***** aibai *****
    # TORA        = 415
    # elevation   = 2060
    # temp_offset = 25
    # slope       = -11
    # surface     = "grass"

    # ***** nomane *****
    # TORA        = 485
    # elevation   = 1810
    # temp_offset = 23
    # slope       = -6
    # surface     = "grass"

    # **** sea level *****
    # TORA        = 500
    # elevation   = 0
    # temp_offset = 0
    # slope       = 0
    # surface     = "grass"

    acf = aircraft.Aircraft()

    acf.S_clean = 26

    for W in np.arange(acf.W_OE, acf.W_MTO + 500, step=500):

        TO_dist = acf.takeoff_ground_run(W, elevation, temp_offset, slope, surface)[0]
        land_dist = acf.landing_ground_distance(W, elevation, temp_offset, -slope, surface, reversible_pitch=True)

        pass_to = "pass" if TO_dist <= TORA else "FAIL"
        pass_ld = "pass" if land_dist <= TORA else "FAIL"

        fuel_and_payload_perc = (W - acf.W_OE) / (acf.W_pl_no_pilot + acf.W_MF) * 100

        print(
            f"""
              *********** W = {W:.0f} ***** S = {acf.S} **********
              Fuel and payload = {fuel_and_payload_perc:>8.2f} %
              Takeoff run      = {TO_dist:>8.0f} m      {pass_to}
              Landing roll     = {land_dist:>8.0f} m      {pass_ld}
              Runway with {TORA}m TORA at {elevation}m ISA+{temp_offset}, {slope}% slope, {surface} runway
              *************************************************
              """
        )


if __name__ == "__main__":
    # plot_variation()
    vary_area()
    # vary_weight()

# we can reduce wing area to 28m^2 if we
# use 50% payload and 50% fuel for 750m ISA+15 takeoff
# decrease MTOW max ROC to 4 m/s
# increase MTOW stall speed to 31 m/s (CS23)
