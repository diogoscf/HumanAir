# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 10:41:53 2024

@author: Alex
"""

from HumanAir.FlightPerformance import aircraft, mission_performance


def calculate_range(
    aircraft=aircraft.Aircraft(),
    guess_tot_range_nm=636,  # [nm]
    better_batteries=False,
    airfield_elevation=750,  # [m]
    cruise_altitude=3000,  # [m]
    loiter_altitude=1200,  # [m]
    temp_offset=18,  # [degC]
    climb_rate=2.5,  # [m/s]
    cruise_speed=60,  # [m/s]
    loiter_duration=75 * 60,  # [s]
    num_legs=4,
    only_electric=False,
):

    reserve = -10

    iterations = 0

    guess = (guess_tot_range_nm - 42 * num_legs) * 1852 / num_legs  # climb and descent add about 42nm per leg

    while (not only_electric and (reserve < 14.95 or reserve > 15.1)) or (
        only_electric and (reserve < 0 or reserve > 1)
    ):
        iterations += 1
        if iterations > 15:
            raise Exception(f"Too many iterations: {iterations}")

        macf = perform_flight(
            aircraft=aircraft,
            better_batteries=better_batteries,
            airfield_elevation=airfield_elevation,
            cruise_altitude=cruise_altitude,
            loiter_altitude=loiter_altitude,
            temp_offset=temp_offset,
            climb_rate=climb_rate,
            cruise_speed=cruise_speed,
            loiter_duration=loiter_duration,
            leg_dist_m=guess,
            num_legs=num_legs,
        )

        req_fuel_res, reserve = macf.calculate_fuel_reserve()

        if not only_electric:

            initial_fuel = macf.f_lst[0] + macf.bat_eq_fuel

            actual_fuel_consumed = initial_fuel - macf.fuel
            target_fuel_consumed = initial_fuel - req_fuel_res

            diff = actual_fuel_consumed / target_fuel_consumed

        else:

            initial_cap = macf.b_lst[0]

            minim = initial_cap * (1 - macf.acf.max_DoD_bat)
            diff = (macf.bat_cap - minim) / (initial_cap - minim)

            reserve = diff * 100

            diff *= -1
            diff += 1
            print(macf.bat_cap)
        guess /= diff

    return macf.ground_distance


def perform_flight(
    aircraft=aircraft.Aircraft(),
    better_batteries=False,
    airfield_elevation=750,  # [m]
    cruise_altitude=3000,  # [m]
    loiter_altitude=1200,  # [m]
    temp_offset=18,  # [degC]
    climb_rate=2.5,  # [m/s]
    cruise_speed=60,  # [m/s]
    loiter_duration=75 * 60,  # [s]
    leg_dist_m=201405,  # gess of cruise distance per leg [m]
    num_legs=4,
):

    macf = mission_performance.MAircraft(airfield_elevation, acf=aircraft)
    macf.dT = temp_offset

    if better_batteries:
        better_battery_factor = 685 / 350
        macf.acf.max_bat_cap = macf.acf.max_bat_cap * better_battery_factor
        macf.bat_cap = macf.acf.max_bat_cap

    #
    # flight
    #

    for i in range(num_legs):
        # whether to perform phases (except cruise) electrically
        perform_electric = num_legs == 1 or i == num_legs - 2

        # takeoff
        print(
            "**************** Takeoff distance: "
            + f"{macf.acf.takeoff_ground_run(macf.W, 750, 18, 0, 'grass')[0]:.2f} m"
            + " *******************"
        )
        macf.take_off(airfield_elevation, "grass", electric=perform_electric)

        V_start_climb = 50
        macf.fly_accelerate_const_alt(V_start_climb, electric=perform_electric)

        # climb
        macf.fly_const_climbrate(climb_rate, cruise_altitude, electric=perform_electric)

        # accelerate to cruise speed
        macf.fly_accelerate_const_alt(cruise_speed, electric=perform_electric)

        # cruise
        macf.fly_const_V_const_alt(
            leg_dist_m, electric=i >= num_legs - 2
        )  # if we have batcap left over this will use it also during fourth cruise

        if i < num_legs - 1:
            # descent
            macf.fly_const_climbrate(-climb_rate, airfield_elevation, electric=perform_electric)

            # land
            macf.land()
        else:
            # descent to loiter alt
            macf.fly_const_climbrate(-climb_rate, loiter_altitude, electric=perform_electric)

            # decelerate or accelerate to loiter speed
            V_loiter = macf.acf.V_Prmin(macf.W, macf.h, macf.dT)
            macf.fly_accelerate_const_alt(V_loiter)

            # loiter
            macf.fly_const_CL_const_alt(loiter_duration)

            # descent
            macf.fly_const_climbrate(-climb_rate, airfield_elevation, electric=False)

            # land
            macf.land()

    # plot and print

    macf.plot_flight()
    macf.print_energy_level()

    return macf


def perform_single_leg_hybrid_flight(
    aircraft=aircraft.Aircraft(),
    better_batteries=False,
    airfield_elevation=750,  # [m]
    cruise_altitude=3000,  # [m]
    loiter_altitude=1200,  # [m]
    temp_offset=18,  # [degC]
    climb_rate=2.5,  # [m/s]
    cruise_speed=60,  # [m/s]
    loiter_duration=75 * 60,  # [s]
):

    macf = mission_performance.MAircraft(airfield_elevation, acf=aircraft)
    macf.payload += 90 * 9.80665
    macf.fuel -= 90 * 9.80665
    macf.dT = temp_offset

    if better_batteries:
        better_battery_factor = 685 / 350
        macf.acf.max_bat_cap = macf.acf.max_bat_cap * better_battery_factor
        macf.bat_cap = macf.acf.max_bat_cap

    #
    # flight
    #

    # takeoff
    print(
        "**************** Takeoff distance:"
        + f" {macf.acf.takeoff_ground_run(macf.W, 750, 18, 0, 'grass')[0]:.2f} m"
        + " *******************"
    )
    macf.take_off(airfield_elevation, "grass", electric=False)

    # acceleration to ideal CL for climbing
    V_start_climb = 50  # macf.acf.V_Prmin(macf.W, macf.h, macf.dT)
    macf.fly_accelerate_const_alt(V_start_climb, electric=False)

    # climb
    macf.fly_const_climbrate(climb_rate, cruise_altitude, electric=False)

    # accelerate to cruise speed
    macf.fly_accelerate_const_alt(cruise_speed, electric=False)

    # cruise
    macf.fly_const_V_const_alt(1852 * 278, electric=False)  # 278 480 610
    macf.fly_const_V_const_alt(1852 * 99, electric=True)  # 99 99 125

    # descent to loiter alt
    macf.fly_const_climbrate(-climb_rate, loiter_altitude, electric=False)

    # decelerate or accelerate to loiter speed
    V_loiter = macf.acf.V_Prmin(macf.W, macf.h, macf.dT)
    macf.fly_accelerate_const_alt(V_loiter)

    # loiter
    macf.fly_const_CL_const_alt(loiter_duration)

    # descent
    macf.fly_const_climbrate(-climb_rate, airfield_elevation, electric=False)

    # land
    macf.land()

    # plot and print

    macf.plot_flight()
    macf.print_energy_level()

    return macf


if __name__ == "__main__":
    # calculate_range()
    # perform_single_leg_hybrid_flight()
    perform_flight()

    # acf = aircraft.Aircraft()
    # acf.W_pl_no_pilot -= 3*90*9.80665
    # acf.W_MF += 32988-acf.W_MTO + 3*90*9.80665
    # calculate_range(acf, num_legs=2)
