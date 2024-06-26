import numpy as np
import matplotlib.pyplot as plt

from HumanAir.FlightPerformance import aircraft

h_max = 7000
h_step = 100


def plot_climb_rate(acf):

    W = acf.W_MTO
    dT = 0

    h_list = np.arange(0, h_max, h_step)
    RC_max = np.zeros(((len(h_list), 2)))
    P_a = np.zeros(len(h_list))
    P_r_min = np.zeros(len(h_list))

    service_ceiling = 0

    for i, h in enumerate(h_list):
        RC_max[i] = acf.RC_max(W, h, dT)
        P_a[i] = acf.P_a(h, dT)
        P_r_min[i] = acf.P_r_min(W, h, dT)

        # service ceiling commonly defined as altitude with 100ft/min (0.508m/s) max climbrate
        # FAA Pilots Handbook of Aeronautical Knowledge
        # also Roskam pt 7 p129
        if RC_max[i][0] < 0.508 and service_ceiling < 0.01 and i > 0:
            service_ceiling = h_list[i - 1]

    print(f"Sea level MTOW max climb rate, clean config: {RC_max[0][0]:.3f} m/s at V={RC_max[0][1]:.2f} m/s")
    print(f"Service ceiling, MTOW, clean config: {service_ceiling:.0f} m")

    plt.figure(figsize=(10, 7))
    plt.plot(RC_max[:, 0], h_list, label="Maximum Climb Rate")
    plt.plot(P_a / 100000, h_list, color="g", label="Power available /100000")
    plt.plot(P_r_min / 100000, h_list, color="r", label="Minimum required power /100000")
    plt.axhline(y=service_ceiling, color="grey", label="Service Ceiling")
    plt.xlabel(r"Climb Rate (m/s)")
    plt.ylabel("Altitude (m)")
    plt.ylim(0, h_max)
    # plt.xlim(0,RC_max[0]*1.05)
    plt.legend()
    plt.text(0, h_max + 100, "RC=0 when CL>1")
    # plt.subplots_adjust(right=0.75)
    # plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.savefig("plots/climb_rate.svg")
    plt.show()


def plot_climb_gradient(acf):

    W = acf.W_MTO
    dT = 0

    h_list = np.arange(0, h_max, h_step)
    max_gradient_clean = np.zeros((len(h_list), 2))
    max_gradient_TO = np.zeros((len(h_list), 2))
    max_gradient_land = np.zeros((len(h_list), 2))

    for i, h in enumerate(h_list):
        max_gradient_clean[i] = acf.climb_slope_max(W, h, dT)
        max_gradient_TO[i] = acf.climb_slope_max(W, h, dT, gear="up", flaps="TO")
        max_gradient_land[i] = acf.climb_slope_max(W, h, dT, gear="down", flaps="land")

    print(
        "Sea level MTOW max climb gradient, clean config: "
        + f"{max_gradient_clean[0][0]:.3f} % at V={max_gradient_clean[0][1]:.2f}"
    )
    print(
        "Sea level MTOW max climb gradient, gear up, take-off flaps: "
        + f"{max_gradient_TO[0][0]:.3f} % at V={max_gradient_TO[0][1]:.2f}"
    )
    print(
        "Sea level MTOW max climb gradient, gear down, landing flaps:"
        + f"{max_gradient_land[0][0]:.3f} % at V={max_gradient_land[0][1]:.2f}"
    )

    plt.figure(figsize=(10, 7))
    plt.plot(max_gradient_clean[:, 0], h_list, color="b", label="MTOW, gear up, flaps up")
    plt.plot(max_gradient_TO[:, 0], h_list, color="purple", label="MTOW, gear up, take-off flaps")
    plt.plot(max_gradient_land[:, 0], h_list, color="lightblue", label="MTOW, gear up, landing flaps")
    plt.xlabel(r"Maximum Climb Gradient with CL<=1 (%)")
    plt.ylabel("Altitude (m)")
    plt.ylim(0, h_max)
    plt.xlim(0, max(max_gradient_clean[0][0], max_gradient_TO[0][0], max_gradient_land[0][0]) * 1.05)
    plt.legend()
    plt.text(0, h_max + 100, "Not the actual max climb gradient since CL is capped at 1")

    # plt.subplots_adjust(right=0.75)
    # plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.savefig("plots/climb_gradient.svg")
    plt.show()


def plot_stall_speed(acf):

    W = acf.W_MTO
    dT = 0

    h_list = np.arange(0, h_max, h_step)
    stall_clean = np.zeros(len(h_list))
    stall_TO = np.zeros(len(h_list))
    stall_land = np.zeros(len(h_list))
    V_max = np.zeros(len(h_list))

    stop_V_max = False
    last_V_max = (0, 0)

    for i, h in enumerate(h_list):
        stall_clean[i] = acf.stall_speed(W, h, dT)
        stall_TO[i] = acf.stall_speed(W, h, dT, flaps="TO")
        stall_land[i] = acf.stall_speed(W, h, dT, flaps="land")

        if not stop_V_max:
            V_max[i] = acf.V_max(W, h, dT)  # after max speed crosses stall speed graph gets weird
            if V_max[i] < stall_clean[i]:
                stop_V_max = True
                last_V_max = ([V_max[i - 1], V_max[i]], [h_list[i - 1], h_list[i]])
                V_max[i] = 0
        else:
            V_max[i] = 0

    print(f"Sea level MTOW stall speed, clean config: {stall_clean[0]:.3f} m/s")
    print(f"Sea level MTOW stall speed, take-off flaps: {stall_TO[0]:.3f} m/s")
    print(f"Sea level MTOW stall speed, landing flaps: {stall_land[0]:.3f} m/s")
    print(f"Sea level MTOW maximum speed, clean config: {V_max[0]:.3f} m/s")

    plt.figure(figsize=(10, 7))
    plt.plot(stall_clean, h_list, color="b", label="Stall, MTOW, flaps up")
    plt.plot(stall_TO, h_list, color="purple", label="Stall, MTOW, take-off flaps")
    plt.plot(stall_land, h_list, color="lightblue", label="Stall, MTOW, landing flaps")
    plt.plot(V_max, h_list, color="red", label="Maximum velocity, MTOW, clean configuration")
    plt.scatter(last_V_max[0], last_V_max[1], s=100, color="red", label="Last V_max before crossing stall speed")
    plt.xlabel("Velocity (m/s)")
    plt.ylabel("Altitude (m)")
    plt.ylim(0, h_max)
    plt.xlim(0, 200)
    # plt.xlim(0,max(stall_clean[0], stall_TO[0], stall_land[0])*1.05)
    plt.legend()

    plt.text(0, h_max + 100, "V_max not accurate for high altitudes with C_L > 1")

    # plt.subplots_adjust(right=0.75)
    # plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.savefig("plots/stall_speed.svg")
    plt.show()


def plot_cruise_power_setting(acf):

    W = acf.W_MTO
    dT = 0

    h_list = np.arange(0, h_max, h_step)
    powersetting = np.zeros(len(h_list))

    for i, h in enumerate(h_list):
        powersetting[i] = acf.power_setting(acf.V_cruise, W, h, dT)

    std_cruise = acf.power_setting(acf.V_cruise, W, acf.h_cruise, acf.dT_default)
    max_cruise = acf.power_setting(acf.V_cruise, W, acf.h_cruise_max, acf.dT_default)
    print(f"Power setting at MTOW, V={acf.V_cruise}, h={acf.h_cruise} ISA+{acf.dT_default}: {std_cruise:.3f}")
    print(f"Power setting at MTOW, V={acf.V_cruise}, h={acf.h_cruise_max} ISA+{acf.dT_default}: {max_cruise:.3f} ")

    plt.figure(figsize=(10, 7))
    plt.plot(powersetting, h_list, color="b", label=f"MTOW, {acf.V_cruise:.1f} m/s")
    plt.xlabel("Power setting (-)")
    plt.ylabel("Altitude (m)")
    plt.ylim(0, h_max)
    plt.xlim(powersetting[0] * 0.9, 1.0)
    plt.legend()
    plt.text(powersetting[0] * 0.9, h_max + 100, "Not accurate for high altitudes with C_L > 1")
    plt.savefig("plots/cruise_power_setting.svg")
    plt.show()


def plot_take_off_performance(acf):
    res = acf.takeoff_ground_run(acf.W_MTO, 750, 18, 0, "grass")
    print(f"750m ISA+18 MTOW grass takeoff ground run: {res[0]:.2f} m with precision {res[1]}")
    # res2 = acf.takeoff_ground_run(acf.W_MTO, 0, 0, -5, "grass")
    # print(f"Takeoff performance: {res2[0]:.2f} with precision {res2[1]}")


def plot_landing_performance(acf):
    dist = acf.landing_ground_distance(acf.W_MTO, 750, 18, 0, "grass", reversible_pitch=True)
    print(f"750m ISA+18 MTOW grass with thrust reverse landing distance: {dist:.2f} m")
    dist = acf.landing_ground_distance(acf.W_MTO, 750, 18, 0, "grass", reversible_pitch=False)
    print(f"750m ISA+18 MTOW grass without thrust reverse landing distance: {dist:.2f} m")


def plot_thrust(acf):
    v_lst = np.arange(0, acf.V_max(acf.W_OE, 0, 0), 1)
    T_lst = []
    for V in v_lst:
        T_lst.append(acf.T(V, 0, 0, use_takeoff_power=True))
    plt.figure(figsize=(10, 7))
    plt.plot(v_lst, T_lst)
    plt.ylim(0, T_lst[0] * 1.1)
    plt.xlim(0, v_lst[-1] * 1.1)
    plt.xlabel("Velocity [m/s]")
    plt.ylabel("Thrust [N]")
    plt.savefig("plots/thrust_speed.svg")
    plt.show()

    # print(f"Max thrust: {T_lst[0]:.2f}")
    # print(f"Min thrust: {T_lst[-1]:.2f}")

    v_lst = np.arange(0, acf.V_max(acf.W_OE, 0, 0) + 1, 1)
    T_lst = []
    for V in v_lst:
        T_lst.append(acf.prop_eff(V, 0, 0, use_takeoff_power=True))
    plt.figure(figsize=(10, 7))
    plt.plot(v_lst, T_lst)
    plt.ylim(0, T_lst[-1] * 1.1)
    plt.xlim(0, v_lst[-1] * 1.1)
    plt.xlabel("Velocity [m/s]")
    plt.ylabel("Prop efficiency [-]")
    plt.savefig("plots/prop_eff.svg")
    plt.show()

    print(f"Max eff: {T_lst[-1]:.2f}")


def plot_RC_F_max(acf):
    h_list = np.arange(0, h_max, h_step)
    V_list1 = np.zeros(len(h_list))
    V_list2 = np.zeros(len(h_list))

    for i, h in enumerate(h_list):
        _, V_list1[i] = acf.RC_max(acf.W_MTO, h, 0)
        V_list2[i] = acf.V_Prmin(acf.W_MTO, h, 0)

    V_lst = np.arange(acf.stall_speed(acf.W_MTO, 0, 0), acf.V_max(acf.W_MTO, 0, 0) + 1, 1)
    P_r_lst = np.zeros(len(V_lst))
    for i, V in enumerate(V_lst):
        P_r_lst[i] = acf.P_r(acf.W_MTO, V, 0, 0)

    plt.figure(figsize=(10, 7))
    _, ax1 = plt.subplots()
    ax1.plot(V_list1, h_list, label="Speed RCmax")
    ax1.plot(V_list2, h_list, label="Speed P_r_min")
    ax1.set_xlabel("Velocity [m/s]")
    ax1.set_ylabel("Altitude [m]")

    ax2 = ax1.twinx()
    ax2.set_ylabel("Power required [W]")
    ax2.plot(V_lst, P_r_lst, label="Power required")

    plt.legend()
    plt.show()


if __name__ == "__main__":
    acf = aircraft.Aircraft()
    # print(acf.RC_max(acf.W_MTO, 750, 18))
    # plot_climb_rate(acf)
    # plot_RC_F_max(acf)
    # plot_climb_gradient(acf)
    # plot_stall_speed(acf)
    plot_cruise_power_setting(acf)
    # plot_take_off_performance(acf)
    # plot_landing_performance(acf)
    # plot_thrust(acf)
