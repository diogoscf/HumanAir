import numpy as np
from HumanAir.FlightPerformance import aircraft
from HumanAir.FlightPerformance.helper import density
import matplotlib.pyplot as plt

g = 9.80665


class MAircraft:
    """
    Mission-Aircraft class.

    Wraps Aircraft object for use for mission analysis.
    """

    def __init__(self, elevation, acf=aircraft.Aircraft()):
        self.acf = acf

        #
        # constants
        #

        self.dT = self.acf.dT_default  # ISA temperature offset used throuhout flight [deg C]
        self.S = self.acf.S()  # clean config

        #
        # flight phase ID's
        #

        self.id_TO = 0
        self.id_acc = 1
        self.id_climb = 2
        self.id_cruise = 3
        self.id_descent = 4
        self.id_loiter = 5
        self.id_land = 6

        # index corresponds to ID
        self.phase_names = ["take-off", "acceleration", "climb", "cruise", "descent", "loiter", "landing"]

        #
        # variables for mission analysis
        #

        self.V = 0  # current true air speed [m/s]
        self.h = elevation  # current geopotential altitude [m]

        self.payload = self.acf.W_pl_no_pilot  # current payload [N]
        self.fuel = self.acf.W_MF  # current fuel weight [N]
        self.bat_cap = self.acf.max_bat_cap  # current battery capacity [Wh]

        # fuel weight in newtons that would be equivalent to battery capacity used
        self.bat_eq_fuel = 0

        self.ground_distance = 0  # distance flown [m]
        self.flight_time = 0  # time since start takeoff [s]

        #
        # lists to keep track of changing parameters
        #

        self.t_lst = []  # time since start takeoff [s]
        self.v_lst = []  # speed [m/s]
        self.h_lst = []  # altitude [m]
        self.s_lst = []  # ground distance [m]
        self.f_lst = []  # fuel weight [N]
        self.b_lst = []  # battery capacity [Wh]
        self.CL_lst = []  # lift coefficient [-]

        #
        # variables to keep track of fuel consumption during different phases
        # index corresponds to flight phase ID
        #

        self.phase_fuel = [0, 0, 0, 0, 0, 0, 0]
        self.phase_eq_fuel = [0, 0, 0, 0, 0, 0, 0]
        self.phase_bat = [0, 0, 0, 0, 0, 0, 0]

        #
        # variables for printing energy consumption per flight part
        # these are first initialized in the first takeoff
        #

        self.last_fuel = 0
        self.last_eq_fuel = 0
        self.last_bat_cap = 0

    def _update_lists(self):
        """
        Updates flight log that will be plotted with plot_flight().

        Returns
        -------
        None.

        """
        # print warnings if neccessary
        if len(self.b_lst) > 0 and self.b_lst[-1] > 0 and self.bat_cap < 0:
            print("***** WARNING: battery capacity below 0 ******")
        if len(self.f_lst) > 0 and self.f_lst[-1] > 0 and self.fuel < 0:
            print("***** WARNING: fuel level below 0 ******")

        # append current values to lists
        self.t_lst.append(self.flight_time)
        self.v_lst.append(self.V)
        self.h_lst.append(self.h)
        self.s_lst.append(self.ground_distance)
        self.f_lst.append(self.fuel)
        self.b_lst.append(self.bat_cap)
        if self.V > 0:
            self.CL_lst.append(self.CL)
        else:
            self.CL_lst.append(self.acf.CL_ground_TO)

    def plot_flight(self):
        """
        Plots the flight logs.

        Returns
        -------
        None.

        """
        fig, axs = plt.subplots(6, 1)

        fig.set_figheight(13)
        fig.set_figwidth(10)

        axs[0].plot(self.t_lst, self.h_lst)

        axs[1].plot(self.t_lst, self.v_lst)

        axs[2].plot(self.t_lst, self.f_lst)

        axs[3].plot(self.t_lst, self.b_lst)

        axs[4].plot(self.t_lst, self.s_lst)

        axs[5].plot(self.t_lst[1:-1], self.CL_lst[1:-1])

        plt.show()

    def _print_update(self, phase_id):
        """
        Prints fuel used for flight phase.

        Parameters
        ----------
        label : string
            A description of the flight segment to be printed with info.

        Returns
        -------
        None.

        """

        label = self.phase_names[phase_id]

        # get and print difference
        fuel_used = self.last_fuel - self.fuel
        eq_fuel_used = self.bat_eq_fuel - self.last_eq_fuel  # starts at zero so order reversed
        bat_cap_used = self.last_bat_cap - self.bat_cap
        print(
            f"** {fuel_used:8.1f} N fuel and"
            + f" {bat_cap_used:8.1f} Wh energy (={eq_fuel_used:5.1f} N fuel)"
            + f" used during {label}"
        )

        # update
        self.last_fuel = self.fuel
        self.last_eq_fuel = self.bat_eq_fuel
        self.last_bat_cap = self.bat_cap

    def print_energy_level(self):
        """
        Prints how much fuel/bat is left over.

        Returns
        -------
        None.

        """
        _, perc = self.calculate_fuel_reserve()
        b_1 = self.b_lst[0]
        eq_fuel = sum(self.phase_eq_fuel) - self.phase_eq_fuel[self.id_loiter]
        print("***************************************************************")
        print(f"{self.fuel:8.1f} N fuel left over ({perc:.1f}% fuel reserve)")
        print(f"{eq_fuel:5.1f} N eq fuel total")
        print(f"{self.bat_cap:8.1f} Wh battery left over ({self.bat_cap/b_1*100:.1f}%)")
        print("***************************************************************")

        for i in range(7):
            print(
                f"Total {self.phase_fuel[i]:8.1f} N fuel, \
{self.phase_bat[i]:8.1f} Wh energy \
({self.phase_eq_fuel[i]:5.1f} N equivalent) \
for {self.phase_names[i]}"
            )

        print("***************************************************************")
        print(f"Ground distance: {self.ground_distance:.1f} m  =  {self.ground_distance/1852:.1f} nm")
        print(f"Flight duration: {self.flight_time:.1f} s")
        print("***************************************************************")

    def calculate_fuel_reserve(self):
        """
        Calculates the required fuel reserve in newtons given (15% of trip fuel,
        in other words all fuel consumed except by loiter).

        Returns
        -------
        fuel_reserve_required: float
            Required fuel reserve to reach 15% trip fuel reserve
        percentage_reserve: float
            Current fuel level as % of trip fuel

        """

        # trip fuel = all fuel except loiter
        trip_fuel = sum(self.phase_fuel) - self.phase_fuel[self.id_loiter]
        trip_fuel += sum(self.phase_eq_fuel) - self.phase_eq_fuel[self.id_loiter]

        fuel_reserve_required = trip_fuel * 0.15

        percentage_reserve = self.fuel / trip_fuel * 100

        return fuel_reserve_required, percentage_reserve

    def _get_W(self):
        """Returns current gross weight [N]"""
        return self.acf.W_OE + self.payload + self.fuel

    W = property(fget=_get_W)

    def _get_rho(self):
        """Returns air density in current conditions"""
        return density(self.h, self.dT)

    rho = property(fget=_get_rho)

    def _get_CL(self):
        """Returns current CL assuming small flight path angle (L = W)"""
        return self.W / (0.5 * self.rho * self.V**2 * self.S)

    CL = property(fget=_get_CL)

    def _get_D(self):
        """Returns current drag force for clean config assuming small flight path angle"""
        return 0.5 * self.rho * self.V**2 * self.S * self.acf.CD(self.CL)

    D = property(fget=_get_D)

    def P_r(self):
        """Power required for steady level flight"""
        return self.D * self.V

    def _consume_energy(self, P_shaft, duration, electric, phase_id):
        """
        Consume fuel or battery capacity for given shaft power and duration

        Parameters
        ----------
        P_shaft : float
            Shaft power required [W].
        duration : float
            Duration of power consumption [s].
        electric : boolean
            Whether the electric motor is used.
        phase_id : integer
            Phase id of flight phase

        Returns
        -------
        None.

        """
        is_takeoff = phase_id == self.id_TO
        fuel_cons = self.acf.fuel_rate(P_shaft, is_takeoff) * duration
        if electric:
            # update battery capacity and equivalent fuel consumption
            bat_cons = self.acf.bat_cap_rate(P_shaft) * duration
            self.bat_cap -= bat_cons
            self.bat_eq_fuel += fuel_cons  # plus as it starts from 0
            self.phase_bat[phase_id] += bat_cons
            self.phase_eq_fuel[phase_id] += fuel_cons
        else:
            # update fuel level
            self.fuel -= fuel_cons
            self.phase_fuel[phase_id] += fuel_cons

    def fly_accelerate_const_alt(self, V_target, electric=False, time_step=0.5):
        """
        This function will either accelerate or decelerate the aircraft.
        depending on whether V_target is greater or smaller than the current
        velocity.

        Acceleration is performed with max const. throttle.

        Deceleration is performed with zero throttle.

        Parameters
        ----------
        V_target : float
            True air speed at end of de/acceleration.
        electric : boolean, optional
            Whether this is performed with the electric motor. The default is
            False.
        time_step : float, optional
            Time step for integration in seconds. The default is .5.

        Returns
        -------
        None.

        """
        # check if we need to change speed at all
        if np.isclose(self.V, V_target):
            return

        # check whether we need to accelerate or decelerate
        if self.V < V_target:
            accelerate = True
        else:
            accelerate = False

        #
        # integration
        #
        while (accelerate and self.V < V_target) or (not accelerate and self.V > V_target):
            if accelerate:
                # max constant thrust
                T = self.acf.T(self.V, self.h, self.dT, electric=electric)
                P_shaft = self.acf.P_shaft(self.h, self.dT, electric=electric)
                self._consume_energy(P_shaft, time_step, electric, self.id_acc)
            else:
                # zero throttle
                T = 0

            # calculate deceleration or acceleration, assumed constant during time_step
            a = g * (T - self.D) / self.W

            # update aircraft parameters
            self.ground_distance += (self.V + a * time_step / 2) * time_step  # avg speed over dt
            self.V += a * time_step
            self.flight_time += time_step
            self._update_lists()

        # done, print energy usage
        self._print_update(self.id_acc)

    def fly_const_climbrate(self, RC, h_target, electric=False, time_step=10):
        """
        Fly at constany climb rate. May be used for descent as well as ascent.
        In climb the speed is kept constant, in descent the speed may increase
        but it will not decrease.

        Parameters
        ----------
        RC : float
            Climb rate [m/s].
        h_target : float
            Target altitude [m].
        electric : boolean, optional
            Whether to use electric power. The default is False.
        time_step : float, optional
            Integration time step. The default is 10.

        Returns
        -------
        None.

        """
        # speed will not decrease below this
        V = self.V

        # climb
        if RC > 0:
            # determine ground speed, stays constant
            V_g = np.sqrt(V**2 - RC**2)  # ground speed (so we are taking into account flight path angle here)

            while h_target > self.h:
                # derived from RoC = (P_a - P_r) / W
                P = RC * self.W + self.V * self.D
                P_shaft = P / self.acf.prop_eff(self.V, self.h, self.dT)
                self._consume_energy(P_shaft, time_step, electric, self.id_climb)

                # update aircraft parameters
                self.ground_distance += V_g * time_step
                self.flight_time += time_step
                self.h += RC * time_step
                self._update_lists()

            # done, print energy consumed
            self._print_update(self.id_climb)

        # descent
        else:
            while h_target < self.h:
                # determine ground speed
                V_g = np.sqrt(V**2 - RC**2)  # ground speed (so we are taking into account flight path angle here)

                # derived from RoC = (P_a - P_r) / W
                P = RC * self.W + self.V * self.D

                if P > 0:
                    P_shaft = P / self.acf.prop_eff(self.V, self.h, self.dT)
                    self._consume_energy(P_shaft, time_step, electric, self.id_descent)
                else:
                    # accelerate
                    a = -P / self.V / self.W
                    self.V += a * time_step

                # update aircraft parameters
                self.ground_distance += V_g * time_step
                self.flight_time += time_step
                self.h += RC * time_step
                self._update_lists()

            # done, print energy consumed
            self._print_update(self.id_descent)

    def fly_const_V_const_alt(self, target_distance, time_step=10, electric=False):
        """
        Fly keeping current altitude and speed constant (good for cruise). When
        electric is set to true it automatically switches to the fuel engine
        when the battery is depleted (reaches max Depth of Discharge).

        Parameters
        ----------
        target_distance : float
            Ground distance to cover [m].
        time_step : float, optional
            Time step for integratin [s]. The default is 60.
        electric : boolean, optional
            Whether to use the electric motor. Once the battery is depleted
            it automatically switches to the fuel engine. The default is False.

        Returns
        -------
        None.

        """
        # print(f"CL before: {self.CL:.2f} CL/CD: {self.CL/self.acf.CD(self.CL):.2f} opt: {self.acf.LDmax:.2f}")
        # stores distance covered during this phase
        distance = 0

        # constant V means constant prop efficiency
        prop_eff = self.acf.prop_eff(self.V, self.h, self.dT)

        # constant alt = constant rho
        rho = self.rho

        #
        # integration
        #
        while distance < target_distance:
            # get required shaft power
            P_shaft = self.V * 0.5 * rho * self.V**2 * self.S * self.CL / 12 / prop_eff

            if electric:
                # make sure DoD is not exceeded
                energy = self.acf.bat_cap_rate(P_shaft) * time_step
                DoD = 1 - (self.bat_cap - energy) / self.acf.max_bat_cap

                if DoD >= self.acf.max_DoD_bat:
                    # switch to fuel
                    electric = False

            self._consume_energy(P_shaft, time_step, electric, self.id_cruise)

            # update ground distance of this segment
            distance += self.V * time_step

            # update aircraft parameters
            self.ground_distance += self.V * time_step
            self.flight_time += time_step
            self._update_lists()

        # print(f"CL after: {self.CL:.2f}")

        # done, print energy consumption
        self._print_update(self.id_cruise)

    def fly_const_CL_const_alt(self, duration, time_step=10, electric=False):
        """
        Loiter at constant CL at constant altitude without covering more ground
        distance. The current CL is maintained.

        Parameters
        ----------
        duration : float
            Duration of this flight segment [s].
        time_step : float, optional
            Time step for integration. The default is 60.
        electric : boolean, optional
            Whether electric power is used. The default is False.

        Returns
        -------
        None.

        """
        # const alt = const rho
        rho = self.rho

        # store CL to keep constant
        CL = self.CL

        # duration of current phase
        t = 0

        #
        # integrate
        #
        while t < duration:
            # get required shaft power, neglecting deceleration to maintain CL
            P_shaft = self.V * self.D / self.acf.prop_eff(self.V, self.h, self.dT)

            self._consume_energy(P_shaft, time_step, electric, self.id_loiter)

            # update speed to maintain CL
            self.V = np.sqrt(self.W / (0.5 * rho * self.S * CL))

            # update time
            t += time_step

            # update aircraft parameters
            # no ground distance covered
            self.flight_time += time_step
            self._update_lists()

        # done, print energy usage
        self._print_update(self.id_loiter)

    def take_off(self, elevation, surface, electric=False):
        """
        Make the aircraft take off.

        At the end the aircraft parameters are:
            h = airfield elevation
            V = liftoff speed
            t = takeoff duration

        Weight and battery capacity are also updated but fuel consumption is
        not taken into account during takeoff run (tiny effects anyway).

        Parameters
        ----------
        elevation : float
            Airfield elevation [m].
        surface : string
            Runway surface ("grass" or "paved).
        electric : boolean, optional
            Whether to use electric power. The default is False.

        Raises
        ------
        Exception
            Raised if gross weight exceeds MTOW.

        Returns
        -------
        None.

        """
        # check if MTOW is not exceeded
        if self.W > self.acf.W_MTO:
            raise Exception(f"Aircraft weight ({self.W:.0f} N) exceeds MTOW ({self.acf.W_MTO:.0f} N)")

        if self.flight_time == 0:
            # perform some variable initialisation
            self._update_lists()
            self.last_fuel = self.fuel
            self.last_eq_fuel = self.bat_eq_fuel
            self.last_bat_cap = self.bat_cap

        # get duration and final speed of takeoff phase
        t, V, _ = self.acf.takeoff_ground_run(self.W, elevation, self.dT, 0, surface, electric=electric, calc_time=True)

        # get max TO shaft power
        P_shaft = self.acf.P_shaft(elevation, self.dT, use_takeoff_power=True, electric=electric)

        self._consume_energy(P_shaft, t, electric, self.id_TO)

        # update aircraft parameters
        # note that there may be multiple takeoffs in a flight so we cannot
        # reset these parameters
        self.flight_time += t
        self.V = V
        self.h = elevation  # neglect what little alt has been gained
        self._update_lists()

        # print energy used
        self._print_update(self.id_TO)

    def land(self):
        """
        Lands aircraft.

        Takes into account fuel consumption by reverse thrust, then sets V to zero.

        Assumes fuel engine is used.
        """
        t = self.acf.landing_ground_distance(self.W, self.h, self.dT, 0, "grass", reversible_pitch=True, calc_time=True)
        P = self.acf.P_shaft(self.h, self.dT) * 0.4
        self._consume_energy(P, t, False, self.id_land)
        self.V = 0
        self.flight_time += t
        self._update_lists()
