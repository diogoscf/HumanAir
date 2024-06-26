# ISA Calculator by Diogo S.C. Fernandes (ST# 5446538)
import math

G_0 = 9.80665  # gravitational acceleration at sea-level
R_air = 287.05  # specific gas constant for air
T_sea = 288.15  # ISA temperature at sea level (15 Celsius)
p_sea = 101325  # ISA pressure at sea level (1 atm)
rho_sea = p_sea / (R_air * T_sea)

ISA = [
    (11000.0, -0.0065),  # troposphere
    (20000.0, 0),  # tropopause
    (32000.0, 0.0010),  # stratosphere I
    (47000.0, 0.0028),  # stratosphere II
    (51000.0, 0),  # stratopause
    (71000.0, -0.0028),  # mesosphere I
    (86000.0, -0.0020),  # mesosphere II
]


# Round to significant figures
def round_sig(x, sig):
    return round(x, sig - int(math.floor(math.log10(abs(x)))) - 1)


# Temperature Calculator
def t1(t0, a, h0, h1):
    return t0 + a * (h1 - h0)


# Pressure Calculator
def p1(p0, t0, a, h0, h1):
    if a != 0:
        T1 = t1(t0, a, h0, h1)
        return p0 * ((T1 / t0) ** (-G_0 / (a * R_air)))
    return p0 * ((math.e) ** ((-G_0 / (R_air * t0)) * (h1 - h0)))


# Pressure Calculator
def rho1(p, t):
    return p / (R_air * t)


# Print Results in User-Readable Format
def results(t, p, rho):  # pragma: no cover
    t, p, rho = round_sig(t, 5), round_sig(p, 5), round_sig(rho, 5)
    t_celsius = round(t - 273.15, 2)
    p_perc, rho_perc = round_sig((100 * p) / p_sea, 2), round_sig((100 * rho) / rho_sea, 2)
    print(f"\nTemperature: {t} K ({t_celsius} *C)")
    print(f"Pressure: {p} Pa ({p_perc}% SL)")
    print(f"Density: {rho} kg/m^3 ({rho_perc}% SL)\n")


# Main ISA logic
def isa(h, delta_T=0):
    t, p, h0 = T_sea, p_sea, 0
    for h_max, a in ISA:
        h1 = min(h_max, h)
        t, p, h0 = t1(t, a, h0, h1), p1(p, t, a, h0, h1), h1
        if h <= h_max:
            return t + delta_T, p, rho1(p, t + delta_T)

    print("ERROR: Height too large - The ISA only goes to a height of 86 km (282152 ft or FL 2821)")  # pragma: no cover
    exit(0)  # pragma: no cover


# Feet and Flight Level converters
def ft_to_m(x):  # pragma: no cover
    return x * 0.3048


def FL_to_m(x):
    return x * 0.3048 * 100.0


def main():  # pragma: no cover
    print("\n* * * ISA calculator * * *\n")
    print("\n1. Calculate ISA for altitude in meters")
    print("2. Calculate ISA for altitude in feet")
    print("3. Calculate ISA for altitude in flight level (FL)\n")
    choice = input("Enter your choice: ")
    if choice == "1":
        results(*isa(float(input("\nEnter altitude [m]: "))))
    elif choice == "2":
        results(*isa(ft_to_m(float(input("\nEnter altitude [ft]: ")))))
    elif choice == "3":
        results(*isa(FL_to_m(float(input("\nEnter altitude [FL]: ")))))
    else:
        print("Not a valid choice, please select one of the options\n")
        main()
        return

    print("Would you like to make another calculation? Select Y for yes, N for no (default: N).")
    continuation = input("Your selection [Y/N]: ")
    if continuation == "Y" or continuation == "y":
        print("\n")
        main()


if __name__ == "__main__":  # pragma: no cover
    main()
