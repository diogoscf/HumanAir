import matplotlib.pyplot as plt

# import aircraft
# import mission_evaluation

g = 9.80665
pl1, pl2, pl3 = 630, 540, 0
extra_payload = (pl1 - pl2) * g


#
# hybrid, four legs
#

r1, r2, r3 = 466.2, 668.9, 819.9

# # point 1, max payload, reduced fuel
# acf = aircraft.Aircraft()
# acf.W_pl_no_pilot += extra_payload
# acf.W_MF -= extra_payload
# r1 = mission_evaluation.calculate_range(acf, guess_tot_range_nm=590) / 1852

# # point 2, design payload, max fuel
# acf = aircraft.Aircraft()
# r2 = mission_evaluation.calculate_range(acf, guess_tot_range_nm=636) / 1852

# # point 3, zero payload, max fuel
# acf = aircraft.Aircraft()
# acf.W_pl_no_pilot = pl3*g
# r3 = mission_evaluation.calculate_range(acf, guess_tot_range_nm=670) / 1852


#
# hybrid, one leg, has to be determined manually
#

# ro1, ro2, ro3 = 431.5, 633.4, 788.7 # TODO: not updated for corrected design.json


#
# electric only, one leg
#

rb1, rb2, rb3 = 108.8, 111.2, 132.5

# # point 1, max payload, no fuel
# acf = aircraft.Aircraft()
# acf.W_pl_no_pilot += extra_payload
# acf.W_MF = 420 # reserve fuel, when changes are made, check if still correct
# rb1 = mission_evaluation.calculate_range(acf, guess_tot_range_nm=103, num_legs=1, only_electric=True) / 1852

# # point 2, design payload, no fuel
# acf = aircraft.Aircraft()
# acf.W_MF = 410 # reserve fuel, when changes are made, check if still correct
# rb2 = mission_evaluation.calculate_range(acf, guess_tot_range_nm=105, num_legs=1, only_electric=True) / 1852

# # point 3, zero payload, no
# acf = aircraft.Aircraft()
# acf.W_pl_no_pilot = pl3*g
# acf.W_MF = 330 # reserve fuel, when changes are made, check if still correct
# rb3 = mission_evaluation.calculate_range(acf, guess_tot_range_nm=110, num_legs=1, only_electric=True) / 1852


print(f"Range max. payload, reduced fuel, four legs: {r1:.1f} nm")
print(f"Range des. payload, max fuel, four legs:     {r2:.1f} nm")
print(f"Range no payload, max fuel, four legs:       {r3:.1f} nm")

# print(f"Range max. payload, reduced fuel, one leg:   {ro1:.1f} nm")
# print(f"Range des. payload, max fuel, one leg:       {ro2:.1f} nm")
# print(f"Range no payload, max fuel, one leg:         {ro3:.1f} nm")


print(f"Range max. payload, no fuel, one leg:        {rb1:.1f} nm")
print(f"Range des. payload, no fuel, one leg:        {rb2:.1f} nm")
print(f"Range no payload, no fuel, one leg:          {rb3:.1f} nm")


plt.figure(figsize=(6, 3))
plt.plot([0, r1, r2, r3], [pl1, pl1, pl2, pl3], color="b", label="Hybrid, four flight legs")
# plt.plot([0, ro1,ro2,ro3],[pl1, pl1,pl2,pl3], color="b", label="Hybrid, one flight leg")
plt.plot([0, rb1, rb2, rb3], [pl1, pl1, pl2, pl3], color="g", label="Electric, one flight leg")
plt.scatter(600, 540, marker="x", color="r", label="Design point")
plt.xlabel("Range [NM]")
plt.ylabel("Payload [kg]")
plt.ylim(0, pl1 * 1.1)
plt.xlim(0, r3 * 1.1)
plt.grid()
plt.legend()
plt.savefig("plots/payload_range.svg")
plt.show()
