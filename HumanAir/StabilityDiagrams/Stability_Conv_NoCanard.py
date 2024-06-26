import HumanAir.LoadingDiagram.Plotting as plt
import matplotlib.pyplot as plot
import numpy as np

SSh = 0.170

xlst = np.arange(-0.8, 1.1, 0.1)
Stability = [
    -0.34337479816698,
    -0.310284854025373,
    -0.277194909883765,
    -0.244104965742158,
    -0.21101502160055,
    -0.177925077458942,
    -0.144835133317335,
    -0.111745189175727,
    -0.0786552450341198,
    -0.0455653008925122,
    -0.0124753567509046,
    0.0206145873907029,
    0.0537045315323105,
    0.0867944756739181,
    0.119884419815526,
    0.152974363957133,
    0.186064308098741,
    0.219154252240348,
    0.252244196381956,
]
Stability_NoMargin = [
    -0.359919770237784,
    -0.326829826096177,
    -0.293739881954569,
    -0.260649937812961,
    -0.227559993671354,
    -0.194470049529746,
    -0.161380105388139,
    -0.128290161246531,
    -0.0952002171049235,
    -0.062110272963316,
    -0.0290203288217084,
    0.00406961531989915,
    0.0371595594615067,
    0.0702495036031143,
    0.103339447744722,
    0.136429391886329,
    0.169519336027937,
    0.202609280169545,
    0.235699224311152,
]
Controllability = [
    1.30436194472079,
    1.19599824956778,
    1.08763455441477,
    0.979270859261763,
    0.870907164108754,
    0.762543468955746,
    0.654179773802737,
    0.545816078649728,
    0.437452383496719,
    0.32908868834371,
    0.220724993190701,
    0.112361298037692,
    0.0039976028846827,
    -0.104366092268326,
    -0.212729787421335,
    -0.321093482574344,
    -0.429457177727353,
    -0.537820872880362,
    -0.646184568033371,
]
Min = 0.25
Max = 0.75
ylst = np.arange(SSh - 0.02, SSh + 0.03, 0.01)

xline = np.linspace(Min, Max, 2)
yline = np.ones(len(xline)) * SSh

plot.figure(figsize=(10, 7))

plt.Ploty(xlst, Stability, "Stability", colour="green", linestyle="solid")
plt.Ploty(xlst, Stability_NoMargin, "Stability_NoMargin", colour="green", linestyle="dashed")
plt.Ploty(xlst, Controllability, "Controllability", colour="blue", linestyle="solid")
plt.Plotx(Min, ylst, "CG Excursion", colour="red", linestyle="solid")
plt.Plotx(Max, ylst, "", colour="red", linestyle="solid")
plot.plot(xline, yline, color="red", linestyle="solid")
plot.fill_between(xlst, 0, Stability, color="red", alpha=0.1)
plot.fill_between(xlst, 0, Controllability, color="red", alpha=0.1)

plot.xlabel(r"$X_{cg}/MAC$", fontsize=12, loc="right")
plot.ylabel(r"$S_h/S$", fontsize=12, loc="top")

plot.axhline(color="black", lw=0.5)
plot.axvline(color="black", lw=0.5)
plot.legend()
plot.xlim((0, 1))
plot.ylim((0, 0.5))
plot.savefig("Stability_Conv_NoCanard.png")
plot.savefig("Stability_Conv_NoCanard.svg")
plot.show()
