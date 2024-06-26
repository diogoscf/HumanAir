import HumanAir.LoadingDiagram.Plotting as plt
import matplotlib.pyplot as plot
import numpy as np

SSh = 0.218

xlst = np.arange(-1.4, 1.3, 0.1)
Stability = [
    0.593573508771592,
    0.55681726681827,
    0.520061024864947,
    0.483304782911624,
    0.446548540958302,
    0.409792299004979,
    0.373036057051656,
    0.336279815098334,
    0.299523573145011,
    0.262767331191688,
    0.226011089238366,
    0.189254847285043,
    0.15249860533172,
    0.115742363378398,
    0.0789861214250749,
    0.0422298794717522,
    0.00547363751842955,
    -0.0312826044348931,
    -0.0680388463882158,
    -0.104795088341538,
    -0.141551330294861,
    -0.178307572248184,
    -0.215063814201507,
    -0.251820056154829,
    -0.288576298108152,
    -0.325332540061475,
    -0.362088782014797,
]
Stability_NoMargin = [
    0.611951629748254,
    0.575195387794931,
    0.538439145841608,
    0.501682903888286,
    0.464926661934963,
    0.42817041998164,
    0.391414178028318,
    0.354657936074995,
    0.317901694121672,
    0.28114545216835,
    0.244389210215027,
    0.207632968261704,
    0.170876726308382,
    0.134120484355059,
    0.0973642424017363,
    0.0606080004484136,
    0.0238517584950909,
    -0.0129044834582318,
    -0.0496607254115545,
    -0.0864169673648772,
    -0.1231732093182,
    -0.159929451271522,
    -0.196685693224845,
    -0.233441935178168,
    -0.270198177131491,
    -0.306954419084813,
    -0.343710661038136,
]
Controllability = [
    0.33560523762027,
    0.316976652397814,
    0.298348067175358,
    0.279719481952902,
    0.261090896730446,
    0.24246231150799,
    0.223833726285534,
    0.205205141063078,
    0.186576555840622,
    0.167947970618166,
    0.14931938539571,
    0.130690800173254,
    0.112062214950798,
    0.0934336297283415,
    0.0748050445058855,
    0.0561764592834295,
    0.0375478740609735,
    0.0189192888385175,
    0.000290703616061477,
    -0.0183378816063945,
    -0.0369664668288505,
    -0.0555950520513065,
    -0.0742236372737626,
    -0.0928522224962186,
    -0.111480807718675,
    -0.130109392941131,
    -0.148737978163587,
]
Min = -0.77
Max = -0.39
ylst = np.arange(SSh - 0.02, SSh + 0.02, 0.005)

xline = np.linspace(Min, Max, 2)
yline = np.ones(len(xline)) * SSh

plot.figure(figsize=(10, 7))

plt.Ploty(xlst, Stability, "Stability", colour="green", linestyle="solid")
plt.Ploty(xlst, Stability_NoMargin, "Stability_NoMargin", colour="green", linestyle="dashed")
plt.Ploty(xlst, Controllability, "Controllability", colour="blue", linestyle="solid")
plt.Plotx(Min, ylst, "CG Excursion", colour="red", linestyle="solid")
plt.Plotx(Max, ylst, "", colour="red", linestyle="solid")
plot.plot(xline, yline, color="red", linestyle="solid")
plot.fill_between(xlst, Stability, 1, color="red", alpha=0.1)
plot.fill_between(xlst, 0, Controllability, color="red", alpha=0.1)

plot.xlabel(r"$X_{cg}/MAC$", fontsize=12, loc="right")
plot.ylabel(r"$S_h/S$", fontsize=12, loc="top")

plot.axhline(color="black", lw=0.5)
plot.axvline(color="black", lw=0.5)
plot.legend()
plot.xlim((-1, 0.5))
plot.ylim((0, 0.5))
plot.savefig("Stability_Conv_WithCanard.png")
plot.savefig("Stability_Conv_WithCanard.svg")
plot.show()
