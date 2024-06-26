# Sources for flight performance evaluation:
#
# ruijgrok
# = "element of airplane performance" (physical book)
#
# gudmundsson:
# = https://www.sciencedirect.com/book/9780123973085/general-aviation-aircraft-design
#
# nicolai
# = https://app-knovel-com.tudelft.idm.oclc.org/kn/resources/kpFAADVIA3/toc?cid=kpFAADVIA3
#


from HumanAir.FlightPerformance import aircraft
import matplotlib.pyplot as plt
import numpy as np

acf = aircraft.Aircraft()


CL_list = np.linspace(0, 3, 50)
CD_list = np.zeros(len(CL_list))
CD_TO = np.zeros(len(CL_list))
CD_ld = np.zeros(len(CL_list))

for i, CL in enumerate(CL_list):
    CD_list[i] = acf.CD(CL)
    CD_TO[i] = acf.CD(CL, flaps="TO")
    CD_ld[i] = acf.CD(CL, flaps="land", gear="down")

plt.figure(figsize=(10, 7))
plt.plot(CL_list, CD_list, label="CL, clean")
plt.plot(CL_list, CD_TO, label="CL, takeoff flap")
plt.plot(CL_list, CD_ld, label="CL, landing flap, gear down")
plt.scatter(0.880069678, 0.056409609, label="CL/CD CFD lg down")
plt.scatter(0.873324273, 0.04590758, label="CL/CD CFD lg up")
plt.xlabel("C_L")
plt.ylabel("C_D")
plt.ylim(0, 0.25)
plt.legend()
plt.show()

print(acf.LDmax)
