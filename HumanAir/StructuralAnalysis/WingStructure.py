import numpy as np
import os
import sys
import copy

# import pandas as pd
# import matplotlib.pyplot as plt
# import re

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from HumanAir.aircraft_data import aircraft_data, airfoil_shape

# from Functions import chord, import_data2


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def interp(x, x_given, y_given):
    i = find_nearest(x_given, x)
    if i == 0:
        return y_given[0]
    elif i == len(x_given):
        return y_given[-1]
    else:
        x1, x2 = x_given[i - 1], x_given[i]
        y1, y2 = y_given[i - 1], y_given[i]
        return y1 + (y2 - y1) * (x - x1) / (x2 - x1)


class WingStructure:
    # TODO: Make it so we can choose between wing and HT
    def __init__(self, ac_data, airfoil_data, nodes=501):
        self.airfoil_data = airfoil_data
        # self.file_path_y = file_path_y  # file path to span
        self.Sw = ac_data["Aero"]["S_Wing"]
        self.taper_ratio = ac_data["Aero"]["Taper_Wing"]
        # self.AoA = AoA  # [deg]
        self.nodes = nodes  # number of nodes
        self.t1_spar = ac_data["Geometry"]["t_spar_tip"]  # [m] thickness at the tip
        self.t2_spar = ac_data["Geometry"]["t_spar_root"]  # [m] thickness at the root
        self.t_skin = ac_data["Geometry"]["t_skin_wing"]  # [m] thickness of skin
        self.spar_pos = np.array(ac_data["Geometry"]["spar_pos"])  # position of spars
        self.stringer_area = ac_data["Geometry"]["wing_stringer_area_m2"]  # [m2] area of stringers
        self.cr = ac_data["Aero"]["c_root_wing"]
        self.ct = self.cr * self.taper_ratio
        self.b = ac_data["Aero"]["b_Wing"]
        self.ypts = np.linspace(-self.b / 2, self.b / 2, self.nodes)
        self.stringer_number = ac_data["Geometry"]["wing_stringer_number"]
        self.stringer_sections = ac_data["Geometry"]["wing_stringer_sections"]
        self.material_rho = ac_data["Materials"][ac_data["Geometry"]["wingbox_material"]]["rho"]
        self.total_structural_weight = ac_data["Structures"]["structural_wing_weight"]
        self.material_E = ac_data["Materials"][ac_data["Geometry"]["wingbox_material"]]["E"]
        self.material_G = ac_data["Materials"][ac_data["Geometry"]["wingbox_material"]]["G"]
        self.w_fuselage = ac_data["Geometry"]["fus_width_m"]

        self.chord_distribution = self.calculate_chord_distribution()
        self.spars = self.chord_distribution.reshape(-1, 1) * self.spar_pos

        # Assume spar thickness is uniformly decreasing from root to tip
        self.t_spar_dist = self.calc_spar_dist()

        self.stringer_dist = self.calc_stringer_dist()

        self.airfoil_division = self.calc_airfoil_division()
        self.y_up, self.y_down = self.y_updown()

        self.hmax_dist = self.calc_hmax_dist()

        self.enclosed_area_dist = self.calc_enclosed_area_dist()

    # def import_data(self):
    #     df = pd.read_csv(self.file_path, sep="\s+", header=None, names=["x", "y"], skiprows=1)
    #     return df

    # def import_data2(self):
    #     data = {}
    #     with open(self.file_path_y, "r") as file:
    #         lines = file.readlines()
    #         angle_file = lines[0]
    #         angles = re.findall(r"VLM1 -\s*-?\d+\.?\d*", angle_file)
    #         angles = [float(re.search(r"VLM1 -\s*(-?\d+\.?\d*)", angle).group(1)) for angle in angles]

    #         for line in lines[1:]:
    #             values = line.split()
    #             for i in range(0, len(values), 2):
    #                 angle_index = i // 2
    #                 angle = angles[angle_index]
    #                 y_positions = float(values[i])
    #                 lift_distribution = float(values[i + 1])

    #                 if angle not in data:
    #                     data[angle] = {"y_span": [], "coefficient": []}

    #                 data[angle]["y_span"].append(y_positions)
    #                 data[angle]["coefficient"].append(lift_distribution)
    #     return data

    def calc_spar_dist(self):
        return self.t1_spar + (self.t2_spar - self.t1_spar) / (self.cr - self.ct) * (self.chord_distribution - self.ct)

    def calc_hmax_dist(self):
        hmax_airfoil = np.max(self.airfoil_division[0][:, 1])
        return hmax_airfoil * self.chord_distribution

    def calc_stringer_dist(self):
        if len(self.stringer_sections) != len(self.stringer_number):
            raise ValueError("The number of stringer sections should be equal to the number of stringers per section")

        stringers = np.ones(self.nodes // 2)
        stringer_section_ends = np.rint(np.cumsum(self.stringer_sections) * (len(stringers))).astype(int)
        stringer_section_starts = np.insert(stringer_section_ends[:-1], 0, 0)
        for i in range(len(self.stringer_sections)):
            start = stringer_section_starts[i]
            end = stringer_section_ends[i]
            # print(start, end)
            stringers[start:end] = self.stringer_number[i]

        middle = [] if self.nodes % 2 == 0 else [self.stringer_number[0]]
        stringers = np.concatenate([np.flip(stringers), middle, stringers])

        return stringers

    def calc_airfoil_division(self):
        df = copy.deepcopy(self.airfoil_data)
        for i in range(len(df)):
            if df["x"][i] == 0:
                idx = i
        df_1 = np.array(df.iloc[:idx, :].values.tolist())
        df_2 = np.array(df.iloc[idx:, :].values.tolist())
        return df_1, df_2

    def calculate_chord_distribution(self):
        chord_length = (
            2
            * self.Sw
            / (1 + self.taper_ratio)
            / self.b
            * (1 - ((1 - self.taper_ratio) / self.b * np.abs(2 * self.ypts)))
        )
        return chord_length

    # def chord(self):
    #     Cl_DATA = import_data2(self.file_path_y)
    #     b = Cl_DATA[self.AoA]['y_span'][-1] * 2
    #     y = np.linspace(Cl_DATA[AoA]['y_span'][0], Cl_DATA[AoA]['y_span'][-1], self.n)

    # def spars(self):
    #     return self.chord_distribution * self.spar_pos

    def y_spar(self, df):
        x_spar = self.spars
        y_spar = []
        for i in range(len(x_spar)):
            for j in range(len(x_spar[i])):
                y_spar.append(
                    interp(x_spar[i][j] / self.chord_distribution[i], df[:, 0], df[:, 1]) * self.chord_distribution[i]
                )
        y_spar = np.array(y_spar)
        return y_spar.reshape((len(x_spar), 2))

    def y_updown(self):
        df_up, df_down = self.airfoil_division
        y_up = self.y_spar(df_up)
        y_down = self.y_spar(df_down)
        return y_up, y_down

    def diff(self, data):
        dx = []
        for i in range(len(data)):
            for j in range(len(data[i])):
                if j == len(data[i]) - 1:
                    continue
                else:
                    dxx = data[i][j + 1] - data[i][j]
                    dx.append(dxx)
        return dx

    def y(self, data):
        dx = []
        for i in range(len(data)):
            for j in range(len(data[i])):
                if j == len(data[i]) - 1:
                    continue
                else:
                    dxx = (data[i][j + 1] + data[i][j]) / 2
                    dx.append(dxx)
        return np.array(dx)

    def d_s1s2(self):
        dx = self.diff(self.spars)
        dy_up = self.diff(self.y_up)
        dy_down = self.diff(self.y_down)
        l_box_up = []
        l_box_down = []
        for i in range(len(dx)):
            l_box_up.append(np.sqrt(dx[i] ** 2 + dy_up[i] ** 2))
            l_box_down.append(np.sqrt(dx[i] ** 2 + dy_down[i] ** 2))
        return np.array(l_box_up), np.array(l_box_down)

    def h_s1s2(self):
        h_mid = np.empty((len(self.y_up), len(self.y_up.T)))
        h_emp = np.empty((len(self.y_up), len(self.y_up.T)))
        for i in range(len(self.y_up)):
            for j in range(len(self.y_up[i])):
                h_mid[i, j] = (self.y_up[i][j] + self.y_down[i][j]) / 2
                h_emp[i, j] = self.y_up[i][j] - self.y_down[i][j]
        return h_mid, h_emp

    def centroid(self):
        y_up, y_down = self.y_updown()
        l_box_up, l_box_down = self.d_s1s2()
        h_mid, h_s1s2 = self.h_s1s2()

        # upper skin
        A_uskin = l_box_up.reshape((len(l_box_up), 1)) * self.t_skin
        y_uskin = self.y(y_up).reshape((len(y_up), 1))

        # lower skin
        A_lskin = l_box_down.reshape((len(l_box_down), 1)) * self.t_skin
        y_lskin = self.y(y_down).reshape((len(y_down), 1))

        # left spar
        A_lspar = self.t_spar_dist * h_s1s2[:, 0].reshape((len(h_s1s2), 1))
        y_lspar = h_mid[:, 0].reshape((len(h_mid), 1))

        # right spar
        A_rspar = self.t_spar_dist * h_s1s2[:, 1].reshape((len(h_s1s2), 1))
        y_rspar = h_mid[:, 1].reshape((len(h_mid), 1))

        return (A_uskin * y_uskin + A_lskin * y_lskin + A_lspar * y_lspar + A_rspar * y_rspar) / (
            A_uskin + A_lskin + A_lspar + A_rspar
        )

    def Ixx(self):
        _, h_s1s2 = self.h_s1s2()
        l_box_up, l_box_down = self.d_s1s2()
        h_frontspar = h_s1s2[:, 0].flatten()
        h_rearspar = h_s1s2[:, 1].flatten()
        htot = h_frontspar + h_rearspar
        h_avemax = htot / 4
        I_stringer = self.stringer_area * self.stringer_dist * h_avemax**2
        I_spar = 1 / 12 * self.t_spar_dist * (h_frontspar**3 + h_rearspar**3)
        I_skin = self.t_skin * (l_box_down + l_box_up) * h_avemax**2
        return I_stringer + I_spar + I_skin

    def shear_centre(self):
        Ixx = self.Ixx()
        V = 1

        _, h_s1s2 = self.h_s1s2()
        h_frontspar = h_s1s2[:, 0]
        h_rearspar = h_s1s2[:, 1]

        spar_distance = self.spars[:, 1] - self.spars[:, 0]

        l_ab = h_rearspar
        l_dc = h_frontspar
        l_ad = l_bc = np.sqrt((spar_distance) ** 2 + ((h_frontspar - h_rearspar) / 2) ** 2)
        beta = np.arctan((h_frontspar - h_rearspar) / (2 * spar_distance))  # in rad

        qab_part = V / Ixx * (-(l_bc**3) * np.sin(beta) / 6 + l_dc * l_bc**2 / 4)

        qad_part = (
            V
            / Ixx
            * (
                l_ab**3 / 4
                - l_ab**3 / 6
                + self.t_skin / self.t_spar_dist * (l_dc * l_ab * l_bc / 2 - l_ab * l_ad**2 * np.sin(beta) / 2)
            )
        )
        qdc_part = (
            V
            / Ixx
            * (-(l_ad**2) * l_ab / 4 - l_ad**3 / 6 * np.sin(beta) + l_dc * l_ad**2 / 2 - l_ad**3 * np.sin(beta) / 2)
        )
        qbc_part = (
            V
            / Ixx
            * (
                -(l_dc**3) / 4
                + l_dc**3 / 6
                + self.t_skin / self.t_spar_dist * (-(l_ad**2) * l_dc * np.sin(beta) + l_dc * (l_dc - l_ab) * l_ad / 2)
            )
        )

        qs0 = -(
            (qab_part + qad_part + qdc_part + qbc_part)
            / (l_ab / self.t_spar_dist + l_ad / self.t_skin + l_dc / self.t_spar_dist + l_bc / self.t_skin)
        )

        A_tra = (l_ab + l_dc) / 2 * (self.spars[:, 1] - self.spars[:, 0])

        # Shear centre
        eta = (
            2 * A_tra * qs0
            + (self.spars[:, 1] - self.spars[:, 0])
            * V
            / Ixx
            * (
                self.t_spar_dist * (l_ab**3 / 4 - l_ab**3 / 6)
                + self.t_skin * (l_dc * l_ad * l_ab / 2 - l_ad**2 * np.sin(beta) * l_ab / 2)
            )
            + V
            / Ixx
            * l_dc
            * self.t_skin
            * np.cos(beta)
            * (-l_ab * l_ad**2 / 4 - l_ad**3 / 6 * np.sin(beta) + l_dc * l_ad**2 / 2 - l_bc**3 * np.sin(beta) / 2)
        ) / V

        return eta + self.spars[:, 0]  # distance from LE

    def torsional_constant(self):
        t_spar = self.t_spar_dist
        t_skin = self.t_skin
        enclosed_area = self.enclosed_area_dist

        _, h_s1s2 = self.h_s1s2()
        l_box_up, l_box_down = self.d_s1s2()

        h_frontspar = h_s1s2[:, 0]
        h_rearspar = h_s1s2[:, 1]
        htot = (h_frontspar + h_rearspar).flatten()
        w_top = l_box_up.flatten()
        w_bottom = l_box_down.flatten()
        wtot = w_bottom + w_top

        J = 4 * (enclosed_area**2) / ((htot / t_spar) + (wtot / t_skin))
        return J

    def weight_dist(self):
        _, h_s1s2 = self.h_s1s2()
        l_box_up, l_box_down = self.d_s1s2()
        h_frontspar = h_s1s2[:, 0].flatten()
        h_rearspar = h_s1s2[:, 1].flatten()
        htot = h_frontspar + h_rearspar
        w_top = l_box_up.flatten()  # width of top "straight" skin, as a function of y
        w_bottom = l_box_down.flatten()  # width of bottom "straight" skin, as a function of y
        wtot = w_bottom + w_top

        W_spar = self.material_rho * self.t_spar_dist * htot
        W_skin = self.material_rho * self.t_skin * wtot
        W_stringers = self.material_rho * self.stringer_area * self.stringer_dist
        weight = W_skin + W_spar + W_stringers
        return weight  # [N/m]

    def calc_enclosed_area_dist(self):
        _, h_s1s2 = self.h_s1s2()
        h_frontspar = h_s1s2[:, 0].flatten()
        h_rearspar = h_s1s2[:, 1].flatten()
        return (h_frontspar + h_rearspar) / 2 * (self.spar_pos[1] - self.spar_pos[0]) * self.chord_distribution


if __name__ == "__main__":  # pragma: no cover
    # INPUT
    # Sw = 34.56  # [m2]
    # taper_ratio = 0.4
    # AoA = -6  # [deg]
    n = 500
    # t1_spar = 0.010  # [m] thickness at the tip
    # t2_spar = 0.025  # [m] thickness at the root
    # t_skin = 0.007  # [m] thickness of skin
    # file_path = "HumanAir/WingBox/airfoil.txt"
    # file_path_y = "HumanAir\WingBox\Cl_DATA.txt"
    # A_str = 0.02
    # Cr = 2.5  # [m] root chord length
    # b = 19.93
    # x_pos = np.array([0.15, 0.5])

    wing_structure_data = WingStructure(aircraft_data, airfoil_shape, nodes=501)

    df = wing_structure_data.airfoil_data
    chord1 = wing_structure_data.chord_distribution
    y = wing_structure_data.ypts
    # plt.plot(df['x']*chord1[1], df['y']*chord1[1])
    # plt.axis('equal')
    # plt.show()
    idx = find_nearest(y, 3.986)
    # print(y[idx])

    # print(torisonal_stiffness.spars())
    h_mid, h_s1s2 = wing_structure_data.h_s1s2()
    # print(h_s1s2[idx])

    # plt.plot(df_down[:,0], df_down[:,1])
    # plt.plot(df_up[:,0], df_up[:,1])
    # plt.plot(x_spar[0]/chord_length[0], y_spar_down[0]/chord_length[0])
    # plt.plot(x_spar[0]/chord_length[0], y_spar_up[0]/chord_length[0])
    # plt.axis('equal')
    # plt.show()
