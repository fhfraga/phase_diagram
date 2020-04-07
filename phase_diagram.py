import matplotlib.pyplot as plt
import numpy as np
from scipy import constants
import pint
import pandas as pd

DF = pd.read_csv('data/data.csv')

FORMULAS = DF['Formula']
NAMES = DF['Name']
CAS = DF['CAS_number']

ureg = pint.UnitRegistry()
ureg.setup_matplotlib(True)

gas_constant = constants.gas_constant * ureg.J/(ureg.mol*ureg.K)


class PhaseDiagram:
    def __init__(self, compound):
        search = DF.loc[:, ['Name', 'Formula', 'CAS_number']].isin([compound])
        mask = search.any()

        if mask.any():
            column = mask[mask == True].index[0]
            self.idx = search[column][search[column] == True].index[0]
        else:
            raise ValueError

        # compound identification
        self.name = DF.iloc[self.idx, 0]
        self.formula = DF.iloc[self.idx, 1]
        self.cas = DF.iloc[self.idx, 2]

        # triple point
        self.TP_temperature = DF.iloc[self.idx, DF.columns.get_loc(
            'TP_temperature')] * ureg.kelvin

        self.TP_pressure = DF.iloc[self.idx, DF.columns.get_loc(
            'TP_pressure')] * ureg.pascal

        # critical point
        self.CP_temperature = DF.iloc[self.idx, DF.columns.get_loc(
            'CP_temperature')] * ureg.kelvin
        self.CP_pressure = DF.iloc[self.idx, DF.columns.get_loc(
            'CP_pressure')] * ureg.pascal

        # fusion
        self.H_melt = DF.iloc[self.idx, DF.columns.get_loc(
            'H_melt')] * ureg.kJ / ureg.mol
        self.V_melt = DF.iloc[self.idx, DF.columns.get_loc(
            'V_melt')] * ureg.cc / ureg.mol
        self.V_melt_calc = DF.iloc[self.idx, DF.columns.get_loc(
            'V_melt_calc')] * ureg.cc / ureg.mol

        # vaporization
        self.H_vap = DF.iloc[self.idx, DF.columns.get_loc(
            'H_vap')] * ureg.kJ / ureg.mol
        self.H_vap_boil = DF.iloc[self.idx, DF.columns.get_loc(
            'H_vap_boil')] * ureg.kJ / ureg.mol

        # sublimation
        self.H_sub = DF.iloc[self.idx, DF.columns.get_loc(
            'H_sub')] * ureg.kJ / ureg.mol

        # Antoine equation. Pressure in mmHg and temperature in Celsius
        self.antoine_A = DF.iloc[self.idx, DF.columns.get_loc('A')]
        self.antoine_B = DF.iloc[self.idx, DF.columns.get_loc('B')]
        self.antoine_C = DF.iloc[self.idx, DF.columns.get_loc('C')]
        self.antoine_Tmin = DF.iloc[self.idx, DF.columns.get_loc('Tmin')]
        self.antoine_Tmax = DF.iloc[self.idx, DF.columns.get_loc('Tmax')]

    def clapeyron_sl(self, temp_range=5):
        # P(T) = P' + (H_melt / V_melt) ln(T / T')
        if np.isnan(self.V_melt.magnitude):
            V_melt = self.V_melt_calc
        else:
            V_melt = self.V_melt

        if V_melt > 0:
            temp_range = -temp_range

        T_arr = np.linspace(self.TP_temperature.magnitude -
                            temp_range, self.TP_temperature.magnitude, 100) * ureg.K
        cte = self.H_melt / V_melt
        P_arr = self.TP_pressure + cte * np.log(T_arr / self.TP_temperature)
        return T_arr, P_arr

    def clapeyron_sv(self, temp_range=60):
        # P(T) = P' exp[ (H_sub / R) (1 / T' - 1 / T) ]
        T_arr = np.linspace(self.TP_temperature.magnitude -
                            temp_range, self.TP_temperature.magnitude, 100) * ureg.K
        cte = self.H_sub / gas_constant
        P_arr = self.TP_pressure * \
            np.exp(cte * (1/self.TP_temperature - 1/T_arr))
        return T_arr, P_arr

    def clapeyron_lv(self):
        # P(T) = P' exp[ (H_vap / R) (1 / T' - 1 / T) ]
        T_arr = np.linspace(self.TP_temperature.magnitude,
                            self.CP_temperature.magnitude, 100) * ureg.K

        if np.isnan(self.H_vap_boil.magnitude):
            H_vap = self.H_vap
        else:
            H_vap = self.H_vap_boil

        cte = H_vap / gas_constant
        P_arr = self.TP_pressure * \
            np.exp(cte * (1/self.TP_temperature - 1/T_arr))
        return T_arr, P_arr

    def antoine_lv(self):
        # log10(P) = A - (B / (C + T))
        T_arr = np.linspace(self.TP_temperature.magnitude,
                            self.CP_temperature.magnitude, 100) * ureg.K

        A = self.antoine_A + np.log10(101325/760)
        B = self.antoine_B
        C = self.antoine_C - 273.15
        Tmin = self.antoine_Tmin + 273.15
        Tmax = self.antoine_Tmax + 273.15

        right_side = A - (B / (C + T_arr.magnitude))

        P_arr = 10**right_side * ureg.Pa

        return T_arr, P_arr, A, B, C, Tmin, Tmax
