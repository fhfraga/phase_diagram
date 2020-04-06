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

        # sublimation
        self.H_sub = DF.iloc[self.idx, DF.columns.get_loc(
            'H_sub')] * ureg.kJ / ureg.mol

        # Antoine equation. Pressure in mmHg and temperature in Celsius
        self.antoine_A = DF.iloc[self.idx, DF.columns.get_loc('A')]
        self.antoine_B = DF.iloc[self.idx, DF.columns.get_loc('B')]
        self.antoine_C = DF.iloc[self.idx, DF.columns.get_loc('C')]
        self.antoine_Tmin = DF.iloc[self.idx, DF.columns.get_loc('Tmin')]
        self.antoine_Tmax = DF.iloc[self.idx, DF.columns.get_loc('Tmax')]

    def clapeyron_sl(self):
        pass

    def clapeyron_sv(self):
        pass

    def clapeyron_lv(self):
        pass

    def antoine_lv(self):
        pass
