import matplotlib.pyplot as plt
import numpy as np
from scipy import constants
import pint
import pandas as pd
import re

DF = pd.read_csv('data/data.csv')

FORMULAS = DF['Formula']
NAMES = DF['Name']
CAS = DF['CAS_number']

ureg = pint.UnitRegistry()
ureg.setup_matplotlib(True)

gas_constant = constants.gas_constant * ureg.J/(ureg.mol*ureg.K)


class PhaseDiagram:
    def __init__(self, compound):
        """Object initialization

        Parameters
        ----------
        compound : string
            A valid formula, name or CAS number

        Raises
        ------
        ValueError
            Not a valid formula, name or CAS number. Compound not in the
            available data.
        """
        search = DF.loc[:, ['Name', 'Formula', 'CAS_number']].isin([compound])
        mask = search.any()

        if mask.any():
            column = mask[mask == True].index[0]
            self.idx = search[column][search[column] == True].index[0]
        else:
            raise ValueError('Not a valid compound.')

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
        """Clausius-Clapeyron solid-liquid line data

        Parameters
        ----------
        temp_range : int, optional
            Temperature range around the triple point, by default 5

        Returns
        -------
        tuple
            Tuple of arrays (temperature, pressure)
        """
        # P(T) = P' + (H_melt / V_melt) ln(T / T') where T' is TP_temperature
        if np.isnan(self.V_melt.magnitude):
            V_melt = self.V_melt_calc
        else:
            V_melt = self.V_melt

        if V_melt > 0:
            temp_range = -temp_range

        T_arr = np.linspace(self.TP_temperature.magnitude,
                            self.TP_temperature.magnitude-temp_range,
                            100) * ureg.K
        cte = self.H_melt / V_melt
        P_arr = self.TP_pressure + cte * np.log(T_arr / self.TP_temperature)
        return T_arr, P_arr

    def clapeyron_sv(self, temp_range=60):
        """Clausius-Clapeyron solid-vapor line data

        Parameters
        ----------
        temp_range : int, optional
            Temperature range around the triple point, by default 60

        Returns
        -------
        tuple
            Tuple of arrays (temperature, pressure)
        """
        # P(T) = P' exp[ (H_sub / R) (1 / T' - 1 / T) ] where T' is TP_temperature
        T_arr = np.linspace(self.TP_temperature.magnitude - temp_range,
                            self.TP_temperature.magnitude,
                            100) * ureg.K
        cte = self.H_sub / gas_constant
        P_arr = self.TP_pressure * \
            np.exp(cte * (1/self.TP_temperature - 1/T_arr))
        return T_arr, P_arr

    def clapeyron_lv(self):
        """Clausius-Clapeyron liquid-vapor line data

        Returns
        -------
        tuple
            Tuple of arrays (temperature, pressure)
        """
        # P(T) = P' exp[ (H_vap / R) (1 / T' - 1 / T) ] where T' is TP_temperature
        T_arr = np.linspace(self.TP_temperature.magnitude,
                            self.CP_temperature.magnitude,
                            100) * ureg.K

        if np.isnan(self.H_vap_boil.magnitude):
            H_vap = self.H_vap
        else:
            H_vap = self.H_vap_boil

        cte = H_vap / gas_constant
        P_arr = self.TP_pressure * \
            np.exp(cte * (1/self.TP_temperature - 1/T_arr))
        return T_arr, P_arr

    def antoine_lv(self):
        """Antoine liquid-vapor line data

        Returns
        -------
        tuple
            A, B and C for SI units. Temperature range (Tmin and Tmax) in Kelvin
            (temperature array, pressure array, A, B, C, Tmin, Tmax)
        """
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

    def format_formula(self):
        """ Display chemical formulas in a proper way

        Returns
        -------
        string
            LaTeX code to display chemical formulas in a proper way
        """
        label_formula = re.sub("([0-9])", "_\\1", self.formula)
        label_formula = '$\mathregular{'+label_formula+'}$'
        return label_formula

    def _plot_params(self, ax=None):
        """Internal function for plot parameters.

        Parameters
        ----------
        ax : Matplotlib axes, optional
            axes where the graph will be plotted, by default None
        """
        linewidth = 2
        size = 12

        # grid and ticks settings
        ax.minorticks_on()
        ax.grid(b=True, which='major', linestyle='--',
                linewidth=linewidth - 0.5)
        ax.grid(b=True, which='minor', axis='both',
                linestyle=':', linewidth=linewidth - 1)
        ax.tick_params(which='both', labelsize=size+2)
        ax.tick_params(which='major', length=6, axis='both')
        ax.tick_params(which='minor', length=3, axis='both')

        # labels and size
        ax.xaxis.label.set_size(size+4)
        ax.yaxis.label.set_size(size+4)
        # ax.title.set_fontsize(size+6)  # not working, don't know why...

        return

    def plot(self, parts=(1, 1, 0, 1), size=(10, 8), ax=None, T_unit='K',
             P_unit='Pa', scale_log=True, legend=False, title=True,
             title_text=''):
        """Plot function

        Parameters
        ----------
        parts : tuple, optional
            which lines will be plotted, by default (1, 1, 0, 1)
            By default, the solid-liquid, solid-vapor and liquid-vapor from
            Antoine equation lines are plotted. This can be changed with 0 and
            1's in a tuple`. 0 means turn off and 1 means turn on. The order in
            the tuple is:
            (solid-liquid Clausius-Clapeyron, solid-vapor Clausius-Clapeyron,
            liquid-vapor Clausius-Clapeyron, liquid-vapor Antoine)
        size : tuple, optional
            plot size, by default (10, 8)
        ax : Matplotlib axes, optional
            axes where the graph will be plotted, by default None
        T_unit : str, optional
            temperature unit, by default 'K'
        P_unit : str, optional
            pressure unit, by default 'Pa'
        scale_log : bool, optional
            logarithmic scale, by default True
        legend : bool, optional
            If a legend will be shown, by default False
        title : bool, optional
            If the plot will have a title, by default True
        title_text : str, optional
            Title text, by default ''

        Returns
        -------
        Matplotlib axes
            axes where the graph will be plotted
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=size, facecolor=(1.0, 1.0, 1.0))

        self._plot_params(ax)

        linewidth = 3.0

        if parts[0] == 1:
            T_clapeyron_sl, P_clapeyron_sl = self.clapeyron_sl()

            # ax.plot(T_clapeyron_sl.to(T_unit), P_clapeyron_sl.to(P_unit))
            # in order to avoid long SL lines, limit the pressure values to
            # those lower than self.CP_pressure
            P_clapeyron_sl = P_clapeyron_sl[P_clapeyron_sl < self.CP_pressure]
            ax.plot(T_clapeyron_sl[:len(P_clapeyron_sl)].to(T_unit),
                    P_clapeyron_sl.to(P_unit),
                    'k-', label='SL boundary', linewidth=linewidth)

        if parts[1] == 1:
            T_clapeyron_sv, P_clapeyron_sv = self.clapeyron_sv()
            ax.plot(T_clapeyron_sv.to(T_unit),
                    P_clapeyron_sv.to(P_unit),
                    'b-', label='SV boundary', linewidth=linewidth)

        if parts[2] == 1:
            T_clapeyron_lv, P_clapeyron_lv = self.clapeyron_lv()
            ax.plot(T_clapeyron_lv.to(T_unit),
                    P_clapeyron_lv.to(P_unit),
                    'g--', label='LV boundary', linewidth=linewidth)

        if parts[3] == 1:
            T_antoine_lv, P_antoine_lv, *_ = self.antoine_lv()
            ax.plot(T_antoine_lv.to(T_unit),
                    P_antoine_lv.to(P_unit),
                    'r-', label='LV boundary - Antoine', linewidth=linewidth)

        if parts[2] == 1 or parts[3] == 1:
            ax.scatter(self.CP_temperature, self.CP_pressure,
                       s=100, label='Critical Point',
                       facecolors='orange', edgecolors='orange', zorder=3)

        ax.scatter(self.TP_temperature, self.TP_pressure,
                   s=100, label='Triple Point',
                   facecolors='m', edgecolors='m', zorder=3)

        if scale_log:
            ax.set_yscale('log')
            ax.set_ylabel('log(Pressure / {:~P})'.format(ureg(P_unit).units))
        else:
            # setting the y-axis to scientific notation and
            # getting the order of magnitude
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            ax.yaxis.major.formatter._useMathText = True
            ax.figure.canvas.draw()  # Update the text
            order_magnitude = ax.yaxis.get_offset_text().get_text().replace('\\times', '')
            ax.yaxis.offsetText.set_visible(False)
            ax.set_ylabel('Pressure / ' + order_magnitude +
                          ' {:~P}'.format(ureg(P_unit).units))

        ax.set_xlabel('Temperature / {:~P}'.format(ureg(T_unit).units))

        if legend:
            ax.legend(loc='best', fontsize=14,
                      title=self.format_formula(), title_fontsize=14)

        if not title:
            pass
        elif title_text == '':
            ax.set_title('Calculated phase diagram - ' + self.format_formula(),
                         fontsize=18)
        else:
            ax.set_title(title_text, fontsize=18)

        return ax
