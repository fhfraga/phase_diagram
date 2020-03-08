import matplotlib.pyplot as plt
import numpy as np
from scipy import constants


def phase_diagram(temperature_triple, pressure_triple, temperature_crit,
                  pressure_crit, enthalpy_melt, volume_melt, enthalpy_sub,
                  enthalpy_vap, A_ant_init, A_ant_fin, B_ant_int, B_ant_fin,
                  C_ant_init, C_ant_fin, gas_constant=constants.gas_constant):
    """[summary]

    Parameters
    temperature_triple : float
        It's the temperature where the three phases of the substance
        (solid, liquid and gas) coexist in thermodynamic equilibrium.
        Dimension: temperature
    pressure_triple : float
        It's the pressure where the three phases of the substance
        (solid, liquid and gas) coexist in thermodynamic equilibrium
        Dimension: pressure
    temperature_crit : float
        It's the final temperature of a phase equilibrium curve
        Dimension: temperature
    pressure_crit : float
        It's the final temperature and pressure of a phase equilibrium curve
        Dimension: pressure
    enthalpy_melt : float
        Melting enthalpy of a substance
        Dimension: energy
    volume_melt : float
        Volume of melting of a substance
        Dimension: Volume
    enthalpy_sub : float
        Sublimation enthalpy of a substance
        Dimension: energia
    enthalpy_vap : float
        Boiling enthalpy of the substance
        Dimension: energia
    A_ant_init : float
        Antoine constant initial number for triple point temperature
        Dimensionless
    A_ant_fin : float
        Final number of Antoine constant for critical point temperature
        Dimensionless
    B_ant_int : float
        Antoine constant initial number for triple point temperature
        Dimensionless
    B_ant_fin : float
        Final number of Antoine constant for critical point temperature
        Dimensionless
    C_ant_init : float
        Antoine constant initial number for triple point temperature
        Dimensionless
    C_ant_fin : float
        Final number of Antoine constant for critical point temperature
        Dimensionless
    gas_constant : float, optional
       Ideal gas constant , by default constants.gas_constant
        Dimension: Expressed in units of energy per temperature increment per
        mole
    """

    name_substance = name_of_substance
    # solid-liquid
    temperature_melt = np.linspace(
        temperature_triple - 5.0, temperature_triple, 100)
    pressure_melt = pressure_triple + \
        ((enthalpy_melt / volume_melt) * 10e5) * \
        np.log(temperature_melt / temperature_triple)

    # solid-gas
    temperature_sub = np.linspace(
        temperature_triple - 60.0, temperature_triple, 100)
    pressure_sub = pressure_triple * \
        np.exp((enthalpy_sub / gas_constant) *
               ((1/temperature_triple) - (1/temperature_sub)))

    # liquid-gas (Equation Clausius-Clapeyron)
    temperature_vap = np.linspace(temperature_triple, temperature_crit, 100)
    pressure_vap = pressure_triple * \
        np.exp((enthalpy_vap/gas_constant) *
               ((1/temperature_triple) - (1/temperature_vap)))

    # líquid-gas (Equation Antoine)
    temperature_vap_antoine = np.linspace(
        temperature_triple, temperature_crit, 100)
    A = np.linspace(A_ant_init, A_ant_fin, 100)
    B = np.linspace(B_ant_int, B_ant_fin, 100)
    C = np.linspace(C_ant_init, C_ant_fin, 100)

    # 10e5 is for conversion from bar to pascal (data in NIST in bar
    pressure_vap_antoine = (
        10 ** (A - (B/(C + temperature_vap_antoine)))) * 10 ** 5

    figure = plt.figure(figsize=(12, 10))
    plt.semilogy(temperature_melt, pressure_melt, linewidth=3.0)
    plt.semilogy(temperature_sub, pressure_sub, color='orange', linewidth=3.0)
    plt.semilogy(temperature_vap, pressure_vap,
                 '--', color='red', linewidth=3.0)
    plt.semilogy(temperature_vap_antoine, pressure_vap_antoine,
                 color='green', linewidth=3.0)
    plt.scatter(temperature_crit, pressure_crit, s=100,
                facecolors='none', edgecolors='k')
    plt.scatter(temperature_triple, pressure_triple,
                s=100, facecolors='none', edgecolors='k')
    plt.xlabel('Temperature [K]')
    plt.ylabel('Pressure [Pa]')
    plt.title('Phase diagram of ' + name_substance)
    plt.savefig('Phase diagram of ' + name_substance)
    return plt.show()


if __name__ == "__main__":
    print()
    print()
    print('#'*78)
    print('# {0:^74} #'.format('Phase diagram'))
    print('#'*78)
    print()

    print('provide the following data.' +
          ' Attention to the units do not forget to convert if necessary.')
    print()

name_of_substance = str(input("name of substance: "))
temperature_triple = eval(
    input('Triple point temperature of substance \ Kelvin: '))
pressure_triple = eval(input('triple point pressure of substance \ Pascal: '))
temperature_crit = eval(
    input('Critical point temperature of substance \ Kelvin: '))
pressure_crit = eval(input('Critical point pressure of substance \ Pascal: '))
enthalpy_melt = eval(input('Melting enthalpy of substance \ Joule/mol: '))
volume_melt = eval(input('Melting volume of substance \ cm3/mol: '))
enthalpy_sub = eval(input('Sublimation enthalpy of substance \ Joule/mol: '))
enthalpy_vap = eval(input('Boiling enthalpy of substance \ Joule/mol: '))
A_ant_init = eval(
    input('Antoine equation parameters - A of the triple point: '))
A_ant_fin = eval(
    input('Antoine equation parameter -  A of the critical point: '))
B_ant_int = eval(
    input('Antoine equation parameters - B of the triple point: '))
B_ant_fin = eval(
    input('Antoine equation parameter -  B of the critical point: '))
C_ant_init = eval(
    input('Antoine equation parameters - C of the triple point: '))
C_ant_fin = eval(
    input('Antoine equation parameter -  C of the critical point: '))


phase_diagram(temperature_triple, pressure_triple, temperature_crit,
              pressure_crit, enthalpy_melt, volume_melt, enthalpy_sub,
              enthalpy_vap, A_ant_init, A_ant_fin, B_ant_int, B_ant_fin,
              C_ant_init, C_ant_fin)