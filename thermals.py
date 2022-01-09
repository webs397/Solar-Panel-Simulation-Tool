from scipy.constants import sigma
from math import cos


def incoming_heat_flow(a_s, plate_area, irradiation_global):
    """Calculates the incoming heat flow from global irradiation, output is in Watts"""
    ihf = a_s * plate_area * irradiation_global
    return ihf


def heat_exchange_coefficient(Nusselt-Number,lamda,lenght):
    """Calculates the heat exchange coefficient (Wärmeübergangskoeffizient)"""
    hec = (Nusselt-Number * lamda)/ lenght 

    return


def convective_heat_flow(heat_exchange_coeff, plate_area, temperature_air):
    """Calculates the convective heat flow across the plate"""
    chf = heat_exchange_coeff * plate_area * temperature_air
    return chf


def sky_temperature(air_temp, dew_temp, standard_time):
    """Calculates the sky temperature as referred to in Solar Engineering, Duffie and Beckman, Page 148
        if clouds are included in the calculation then the sky temperature must be higher than a cloudless sky"""
    sky_temp = air_temp * pow(0.711 + 0.0056 * dew_temp + 0.000073 * pow(dew_temp, 2) + 0.013 * cos(15 * standard_time),
                              1 / 4)
    return sky_temp


def radiation_heat_flow(em_deg, area, plate_temp, sky_temp):    # Emission degree of Plate is user defined or given by chosen material
    """Calculates the heat flow resulting from radiation exchange with the sky, WUE P.181"""
    rhf = em_deg * sigma * area * (pow(plate_temp, 4) - pow(sky_temp, 4))
    return rhf


def solar_heat_flow(solar_absorbtion_coeff, area,
                    irradiation_g):  # Solar absorbtion coefficient is a user input or material default
    """Calculates the heat flow resulting from solar radiation"""
    shf = solar_absorbtion_coeff * area * irradiation_g
    return shf
