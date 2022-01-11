from scipy.constants import sigma
from math import cos, sqrt
from numpy import log

LAMBDA = 0.025958010
dynamic_viscosity_g = 1.8264E-05
specific_isobaric_heat_capacity_g = 1.005470208
Temperature_surface = 294,2 


def incoming_heat_flow(a_s, plate_area, irradiation_global):
    """Calculates the incoming heat flow from global irradiation, output is in Watts"""
    ihf = a_s * plate_area * irradiation_global
    return ihf


def Prandtl(dynamic_viscosity_g, specific_isobaric_heat_capacity_g, LAMBDA):
    Pr = (dynamic_viscosity_g * specific_isobaric_heat_capacity_g) / LAMBDA
    return Pr


def Reynold_m(wind_velocity, length, dynamic_viscosity_g, density_g):
    """ mittlere Reynolds-Zahl - (8.17) S.231 """
    Re_m = (wind_velocity * length * density_g) / dynamic_viscosity_g
    return Re_m


def nusselt_number_lam(Re_m, Pr):
    """ mittlere Nu-Zahl - (8.16) S.231 """
    N_lam = 0.664 * pow(Re_m, 1 / 2) * pow(Pr, 1 / 3)
    return N_lam


def heat_exchange_coefficient_lam(N_lam, LAMBDA, length):
    hec_lam = (N_lam * LAMBDA) / length
    return hec_lam


def nusselt_number_turb(Re_m, Pr):
    """ mittlere Nu-Zahl - (8.21) S.232 """
    N_turb = (0.037 * pow(Re_m, 0.8) * Pr) / (1 + 2.443 * pow(Re_m, -0.1) * (pow(Pr, 2 / 3) - 1))
    return N_turb


def heat_exchange_coefficient_turb(nusselt_number_turb, LAMBDA, length):
    hec_turb = (nusselt_number_turb * LAMBDA) / length
    return hec_turb


def nusselt_number_0(nusselt_number_lam, nusselt_number_turb):
    N_0 = sqrt(pow(nusselt_number_lam, 2) + pow(nusselt_number_turb, 2))
    return N_0


def correction_factor(Temperature_surface, Temperature_reference):
     """Temperature_surface = Bezugstemperatur """
    cf = pow(Temperature_reference / Temperature_surface, 0.12)
    return cf


def nusselt_number_m(cf, N_0):
    N_m = cf * N_0
    return N_m


def heat_exchange_coefficient(N_m, LAMBDA, length):
    """Calculates the heat exchange coefficient (Wärmeübergangskoeffizient), WUE P.232"""
    hec = (N_m * LAMBDA) / length
    return hec


def convective_heat_flow(heat_exchange_coefficient, plate_area, air_temp):
    """Calculates the convective heat flow across the plate"""
    chf = heat_exchange_coefficient * plate_area * air_temp
    return chf


def dew_temperature(rel_humidity, air_temp):
    """Calculates dewpoint temperature as referred to in The Relationship between Relative
        Humidity and the Dewpoint Temperature in Moist Air A Simple Conversion and Applications
        BY MARK G. LAWRENCE 2005 """
    b = 243.04
    a = 17.625
    dew_temp = b * (log(rel_humidity / 100) + (a * air_temp / (b + air_temp))) / (
            a - log(rel_humidity / 100) - ((a * air_temp) / (b + air_temp)))
    return dew_temp


def sky_temperature(air_temp, dew_temp, standard_time):
    """Calculates the sky temperature as referred to in Solar Engineering, Duffie and Beckman, Page 148
        if clouds are included in the calculation then the sky temperature must be higher than a cloudless sky"""
    sky_temp = air_temp * pow(0.711 + 0.0056 * dew_temp + 0.000073 * pow(dew_temp, 2) + 0.013 * cos(15 * standard_time),
                              1 / 4)
    return sky_temp


def radiation_heat_flow(em_deg, area, plate_temp,
                        sky_temp):  # Emission degree of Plate is user defined or given by chosen material
    """Calculates the heat flow resulting from radiation exchange with the sky, WUE P.181"""
    rhf = em_deg * sigma * area * (pow(plate_temp, 4) - pow(sky_temp, 4))
    return rhf


def solar_heat_flow(solar_absorbtion_coeff, area,
                    irradiation_g):  # Solar absorbtion coefficient is a user input or material default
    """Calculates the heat flow resulting from solar radiation"""
    shf = solar_absorbtion_coeff * area * irradiation_g
    return shf


def plateTemp(em_deg, length, width, solar_absorbtion_coeff, air_temp, irradiation_g, plate_temp, rel_humidity):
    area = length * width
    hf_sol = solar_heat_flow(solar_absorbtion_coeff, area, irradiation_g)
    dew_temp = dew_temperature(rel_humidity, air_temp)
    sky_temp = sky_temperature(air_temp, dew_temp)
    hf_rad = radiation_heat_flow(em_deg, area, plate_temp, sky_temp)
    cf = correction_factor(plate_temp, )
    N_m = nusselt_number_m(cf,N_0)
    heat_ex_coeff = heat_exchange_coefficient(N_m, LAMBDA, length)
    hf_conv = convective_heat_flow(heat_ex_coeff, area, air_temp)


