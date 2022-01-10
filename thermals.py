from scipy.constants import sigma
from math import cos
from scipy.constants import lambda



def incoming_heat_flow(a_s, plate_area, irradiation_global):
    """Calculates the incoming heat flow from global irradiation, output is in Watts"""
    ihf = a_s * plate_area * irradiation_global
    return ihf

def Prandtl(dynamic_viscosity_g, specific_isobaric_heat_capacity_g , lambda)
    Pr = (dynamic_viscosity_g  * specific_isobaric_heat_capacity_g) / lambda
    return Rr

def Reynold_m(wind_velocity, length, dynamic_viscosity_g, density_g)
    """ mittlere Reynolds-Zahl - (8.17) S.231 """
    Re_m = (wind_velocity * length * density_g) / dynamic_viscosity_g
    return Re_m

def nusselt_number_lam(Reynold_m, Prandtl)
    """ mittlere Nu-Zahl - (8.16) S.231 """
    N_lam = 0,664 * Reynold_m^1/2 * Prandtl^1/3
    return N_lam

def heat_exchange_coefficient_lam(nusselt_number_lam, lambda, length)
    hec_lam = (nusselt_number_lam * lambda)/ length 
    return hec_lam

def nusselt_number_turb(Reynold_m, Prandtl)
""" mittlere Nu-Zahl - (8.21) S.232 """
    N_turb = (0,037 * Reynold_m^0,8 * Prandtl)/(1 + 2,443 * Reynold_m^-0,1 * (Prandtl^(2/3)-1))
    return N_turb

def heat_exchange_coefficient_turb(nusselt_number_turb, lambda, length)
    hec_turb = (nusselt_number_turb * lambda)/ length 
    return hec_turb

def nusselt_number_0(nusselt_number_lam,nusselt_number_turb)
    N_0 = sqrt((nusselt_number_lam)^2+(nusselt_number_turb)^2)
    return N_0

def correction_factor(Temperature_surface, Temperature_reference)
    cf = (Temperature_reference/Temperature_surface)^0,12
    return cf

def nusselt_number_m(correctionfactor, nusselt_number_0)
    N_m = correctionfactor * nusselt_number_o
    return N_m
    

def heat_exchange_coefficient(nusselt_number_m, lambda, length):
    """Calculates the heat exchange coefficient (Wärmeübergangskoeffizient), WUE P.232"""
    hec = (nusselt_number_m * lambda)/ length 
    return hec

                     
def convective_heat_flow(heat_exchange_coefficient, plate_area, temperature_air):
    """Calculates the convective heat flow across the plate"""
    chf = heat_exchange_coefficient* plate_area * temperature_air
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
