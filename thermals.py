from scipy.constants import sigma
from math import cos, sqrt
from numpy import log

LAMBDA = 0.0259580
dynamic_viscosity_g = 1.8264e-05
specific_isobaric_heat_capacity_g = 1.005470208
density_g = 1.420943445
density_l = 1.16110661
cp = 1.00583148
eta = 1.855E-05


def angle_to_vertical(angle):
    angle_vert = 90 - angle
    return angle_vert


def incoming_heat_flow(a_s, plate_area, irradiation_global):
    """Calculates the incoming heat flow from global irradiation, output is in Watts"""
    ihf = a_s * plate_area * irradiation_global
    return ihf


def Prandtl(dynamic_viscosity_g, specific_isobaric_heat_capacity_g, LAMBDA):
    Pr = (dynamic_viscosity_g * specific_isobaric_heat_capacity_g) / LAMBDA
    return Pr


print(Prandtl(dynamic_viscosity_g, specific_isobaric_heat_capacity_g, LAMBDA))


def Reynold_m(wind_velocity, length, dynamic_viscosity_g, density_g):
    """ mittlere Reynolds-Zahl - (8.17) S.231 """
    Re_m = (wind_velocity * length * density_g) / dynamic_viscosity_g
    return Re_m


def nusselt_number_lam(Re_m, Pr):
    """ mittlere Nu-Zahl - (8.16) S.231 """
    Nu_lam = 0.664 * pow(Re_m, 1 / 2) * pow(Pr, 1 / 3)
    return Nu_lam


def heat_exchange_coefficient_lam(Nu_lam, LAMBDA, length):
    hec_lam = (Nu_lam * LAMBDA) / length
    return hec_lam


def nusselt_number_turb(Re_m, Pr):
    """ mittlere Nu-Zahl - (8.21) S.232 """
    N_turb = (0.037 * pow(Re_m, 0.8) * Pr) / (1 + 2.443 * pow(Re_m, -0.1) * (pow(Pr, 2 / 3) - 1))
    return N_turb


def heat_exchange_coefficient_turb(nusselt_number_turb, LAMBDA, length):
    hec_turb = (nusselt_number_turb * LAMBDA) / length
    return hec_turb


def nusselt_number_erzw(nusselt_number_lam, nusselt_number_turb):
    Nu_erzw = sqrt(pow(nusselt_number_lam, 2) + pow(nusselt_number_turb, 2))
    return Nu_erzw


def reference_temperature(plate_temp, air_temp):
    ref_temp = (plate_temp + air_temp) / 2
    return ref_temp


def correction_factor(plate_temp, temp_ref):
    cf = pow(temp_ref / plate_temp, 0.12)
    return cf


def nusselt_number_erzw_corrected(cf, Nu_erzw):
    Nu_erzw_corrected = cf * Nu_erzw
    return Nu_erzw_corrected


def Beta_gas(temp_ref):
    Beta_g = 1 / temp_ref
    return Beta_g


def Rayleigh_number(Beta_gas, plate_temp, air_temp, length, density_l, cp, eta, LAMBDA):
    Ra = 1000 * (9.80665 * Beta_gas * (plate_temp - air_temp) * pow(length, 3) * pow(density_l, 2) * cp) / (
                eta * LAMBDA)
    return Ra


def Rayleigh_number_critical(angle_vert):
    Ra_c = pow(10, 8.9 - 0.00178 * pow(angle_vert, 1.82))
    return Ra_c


def nusselt_number_free(Ra_c, angle_vert, Ra):
    Nu_free = 0.56 * pow(Ra_c * cos(angle_vert), 1 / 4) + 0.13 * (pow(Ra, 1 / 3) - pow(Ra_c, 1 / 3))
    return Nu_free


def nusselt_number_mix(Nu_erz_corrected, Nu_free, angle):
    # Might need to add if statement for angle
    Nu_mix = pow(pow(Nu_erz_corrected, 3) + pow(Nu_free, 3), 1 / 3)
    return Nu_mix


def heat_exchange_coefficient(Nu_mix, LAMBDA, length):
    """Calculates the heat exchange coefficient (Wärmeübergangskoeffizient), WUE P.232"""
    hec = (Nu_mix * LAMBDA) / length
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


def radiation_heat_flow(em_fac, area, plate_temp,
                        sky_temp):  # Emission degree of Plate is user defined or given by chosen material
    """Calculates the heat flow resulting from radiation exchange with the sky, WUE P.181"""
    rhf = em_fac * sigma * area * (pow(plate_temp, 4) - pow(sky_temp, 4))
    return rhf


def solar_heat_flow(ab_fac, area,
                    irradiation_g):  # Solar absorbtion coefficient is a user input or material default
    """Calculates the heat flow resulting from solar radiation"""
    shf = ab_fac * area * irradiation_g
    return shf


def plateTemp(em_fac, length, width, ab_fac, air_temp, irradiation_g, plate_temp, rel_humidity, angle,
              wind_vel):
    angle_vert = angle_to_vertical(angle)
    area = length * width
    hf_sol = solar_heat_flow(ab_fac, area, irradiation_g)
    dew_temp = dew_temperature(rel_humidity, air_temp)
    sky_temp = sky_temperature(air_temp, dew_temp)
    hf_rad = radiation_heat_flow(em_fac, area, plate_temp, sky_temp)
    temp_ref = reference_temperature(plate_temp, air_temp)
    cf = correction_factor(plate_temp, temp_ref)
    Pr = Prandtl(dynamic_viscosity_g, specific_isobaric_heat_capacity_g, LAMBDA)
    Re_m = Reynold_m(wind_vel, length, dynamic_viscosity_g, density_g)
    Nu_lam = nusselt_number_lam(Re_m, Pr)
    Nu_turb = nusselt_number_turb(Re_m, Pr)
    Ra_c = Rayleigh_number_critical(angle_vert)
    Ra = Rayleigh_number(Beta_gas(temp_ref), plate_temp, air_temp, length, density_l, cp, eta, LAMBDA)
    Nu_free = nusselt_number_free(Ra_c, angle_vert, Ra)
    Nu_erz = nusselt_number_erzw(Nu_lam, Nu_turb)
    Nu_erz_corrected = nusselt_number_erzw_corrected(cf, Nu_erz)
    N_mix = nusselt_number_mix(Nu_erz_corrected, Nu_free, angle)
    heat_ex_coeff = heat_exchange_coefficient(N_mix, LAMBDA, length)
    hf_conv = convective_heat_flow(heat_ex_coeff, area, air_temp)
