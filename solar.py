from math import sin, cos, pi, radians, degrees



def declination(day):
    """Declination between the earth and sun"""
    # Could get data from thesunlive.com instead?
    b = radians(360) * (day - 1) / 365
    delta = (radians(180) / pi) * (0.006918 - 0.399912 * cos(b) + 0.070257 * sin(b) - 0.006758 *
                                   cos(2 * b) + 0.000907 * sin(2 * b) - 0.002697 * cos(3 * b) +
                                   0.00148 * sin(3 * b))
    return degrees(delta)


def equation_of_time(day_number):
    """Difference between true time and mean solar time"""
    b = (day_number - 1) * (360 / 365)
    b = radians(b)
    eot = 229.18 * (0.000075 + 0.001868 * cos(b) - 0.032077 * sin(b) - 0.014615 * cos(2 * b) - 0.04089 * sin(2 * b))

    return eot


def standard_meridian(timezone):
    std_mer = timezone * 15

    return std_mer


def mean_solar_time_baer(standard_time, longitude, day_number):
    """Calculate location based DST using Baer Formula (only valid for Europe)"""
    time_offset = standard_time - 12 + ((longitude / 15) - 1) + equation_of_time(day_number) / 60
    st = standard_time + time_offset
    return st


def mean_solar_time_duffie_beckman(standard_time, timezone, longitude, day_number):  # timezone is UTC
    """Calculates solar time according to John A. Duffie, William A. Beckman - Solar Engineering of Thermal Processes
    2013, output is in hours"""
    st = standard_time + (4 * (standard_meridian(timezone) - longitude) + equation_of_time(day_number)) / 60
    return st


def true_solar_time_noaa(standard_time, timezone, longitude, day_number):  # timezone is UTC
    """Calculates solar time according to NOAA: https://gml.noaa.gov/grad/solcalc/solareqns.PDF"""
    time_offset = (equation_of_time(day_number) / 60 + 4 * longitude) / 60 - timezone
    tst = standard_time + time_offset
    return tst


def solar_hour_angle(true_solar_time):
    """Calculates the solar hour angle in degrees [From Feuerriegel WUE Skript]"""
    sha = 15 * true_solar_time
    return sha


def angle_sun_plate():  # NEED TO DO THIS PAGE 175 WUE SKRIPT
    return


def air_mass_factor(m, m0):
    amf = m / m0
    return amf


def relative_optic_air_mass(theta_z):
    roam = 1 / cos(theta_z)
    return roam







"""
print("NOAA: ", mean_solar_time_noaa(12.18, 6, 89.4, 34))
print("Duffie,Beckman: ", mean_solar_time_duffie_beckman(10.5, -6, 89.4, 34))
print("Baer: ", mean_solar_time_baer(11, 13.35, 244))
print(difference_checker(12, 2, 13.35, 244))

# NOAA has same result as Baer example with time = 11
# NOAA has virtually same result as Duffie and Beckman (0.1 off)
"""
