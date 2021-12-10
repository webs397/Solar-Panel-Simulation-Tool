import scipy


def incoming_heat_flow(a_s, plate_area, irradiation_global):
    """Calculates the incoming heat flow from global irradiation, output is in Watts"""
    ihf = a_s * plate_area * irradiation_global
    return ihf









def convective_heat_flow(heat_exchange_coefficient, plate_area, temperature_air):
    """Calculates the convective heat flow across the plate"""
    chf = heat_exchange_coefficient * plate_area * temperature_air
    return chf




def radiation_exchange():