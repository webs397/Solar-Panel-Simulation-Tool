import requests
import json


def jprint(obj):
    text = json.dumps(obj, sort_keys=True, indent=4)
    print(text)


def get_daily_data(lat, lon, month, angle, azimuth, b_global, b_clearsky, b_temp, b_localtime):
    """Contacts the PVGIS API  and requests their average daily irradiance data for the given information"""
    response = requests.get("https://re.jrc.ec.europa.eu/api/DRcalc", params=set_parameters(lat, lon, month, angle, azimuth, b_global, b_clearsky, b_temp, b_localtime))
    # print(response.url)
    jprint(response.json())


def set_parameters(lat, lon, month, angle, azimuth, b_global, b_clearsky, b_temp, b_localtime):
    parameters = {

        "lat": lat,
        "lon": lon,
        "month": month,
        "angle": angle,
        "aspect": azimuth,
        "global": b_global,
        "clearsky": b_clearsky,
        "showtemperature": b_temp,
        "localtime": b_localtime,
        "outputformat": 'json'

    }
    return parameters


get_daily_data(50.77534, 6.0838868, 1, 0, 0, 1, 0, 0, 0)

"""+ lat + "&lon=" + lon + "&month=" + month
                            + "&angle=" + angle + "&aspect=" + azimuth + "&global=" + b_global + "&clearsky=" + b_clearsky
                            + "&showtemperatures=" + b_temp + "&localtime=" + b_localtime + "&outputformat=json")"""
