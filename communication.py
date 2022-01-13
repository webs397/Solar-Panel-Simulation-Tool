import requests
import json


def jprint(obj):
    text = json.dumps(obj, sort_keys=True, indent=4)
    #print(text)


def get_daily_data(lat, lon, month, angle, azimuth, b_global, b_clearsky, b_temp, b_localtime):
    """Contacts the PVGIS API  and requests their average daily irradiance data for the given information"""
    response = requests.get("https://re.jrc.ec.europa.eu/api/DRcalc", params=set_parameters(lat, lon, month, angle, azimuth, b_global, b_clearsky, b_temp, b_localtime))
    #print(response.url)
    # jprint(response.json())
    data = json.loads(response.text)
    hours = []
    for i in range(0, len(data['outputs']['daily_profile'])):
        hours.append(data['outputs']['daily_profile'][i]["time"])

    temps = []
    for i in range(0, len(data['outputs']['daily_profile'])):
        temps.append(data['outputs']['daily_profile'][i]["T2m"])

    i_global = []
    for i in range(0, len(data['outputs']['daily_profile'])):
        #print("Global: ", data['outputs']['daily_profile'][i]["G(i)"])
        i_global.append(data['outputs']['daily_profile'][i]["G(i)"])

    i_clearsky = []
    if b_clearsky:
        for i in range(0, len(data['outputs']['daily_profile'])):
            #print("Global: ", data['outputs']['daily_profile'][i]["Gcs(i)"])
            i_clearsky.append(data['outputs']['daily_profile'][i]["Gcs(i)"])


    proc_data = (hours, temps, i_global, i_clearsky)
    return proc_data


def set_parameters(lat, lon, month, angle, azimuth, b_global, b_clearsky, b_temp, b_localtime):
    parameters = {
        "lat": lat,
        "lon": lon,
        "month": month,
        "angle": angle,
        "aspect": azimuth,
        "global": b_global,
        "clearsky": b_clearsky,
        "showtemperatures": b_temp,
        "localtime": b_localtime,
        "outputformat": 'json'

    }
    return parameters


get_daily_data(50.77534, 6.0838868, 1, 0, 0, 1, 1, 1, 1)

