import io
import sys

import folium
import pyqtgraph as pg
import scipy
from PyQt5.QtCore import *
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import *

from ModelPanel import SolarPanel
from communication import get_daily_data
from thermals import f, solar_heat_flow, radiation_heat_flow, convective_heat_flow, heat_exchange_coefficient, \
    temp_kelvin, dew_temperature, sky_temperature, reference_temperature, correction_factor, Reynold_m, \
    dynamic_viscosity_l, density_l, nusselt_number_lam, prandtl, nusselt_number_turb, Rayleigh_number, gravity, \
    kinematik_visc_l, temp_conductivity_l, nusselt_number_free_24, nusselt_number_erzw, nusselt_number_erzw_corrected, \
    nusselt_number_mix, Beta_gas


class TableView(QTableWidget):
    def __init__(self, data, *args):
        QTableWidget.__init__(self, *args)
        self.data = data
        self.setData()
        self.resizeColumnsToContents()
        self.resizeRowsToContents()
        self.setSortingEnabled(False)

    def setData(self):
        verHeaders = []
        for n, key in enumerate((self.data.keys())):
            verHeaders.append(key)
            for m, item in enumerate(self.data[key]):
                newitem = QTableWidgetItem(item)
                self.setItem(n, m, newitem)
        self.setVerticalHeaderLabels(verHeaders)


class NumberLineEdit(QLineEdit):
    sendNumber = pyqtSignal(float)

    def __init__(self, default=0):
        super().__init__()

        self.default = default
        self.setText(str(default))
        self.lastInput = str(default)
        self.editingFinished.connect(self.onEditingFinished)

    def onEditingFinished(self):
        try:
            float(self.text())
            self.lastInput = self.text()
        except ValueError:
            self.setText(self.lastInput)

        self.sendNumber.emit(float(self.text()))

    def setText(self, a0: str) -> None:
        return super().setText(str(a0))


class CoordLineEdit(QLineEdit):
    def __init__(self, web_view, default=(50.77555, 6.083611)):
        super().__init__()

        self.web_view = web_view
        self.default = default
        self.setText(str(default[0]) + ", " + str(default[1]))
        self.lastInput = str(default[0]) + ", " + str(default[1])
        self.editingFinished.connect(self.onEditingFinished)

    def onEditingFinished(self):
        try:
            coords = tuple(float(x) for x in self.text().split(','))
            self.lastInput = self.text()
        except ValueError:
            self.setText(self.lastInput)
            coords = tuple(float(x) for x in self.text().split(','))

        m = folium.Map(zoom_start=13, location=coords)
        map_data = io.BytesIO()
        m.save(map_data, close_file=False)
        self.web_view.setHtml(map_data.getvalue().decode())


def monthToInt(month):
    if month == "January":
        return 1
    if month == "February":
        return 2
    if month == "March":
        return 3
    if month == "April":
        return 4
    if month == "May":
        return 5
    if month == "June":
        return 6
    if month == "July":
        return 7
    if month == "August":
        return 8
    if month == "September":
        return 9
    if month == "October":
        return 10
    if month == "November":
        return 11
    if month == "December":
        return 12


def stringify(list):
    for i in range(len(list)):
        list[i] = str(round(list[i], 2))
    return list


def openWindow(self, data, inputs):
    self.w = outputWindow()
    output_layout = QVBoxLayout()
    graph_layout = QHBoxLayout()
    hours = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]

    # Irradiance Graph
    irrad_graph = pg.PlotWidget()
    irrad_graph.setBackground('w')
    irrad_graph.addLegend()
    global_pen = pg.mkPen(color=(0, 0, 0), width=3)
    clearsky_pen = pg.mkPen(color=(0, 191, 255), width=3)
    irrad_graph.setXRange(0, 23, padding=0)
    irrad_graph.setLimits(xMin=0, xMax=23)
    if inputs[9]:
        irrad_graph.plot(hours, data[3], pen=clearsky_pen, name="Clear Sky Irradiance")

    else:
        irrad_graph.plot(hours, data[2], pen=global_pen, name="Global Irradiance")
    irrad_graph.setXRange(0, 23, padding=0)
    irrad_graph.setLimits(xMin=0, xMax=23)
    irrad_graph.setYRange(0, 1200)
    irrad_graph.setLimits(yMin=0, yMax=1200)
    irrad_graph.setTitle("Irradiance")
    irrad_graph.setLabel('left', "<span <p>W/m <sup>2</sup></p> </span>")
    irrad_graph.setLabel('bottom', 'Hours')


    # Temperature Graph
    temperature_graph = pg.PlotWidget()
    temperature_graph.setBackground('w')
    temperature_graph.addLegend()
    airtemp_pen = pg.mkPen(color=(255, 0, 0), width=3)
    temp_pen = pg.mkPen(color=(0, 0, 0), width=3)
    temperature_graph.setXRange(0, 23, padding=0)
    temperature_graph.setLimits(xMin=0, xMax=23)
    if inputs[10]:
        temperature_graph.plot(hours, data[1], name="Air Temperature", pen=airtemp_pen)

    plate_temps = []
    if inputs[9]:  # This checks for clear sky bool
        for i in range(0, len(data[1])):
            #print("Inputs: ", inputs[0], inputs[1], inputs[2], inputs[3], data[1][i], data[3][i], inputs[4], inputs[5],
                  #inputs[6], hours[i])
            a = scipy.optimize.newton(f, 0, args=[inputs[0], inputs[1], inputs[2], inputs[3], data[1][i], data[3][i],
                                                  inputs[4], inputs[5], inputs[6], hours[i]], maxiter=10000)
            #print("Hour: ", hours[i], " Plate Temp: ", a.real)
            plate_temps.append(a.real)
        temperature_graph.plot(hours, plate_temps, name="Solar Panel Temperature", pen=temp_pen)
    else:
        for i in range(0, len(data[1])):
            #print("Inputs: ", inputs[0], inputs[1], inputs[2], inputs[3], data[1][i], data[2][i], inputs[4], inputs[5],
                  #inputs[6], hours[i])
            a = scipy.optimize.newton(f, 0, args=[inputs[0], inputs[1], inputs[2], inputs[3], data[1][i], data[2][i],
                                                  inputs[4], inputs[5], inputs[6], hours[i]], maxiter=10000)
            #print("Hour: ", hours[i], " Plate Temp: ", a.real)
            plate_temps.append(a.real)
        temperature_graph.plot(hours, plate_temps, name="Solar Panel Temperature", pen=temp_pen)
    temperature_graph.setYRange(-10, 70)
    temperature_graph.setLimits(yMin=-10, yMax=70)
    temperature_graph.setTitle("Temperatures")
    temperature_graph.setLabel('left', '°C')
    temperature_graph.setLabel('bottom', 'Hours')


    # These calculations are based on a LG300S1V-A5 Solar Panel
    # Efficiency Graph
    eff_graph = pg.PlotWidget()
    eff_graph.setBackground('w')
    eff_graph.addLegend()
    eff_theo_pen = pg.mkPen(color=(0, 0, 0), width=3, style=Qt.PenStyle.DashDotLine)
    eff_real_pen = pg.mkPen(color=(0, 0, 0), width=3)
    eff_theo = []
    eff_real = []
    for i in range(0, len(data[1])):
        eff_calc_theo = -0.41 * plate_temps[i] + 110.25
        eff_theo.append(eff_calc_theo)
        eff_calc_real = -0.41 * plate_temps[i] + 110.25
        if eff_calc_real >= 100:
            eff_calc_real = 100
        eff_real.append(eff_calc_real)
    eff_graph.plot(hours, eff_theo, name="Theoretical Efficiency", pen=eff_theo_pen)
    eff_graph.plot(hours, eff_real, name="Capped Efficiency", pen=eff_real_pen)
    eff_graph.setXRange(0, 23)
    eff_graph.setLimits(xMin=0, xMax=23)
    eff_graph.setYRange(75, 120)
    eff_graph.setLimits(yMin=75, yMax=120)
    eff_graph.setTitle("Efficiency")
    eff_graph.setLabel('left', '%')
    eff_graph.setLabel('bottom', 'Hours')


    # Power Graph
    power_graph = pg.PlotWidget()
    power_graph.setBackground('w')
    power_graph.addLegend()
    power_theo_pen = pg.mkPen(color=(255, 215, 0), style=Qt.PenStyle.DashDotLine)
    power_real_pen = pg.mkPen(color=(255, 215, 0), width=3)
    power_theo = []
    power_real = []
    if inputs[9]:
        for i in range(0, len(data[1])):
            # eff_dif_real = eff_real[i]-100
            # eff_dif_theo = eff_theo[i]-100
            power_calc_theo = (0.175 * eff_theo[i] / 100) * data[3][i] * inputs[1] * inputs[2]
            power_calc_real = (0.175 * eff_real[i] / 100) * data[3][i] * inputs[1] * inputs[2]
            if power_calc_theo >= 300:
                power_calc_theo = 300
            power_theo.append(power_calc_theo)
            if power_calc_real >= 300:
                power_calc_real = 300
            power_real.append(power_calc_real)
        power_graph.plot(hours, power_theo, name="Theoretical Power", width=3, pen=power_theo_pen)
        power_graph.plot(hours, power_real, name="Real Power", pen=power_real_pen)
    else:
        for i in range(0, len(data[1])):
            # eff_dif_real = eff_real[i]-100
            # eff_dif_theo = eff_theo[i]-100
            power_calc_theo = (0.175 * eff_theo[i] / 100) * data[2][i] * inputs[1] * inputs[2]
            power_calc_real = (0.175 * eff_real[i] / 100) * data[2][i] * inputs[1] * inputs[2]
            if power_calc_theo >= 300:
                power_calc_theo = 300
            power_theo.append(power_calc_theo)
            if power_calc_real >= 300:
                power_calc_real = 300
            power_real.append(power_calc_real)
        power_graph.plot(hours, power_theo, name="Theoretical Power", pen=power_theo_pen)
        power_graph.plot(hours, power_real, name="Real Power", pen=power_real_pen)
    power_graph.setXRange(0, 23)
    power_graph.setLimits(xMin=0, xMax=23)
    power_graph.setYRange(0, 300)
    power_graph.setLimits(yMin=0, yMax=300)
    power_graph.setTitle("Power")
    power_graph.setLabel('left', 'Watts')
    power_graph.setLabel('bottom', 'Hours')


    # Table

    def calc_heat_flows(x, em_fac, length, width, ab_fac, air_temp, irradiation_g, rel_humidity, angle,
                        wind_vel, time):
        x = temp_kelvin(x)
        air_temp = temp_kelvin(air_temp)
        plate_temp = x
        # angle_vert = angle_to_vertical(angle)
        area = length * width
        hf_sol = solar_heat_flow(ab_fac, area, irradiation_g)
        dew_temp = dew_temperature(rel_humidity, air_temp)
        sky_temp = sky_temperature(air_temp, dew_temp, time)
        hf_rad = radiation_heat_flow(em_fac, area, plate_temp, sky_temp)
        temp_ref = reference_temperature(plate_temp, air_temp)
        cf = correction_factor(plate_temp, temp_ref)
        Re_m = Reynold_m(wind_vel, length, dynamic_viscosity_l, density_l)
        Nu_lam = nusselt_number_lam(Re_m, prandtl)
        Nu_turb = nusselt_number_turb(Re_m, prandtl)
        # Ra_c = Rayleigh_number_critical(angle_vert)
        Ra = Rayleigh_number(gravity, Beta_gas(temp_ref), plate_temp, air_temp, length, kinematik_visc_l,
                             temp_conductivity_l)
        Nu_free = nusselt_number_free_24(Ra, angle)
        Nu_erz = nusselt_number_erzw(Nu_lam, Nu_turb)
        Nu_erz_corrected = nusselt_number_erzw_corrected(cf, Nu_erz)
        N_mix = nusselt_number_mix(Nu_erz_corrected, Nu_free)
        heat_ex_coeff = heat_exchange_coefficient(N_mix, length)
        hf_conv = convective_heat_flow(heat_ex_coeff, area, air_temp, plate_temp)
        return [hf_sol.real, hf_conv.real, hf_rad.real]

    heat_flow_sol = []
    heat_flow_conv = []
    heat_flow_rad = []

    if inputs[9]:
        for i in range(len(hours)):
            v = calc_heat_flows(plate_temps[i], inputs[0], inputs[1], inputs[2], inputs[3], data[1][i], data[3][i],
                                inputs[4], inputs[5], inputs[6], hours[i])
            heat_flow_sol.append(v[0])
            heat_flow_conv.append(v[1])
            heat_flow_rad.append(v[2])

        data = {'Irradiance [W/m^2]': stringify(data[3]),
                'Air Temperature [°C]': stringify(data[1]),
                'Solar Heat Flow [W]': stringify(heat_flow_sol),
                'Convective Heat Flow [W]': stringify(heat_flow_conv),
                'Radiation Heat Flow [W]': stringify(heat_flow_rad),
                'Plate Temperature [°C]': stringify(plate_temps),
                'Efficiency Theoretical [%]': stringify(eff_theo),
                'Efficiency Real [%]': stringify(eff_real),
                'Power Theoretical [W]': stringify(power_theo),
                'Power Real [W]': stringify(power_real)}
    else:
        for i in range(len(hours)):
            v = calc_heat_flows(plate_temps[i], inputs[0], inputs[1], inputs[2], inputs[3], data[1][i], data[2][i],
                                inputs[4], inputs[5], inputs[6], hours[i])
            heat_flow_sol.append(v[0])
            heat_flow_conv.append(v[1])
            heat_flow_rad.append(v[2])

        data = {'Irradiance [W/m^2]': stringify(data[2]),
                'Air Temperature [°C]': stringify(data[1]),
                'Solar Heat Flow [W]': stringify(heat_flow_sol),
                'Convective Heat Flow [W]': stringify(heat_flow_conv),
                'Radiation Heat Flow [W]': stringify(heat_flow_rad),
                'Plate Temperature [°C]': stringify(plate_temps),
                'Efficiency Theoretical [%]': stringify(eff_theo),
                'Efficiency Real [%]': stringify(eff_real),
                'Power Theoretical [W]': stringify(power_theo),
                'Power Real [W]': stringify(power_real)}

    tableWidget = TableView(data, 10, 24)
    tableWidget.resizeRowsToContents()
    graph_layout.addWidget(irrad_graph)
    graph_layout.addWidget(temperature_graph)
    graph_layout.addWidget(eff_graph)
    graph_layout.addWidget(power_graph)
    output_layout.addLayout(graph_layout)
    output_layout.addWidget(tableWidget)

    self.w.setLayout(output_layout)
    self.w.show()


class outputWindow(QWidget):
    def __init__(self, parent=None):
        super(outputWindow, self).__init__(parent)
        self.setWindowTitle("Solar Panel Simulation Tool - Results")
        # Create the different layouts
        # output_layout = QHBoxLayout()
        # left_layout = QVBoxLayout()
        # geo_data_layout = QHBoxLayout()
        # lat_lon_layout = QVBoxLayout()

        # efficiency_graph = pg.PlotWidget()
        # output_layout.addWidget(efficiency_graph)

    # power_graph = pg.PlotWidget()
    # output_layout.addWidget(power_graph)


class Window(QWidget):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        self.setWindowTitle("Solar Panel Simulation Tool")
        # Create the different layouts
        outer_layout = QHBoxLayout()
        left_layout = QVBoxLayout()
        geo_data_layout = QHBoxLayout()
        lat_lon_layout = QVBoxLayout()

        right_layout = QVBoxLayout()
        right_h_data_layout = QHBoxLayout()
        plate_data_layout = QVBoxLayout()
        meteo_data_layout = QVBoxLayout()

        # Setup for Meteorological Data

        meteo_data_layout.addWidget(QLabel("Month"))
        cb_month = QComboBox()
        cb_month.addItems(["January", "February", "March", "April", "May", "June", "July", "August", "September",
                           "October", "November", "December"])
        # self.cb_month.currentIndexChanged.connect(self.selectionchange)
        meteo_data_layout.addWidget(cb_month)
        meteo_data_layout.addWidget(QLabel("Wind Speed [m/s]"))
        # Need to add Numbers only here
        line_edit_windspeed = QLineEdit()
        meteo_data_layout.addWidget(line_edit_windspeed)
        # meteo_data_layout.addWidget(QLabel("Wind Direction [°]"))
        # line_edit_winddirection = QLineEdit()
        # meteo_data_layout.addWidget(line_edit_winddirection)
        meteo_data_layout.addWidget(QLabel("Relative Humidity [%]"))
        line_edit_rel_humidity = QLineEdit()
        meteo_data_layout.addWidget(line_edit_rel_humidity)
        # check_g_irrad = QCheckBox("Global Irradiance")
        # meteo_data_layout.addWidget(check_g_irrad)
        check_clearsky = QCheckBox("Clear Sky")
        meteo_data_layout.addWidget(check_clearsky)
        check_temp = QCheckBox("Show Air Temperature")
        meteo_data_layout.addWidget(check_temp)

        # Setup for Plate Data

        sp = SolarPanel()
        # ADD  NUMBERS VALIDATOR
        plate_data_layout.addWidget(QLabel("Panel Length [m]"))
        line_edit_length = QLineEdit()
        plate_data_layout.addWidget(line_edit_length)
        plate_data_layout.addWidget(QLabel("Panel Width [m]"))
        line_edit_width = QLineEdit()
        plate_data_layout.addWidget(line_edit_width)
        plate_data_layout.addWidget(QLabel("Angle [°]"))
        line_edit_angle = NumberLineEdit()
        line_edit_angle.sendNumber.connect(sp.setAlpha)
        plate_data_layout.addWidget(line_edit_angle)
        plate_data_layout.addWidget(QLabel("Azimuth [°]"))
        plate_data_layout.addWidget(QLabel("0=South, 90=West, -90=East"))
        line_edit_azimuth = NumberLineEdit()
        line_edit_azimuth.sendNumber.connect(sp.setGamma)
        plate_data_layout.addWidget(line_edit_azimuth)
        plate_data_layout.addWidget(QLabel("Absorption Factor"))
        line_edit_absorb_factor = QLineEdit()
        plate_data_layout.addWidget(line_edit_absorb_factor)
        plate_data_layout.addWidget(QLabel("Emission Factor"))
        line_edit_em_factor = QLineEdit()
        plate_data_layout.addWidget(line_edit_em_factor)

        # verticalSpacer = QSpacerItem(10, 10)
        # plate_data_layout.addSpacerItem(verticalSpacer)

        # Setup for Map Side Layout
        web_view = QWebEngineView()
        coordinates = (50.77555, 6.083611)
        m = folium.Map(zoom_start=13, location=coordinates)
        map_data = io.BytesIO()
        m.save(map_data, close_file=False)
        web_view.setHtml(map_data.getvalue().decode())
        left_layout.addWidget(web_view)
        left_layout.addLayout(geo_data_layout)

        # Setup for Geographic Data
        lat_lon_layout.addWidget(QLabel("Latitude, Longitude"))
        coord_line_edit = CoordLineEdit(web_view)
        lat_lon_layout.addWidget(coord_line_edit)
        left_layout.addLayout(lat_lon_layout)

        def calculateClicked():
            inputs = [float(line_edit_em_factor.text()), float(line_edit_length.text()), float(line_edit_width.text()),
                      float(line_edit_absorb_factor.text()), float(line_edit_rel_humidity.text()),
                      float(line_edit_angle.text()),
                      float(line_edit_windspeed.text()), float(line_edit_azimuth.text()), cb_month.currentText(),
                      check_clearsky.isChecked(), check_temp.isChecked()]

            inputs[8] = monthToInt(inputs[8])
            coords = tuple(float(x) for x in coord_line_edit.text().split(','))
            # if check_g_irrad.isChecked():
            b_global = 1
            # else:
            #    b_global = 0
            if check_clearsky.isChecked():
                b_clearsky = 1
            else:
                b_clearsky = 0

            pvgis_data = get_daily_data(coords[0], coords[1], inputs[8], inputs[4], inputs[7], b_global,
                                        b_clearsky, 1, 1)

            openWindow(self, pvgis_data, inputs)
            return pvgis_data

        # Nest the inner layouts into the outer layout
        right_h_data_layout.addLayout(meteo_data_layout)
        right_h_data_layout.addLayout(plate_data_layout)
        right_layout.addLayout(right_h_data_layout)
        right_layout.addWidget(sp)
        calculate_button = QPushButton("Calculate")
        calculate_button.clicked.connect(calculateClicked)
        right_layout.addWidget(calculate_button)
        outer_layout.addLayout(left_layout)
        outer_layout.addLayout(right_layout)
        self.setLayout(outer_layout)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())
