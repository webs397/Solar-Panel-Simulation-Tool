import sys
import io
from ModelPanel import SolarPanel
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
import folium
from PyQt5.QtWebEngineWidgets import QWebEngineView
import pyqtgraph as pg
from communication import get_daily_data


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




def openWindow(self, data):
    self.w = outputWindow()
    output_layout = QHBoxLayout()
    temperature_graph = pg.PlotWidget()
    hours = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
    temperature_graph.setBackground('w')
    airtemp_pen = pg.mkPen(color=(255, 0, 0), width=3)
    temp_pen = pg.mkPen(color=(0, 0, 0), width=3)
    temperature_graph.setXRange(0, 23, padding=0)
    temperature_graph.setLimits(xMin=0, xMax=23)
    temperature_graph.plot(hours, data[1], pen=airtemp_pen)

    # PLATE TEMPERATURE IN HERE
    #temperature_graph.plot(hours, , pen=temp_pen)
    output_layout.addWidget(temperature_graph)

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
        meteo_data_layout.addWidget(QLabel("Wind Direction [°]"))
        line_edit_winddirection = QLineEdit()
        meteo_data_layout.addWidget(line_edit_winddirection)
        check_g_irrad = QCheckBox("Global Irradiance")
        meteo_data_layout.addWidget(check_g_irrad)
        check_clearsky = QCheckBox("Clear Sky")
        meteo_data_layout.addWidget(check_clearsky)
        check_temp = QCheckBox("Show Temperature")
        meteo_data_layout.addWidget(check_temp)


        # Setup for Plate Data

        sp = SolarPanel()
        plate_data_layout.addWidget(QLabel("Angle [°]"))
        line_edit_angle = NumberLineEdit()
        line_edit_angle.sendNumber.connect(sp.setAlpha)
        plate_data_layout.addWidget(line_edit_angle)
        plate_data_layout.addWidget(QLabel("Azimuth [°]"))
        plate_data_layout.addWidget(QLabel("0=South, 90=West, -90=East"))
        line_edit_azimuth = NumberLineEdit()
        line_edit_azimuth.sendNumber.connect(sp.setGamma)
        plate_data_layout.addWidget(line_edit_azimuth)
        #verticalSpacer = QSpacerItem(10, 10)
        #plate_data_layout.addSpacerItem(verticalSpacer)


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
            inputs = [cb_month.currentText(), line_edit_angle.text(), line_edit_azimuth.text(),
                      line_edit_windspeed.text(), line_edit_winddirection.text(), ]
            inputs[0] = monthToInt(inputs[0])
            coords = tuple(float(x) for x in coord_line_edit.text().split(','))
            if check_g_irrad.isChecked():
                b_global = 1
            else:
                b_global = 0
            if check_clearsky.isChecked():
                b_clearsky = 1
            else:
                b_clearsky = 0
            if check_temp.isChecked():
                b_temp = 1
            else:
                b_temp = 0
            pvgis_data = get_daily_data(coords[0], coords[1], inputs[0], inputs[1], inputs[2], b_global,
                                        b_clearsky, b_temp, 0)

            openWindow(self, pvgis_data)

            print(coords)
            print(b_global)
            print(b_clearsky)
            print(b_temp)
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
