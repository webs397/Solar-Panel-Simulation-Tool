import sys
from PyQt5.QtGui     import *
from PyQt5.QtCore    import *
from PyQt5.QtWidgets import *
from math import pi, cos, sin
import numpy as np


class SolarPanel(QWidget):
    def __init__(self):
        super().__init__()

        #layout = QVBoxLayout()
        #self.setLayout(layout)

        self.alpha = 0
        self.gamma = 0

        self.panel_pos = (0, 0)
        self.panel_size = (100, 100)
        self.setPoints()

    def setPoints(self):
        self.rect_points = np.array([[self.panel_pos[0], self.panel_pos[1], 0], 
                                    [self.panel_pos[0] + self.panel_size[0], self.panel_pos[1], 0], 
                                    [self.panel_pos[0] + self.panel_size[0], self.panel_pos[1] + self.panel_size[1], 0], 
                                    [self.panel_pos[0], self.panel_pos[1] + self.panel_size[1], 0]])
        self.gray_lines = np.array([[[self.panel_pos[0], self.panel_pos[1] + self.panel_size[1] * 0.25, 0], [self.panel_pos[0] + self.panel_size[0], self.panel_pos[1] + self.panel_size[1] * 0.25, 0]],
                                    [[self.panel_pos[0], self.panel_pos[1] + self.panel_size[1] * 0.50, 0], [self.panel_pos[0] + self.panel_size[0], self.panel_pos[1] + self.panel_size[1] * 0.50, 0]],
                                    [[self.panel_pos[0], self.panel_pos[1] + self.panel_size[1] * 0.75, 0], [self.panel_pos[0] + self.panel_size[0], self.panel_pos[1] + self.panel_size[1] * 0.75, 0]],
                                    [[self.panel_pos[0] + self.panel_size[0] * 0.25, self.panel_pos[1], 0], [self.panel_pos[0] + self.panel_size[0] * 0.25, self.panel_pos[1] + self.panel_size[1], 0]],
                                    [[self.panel_pos[0] + self.panel_size[0] * 0.50, self.panel_pos[1], 0], [self.panel_pos[0] + self.panel_size[0] * 0.50, self.panel_pos[1] + self.panel_size[1], 0]],
                                    [[self.panel_pos[0] + self.panel_size[0] * 0.75, self.panel_pos[1], 0], [self.panel_pos[0] + self.panel_size[0] * 0.75, self.panel_pos[1] + self.panel_size[1], 0]]])

    def resizeEvent(self, a0: QResizeEvent) -> None:
        self.panel_pos = (0.25 * a0.size().width(), 0.25 * a0.size().height())
        self.panel_size = (0.2 * a0.size().width(), 0.8 * a0.size().height())
        self.setPoints()
        self.update()
        return super().resizeEvent(a0)                                

    def paintEvent(self, event):
        painter = QPainter(self)
        
        painter.setPen(QPen(Qt.black, 2, Qt.SolidLine))
        painter.setBrush(QBrush(Qt.black))

        points = np.array([self.rotPoint(p) for p in self.rect_points])
        qpoints = [QPoint(*self.orthographicProjection(p)) for p in points]
        painter.drawPolygon(QPolygon(qpoints), fillRule=Qt.FillRule.OddEvenFill)

        if self.isFacingForward(points[0], points[1], points[3]):
            painter.setPen(QPen(Qt.gray, 2, Qt.SolidLine))
            for line in self.gray_lines:
                p1 = self.orthographicProjection(self.rotPoint(line[0]))
                p2 = self.orthographicProjection(self.rotPoint(line[1]))

                painter.drawLine(*p1, *p2)

            painter.setPen(QPen(Qt.black, 2, Qt.SolidLine))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawPolygon(QPolygon(qpoints), fillRule=Qt.FillRule.OddEvenFill)

    def rotPoint(self, point):
        if not type(point) == np.ndarray:
            point = np.array(point)

        alpha_radians = (self.alpha / 360) * 2 * pi
        gamma_radians = (self.gamma / 360) * 2 * pi

        trans_to_center = np.array([-(self.panel_pos[0] + self.panel_size[0] * 0.50), -(self.panel_pos[1] + self.panel_size[1] * 0.50), 0])
        trans_back = np.array([self.panel_pos[0] + self.panel_size[0] * 0.50, self.panel_pos[1] + self.panel_size[1] * 0.50, 0])
        rot_x = np.array([[1, 0, 0], 
                        [0, cos(alpha_radians), -sin(alpha_radians)], 
                        [0, sin(alpha_radians), cos(alpha_radians)]])
        rot_y = np.array([[cos(gamma_radians), 0, sin(gamma_radians)], 
                        [0, 1, 0], 
                        [-sin(gamma_radians), 0, cos(gamma_radians)]])

        return trans_back + rot_y.dot(rot_x.dot(trans_to_center + point))

    def isFacingForward(self, p0, p1, p2):
        return np.cross(p1 - p0, p2 - p0)[2] > 0

    def orthographicProjection(self, point):
        return point[0:2]

    def setAlpha(self, val):
        self.alpha = float(val)
        self.update()

    def setGamma(self, val):
        self.gamma = float(val)
        self.update()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = SolarPanel()
    ex.resize(399, 399)
    ex.show()
    sys.exit(app.exec_())
