import sys
from PyQt5 import QtCore, QtGui, QtWidgets


class Demo(QtWidgets.QGraphicsView):
    def __init__(self):
        super(Demo, self).__init__()
        self._scene = QtWidgets.QGraphicsScene()
        self._scene.setSceneRect(0, 0, 300, 300)
        self.setScene(self._scene)
        self.rect0 = self._scene.addRect(
            10, 30, 300, 300, brush=QtGui.QBrush(QtGui.QColor("grey"))
        )
        self.rect1 = QtWidgets.QGraphicsRectItem(
            100, 30, 100, 30) #, brush=QtGui.QBrush(QtGui.QColor("red"))
        # ~ )
        self.rect1.setFlag(QtWidgets.QGraphicsItem.ItemIgnoresTransformations)

        self.rect2 = QtWidgets.QGraphicsRectItem(
            200, 30, 100, 30) #, brush=QtGui.QBrush(QtGui.QColor("green"))
        # ~ )

        self.rect1.setParentItem(self.rect0)
        self.rect2.setParentItem(self.rect0)

        self.rotate(50)
        self._use_deviceTransform = False

    def mousePressEvent(self, event):
        sp = self.mapToScene(event.pos())
        item = self._scene.itemAt(
            sp,
            self.viewportTransform()
            if self._use_deviceTransform
            else QtGui.QTransform(),
        )
        print(item)

    def set_use_deviceTransform(self, t):
        self._use_deviceTransform = t


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    radiobutton = QtWidgets.QRadioButton("use deviceTransform")
    demo = Demo()
    radiobutton.toggled.connect(demo.set_use_deviceTransform)
    w = QtWidgets.QWidget()
    lay = QtWidgets.QVBoxLayout(w)
    lay.addWidget(radiobutton)
    lay.addWidget(demo)
    w.show()
    w.resize(640, 480)
    sys.exit(app.exec_())
