"""
Figure and control panel for display of DDdata

TRN, PAIR, interpolation, etc.

Figure keeps the data, possibly many t2 frames (or DAS)

Controler has the UI elements, can be shared with multiple figures.

"""
import inspect
import numpy as np
from PyQt5 import QtWidgets
import logging

from common.interpolator import interpolate
from pqr import pqr2 as pqr


def midunwrap(a, axis=-1):
    N = a.shape[axis]//2
    t = a.take(N, axis=axis)
    a = np.unwrap(a, axis=axis)
    s = np.array(a.shape)
    s[axis] = 1
    t -= a.take(N, axis=axis)
#    print(a.shape, axis, t.shape, s)
    a += t.reshape(s) #but we likely need some magic here
    return a


def midunwrap2d(a):
    return midunwrap(midunwrap(a, axis=0))


def _AmplToNoise(x):
    print("_AmplToNoise", x.shape)
    N = 2
    sr = np.empty((x.shape[1] - N, x.shape[2] - N), float)
    si = np.empty((x.shape[1] - N, x.shape[2] - N), float)
    S = []
    res = np.abs(x)
    for k in range(x.shape[0]):
        data = x[k]
        for i in range(sr.shape[0]):
            for j in range(sr.shape[1]):
                sr[i, j] = data.real[i:i + N, j:j + N].std()
                si[i, j] = data.imag[i:i + N, j:j + N].std()

        def filter(y, M=1):
            for i in range(M):
                c, edges = np.histogram(y, bins=20, range=(np.nanmin(y), np.nanmax(y)))
                ci = c.argmax()
                my = 0.5 * (edges[ci] + edges[ci + 1])
                sy = np.nanmean(y)
                y[y > my + 1.5 * sy] = np.NaN

        filter(sr, 5)
        filter(si, 5)

        s = (np.nanmean(sr)**2+np.nanmean(si)**2)**0.5
        print("_AmplToNoise", k, s)
        res[k] /= s
    return res


#TODO: unwrap phase from maximum amplitude
PAIR={'Real': lambda x: x.real,
      'Imag': lambda x: x.imag,
      'Ampl': lambda x: np.abs(x),
      'Phase': lambda x: np.angle(x),
      None:lambda x:x,
      "PhaseA": lambda x: (np.angle(x), np.abs(x)),
      "PhaseU": lambda x: midunwrap2d(np.angle(x)),
      "PhaseUA": lambda x: (midunwrap2d(np.angle(x)), np.abs(x)),
      "Ampl/Noise": _AmplToNoise,
      }


class Node:
    """
    Have data source.
    Have settable controls.
    If any of these two change, process the data by using the controls.
    Outputs processed data.
    """
    _required_controls = tuple()

    def __init__(self):
        super().__init__()
        self._input_data = None
        self._output_data = None
        self._sinks = []
        self._controls = {}

    def add_sink(self, sink):
        """
        sink will be called with processed data
        :param sink:
        :return:
        """
        self._sinks.append(sink)

    def set_input_data(self, data):
        logging.debug("node %s set_input_data()", self)
        self._input_data = data
        self.process()

    def set_controls(self, **kw):
        self._controls.update(**kw)
        self.process()

    def output_data(self):
        return self._output_data

    def _process(self, *args):
        """
        This is to be overloaded for derived classes. Return processed data.
        self._input_data and self._required_controls are guaranteed to be set (but not checked for acceptable value)
        :return:
        """
        return self._input_data

    def process(self):
        logging.debug("node %s process(): data ready %s, controls %s", self, self._input_data is not None, self._controls)
        if self._input_data is None:
            return

        logging.debug("\trequired controls %s", self._required_controls)
        for control in self._required_controls:
            if control not in self._controls:
                return

        self._output_data = self._process(*(self._controls[it] for it in self._required_controls))
        logging.debug("\tsinks %s", self._sinks)
        for sink in self._sinks:
            sink(self._output_data)

# the bad thing is that we cannot hint controls, they are hidden inside _process()
#  possibly can be automated by inspecting derived class's _process signature
#  as it is the order of _required_controls and arguments in _process has to be the same


def automation(cls):
    sig = inspect.signature(cls._process)
    cls._required_controls = tuple(sig.parameters.keys())[1:]  # skip self
    # todo: we can create and illusion of custom set_controls() with the same signature
    return cls

@automation
class SelectionNode(Node):
    def _process(self, pair, trn):
        w1, t2, w3, frames = self._input_data
        data = PAIR[pair](frames[trn])

        if pair == "PhaseA" or pair == "PhaseUA":
            data, temp = data
            # toto transformace by mohla byt trochu nelinearni, utlumit jen male maplitudy, kde se fazi neda verit
            # alpha = temp * (255/temp.reshape((temp.shape[0], -1)).max(axis=1).reshape((-1, 1, 1)))
            temp *= 255 / temp.reshape((temp.shape[0], -1)).max(axis=1).reshape((-1, 1, 1))
            # ~ alpha = (temp > 50) * 255
            # ~ del temp
            alpha = temp
        else:
            alpha = None

        # if pair.startswith("Phase"):
        #     # TODO: use colorcet linear perception colormaps here, allow user selection
        #     self.plotMain.plot().y3.colorMap = pqr.cmRGBCycle  # this should be colorAxis
        # else:
        #     self.plotMain.plot().y3.colorMap = pqr.cmDeepBlueRed

        # for dataRange we need to keep infs and nans out
        temp = data[np.isfinite(data)]
        if temp.size == 0:
            dataRange = (-1, 1)
        else:
            dataRange = (temp.min(), temp.max())  # global range

        return w1, t2, w3, data, dataRange, alpha


# TODO: there are actually other tricks in PAIR, see tabTemplate
@automation
class InterpolationNode(Node):
    def _process(self, interpolation1, interpolation3):
        w1, t2, w3, data, dataRange, alpha = self._input_data
        if interpolation1 > 0:
            data, w1 = interpolate(data, interpolation1, w1, axis=1)
        if interpolation3 > 0:
            data, w3 = interpolate(data, interpolation3, w3, axis=2)
        return w1, t2, w3, data, dataRange, alpha


class DDdataFigureController(QtWidgets.QWidget):
    def __init__(self, *args):
        super().__init__(*args)

        self.cb_trn = QtWidgets.QComboBox()
        self.cb_trn.addItems(("total", "rephasing", "nonrephasing"))
        self.cb_pair = QtWidgets.QComboBox()
        self.cb_pair.addItems(tuple(PAIR.keys()))
        self.check_global_colors = QtWidgets.QCheckBox("Global")
        self.check_fix_colors = QtWidgets.QCheckBox("Fix")

        def setup_spin():
            spin = QtWidgets.QSpinBox()
            spin.setMaximum(1000)
            spin.setSingleStep(100)
            return spin

        self.spin_interpolation1 = setup_spin()
        self.spin_interpolation1.setToolTip("Interpolation of ω1 axis.")
        self.spin_interpolation3 = setup_spin()
        self.spin_interpolation3.setToolTip("Interpolation of ω3 axis.")

        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(self.cb_trn)
        layout.addWidget(self.cb_pair)
        layout.addWidget(self.check_global_colors)
        layout.addWidget(self.check_fix_colors)
        layout.addWidget(QtWidgets.QLabel("Interpolation"))
        layout.addWidget(self.spin_interpolation1)
        layout.addWidget(self.spin_interpolation3)
        self._layout = layout

        self.setLayout(layout)

    def fix_colors(self):
        return self.check_fix_colors.isChecked()

    def global_colors(self):
        return self.check_global_colors.isChecked()


class DDdataFigure(pqr.FigureWidget):
    """

    Data processing sequence is

    raw_data
    selection TRN, PAIR
    interpolation

    """
    def __init__(self, controller):
        super().__init__()

        self._controller = controller

        plot = self.plot()
        xaxis = plot.addAxis(pqr.NumEnumTransform(), pqr.AxisPosition.B, label='ω₁ / cm⁻¹', tickFormat="{:.0f}")
        yaxis = plot.addAxis(pqr.NumEnumTransform(), pqr.AxisPosition.L, label='ω₃ / cm⁻¹', tickFormat="{:.0f}")
        xaxis2 = plot.addAxis(pqr.AOverTransform(1e7, xaxis._transform), pqr.AxisPosition.top, label="λ₁ / nm")
        yaxis2 = plot.addAxis(pqr.AOverTransform(1e7, yaxis._transform), pqr.AxisPosition.right, label="λ₃ / nm")

        xaxis.zoomRangeLChanged.connect(xaxis2._setZoomRange)
        xaxis2.zoomRangeLChanged.connect(xaxis._setZoomRange)
        yaxis.zoomRangeLChanged.connect(yaxis2._setZoomRange)
        yaxis2.zoomRangeLChanged.connect(yaxis._setZoomRange)

        plotter = plot.addDDPlotter(xaxis, yaxis)
        tracer = plot.addDecorator(pqr.TracerDecorator(plot, xaxis, yaxis))

        # data pipeline
        self._processed_data = None
        self._index = None

        self._input_node = Node()
        self._selection_node = SelectionNode()
        self._interpolation_node = InterpolationNode()

        self._input_node.add_sink(self._selection_node.set_input_data)
        self._selection_node.add_sink(self._interpolation_node.set_input_data)
        self._interpolation_node.add_sink(self._set_processed_data)

        self._controller.spin_interpolation1.valueChanged.connect(lambda x: self._interpolation_node.set_controls(interpolation1=x))
        self._controller.spin_interpolation3.valueChanged.connect(lambda x: self._interpolation_node.set_controls(interpolation3=x))

        self._controller.cb_trn.currentTextChanged.connect(lambda x: self._selection_node.set_controls(trn=x))
        self._controller.cb_pair.currentTextChanged.connect(lambda x: self._selection_node.set_controls(pair=x))

        self._controller.check_global_colors.stateChanged.connect(self._update_plot)

        self._interpolation_node.set_controls(interpolation1=self._controller.spin_interpolation1.value(),
                                              interpolation3=self._controller.spin_interpolation3.value())
        self._selection_node.set_controls(pair=self._controller.cb_pair.currentText(),
                                          trn=self._controller.cb_trn.currentText())


    def set_index(self, index):
        self._index = index
        self._update_plot()

    def set_data(self, data_frames, t2, w1, w3):
        """
        :param data_frames: {"total":total_data[t2, w1, w3], "rephasing":rephasing_data[t2, w1, w3], "non-rephasing":...}
        :param t2: might be complex for DAS, imaginary part gives the frequency
        :param w1:
        :param w3:
        :return:
        """
        self._input_node.set_input_data((w1, t2, w3, data_frames))

    def _set_processed_data(self, data):
        self._processed_data = data
        self._update_plot()

    def _update_plot(self):
        """
        new data or new frame selection
        :return:
        """
        logging.debug("data ready %s, index %s",

                      self._processed_data is not None, self._index)
        if self._processed_data is None:
            return

        if self._index is None:
            return

        w1, t2, w3, data, dataRange, alpha = self._processed_data

        if self._index<0 or self._index>=len(t2):
            return

        frame = data[self._index]

        if alpha is not None: alpha = alpha[self._index].T
        # ~ print("\t\tframe", frame)
        # ~ print("\t\tframe all nan", np.all(np.isnan(frame)))

        # ~ print("\t", index, frame)
        if self._controller.global_colors():
            a = max(abs(dataRange[0]), abs(dataRange[1]))
        else:
            # for dataRange we need to keep infs and nans out
            temp = frame[np.isfinite(frame)]
            if temp.size == 0:
                a = 1
            else:
                a = max(abs(temp.min()), abs(temp.max()))

        dataRange = (-a, a)
        # ~ print("\t\tdata range", -a, a)

        # TODO: contours do not work properly from some time (not sure what happend), they display at 0,0 not at axis coordinates
        plot = self.plot()
        plotter = plot.plotter(0)

        plotter._xaxis._transform.stops = w1
        plotter._yaxis._transform.stops = w3

        fix_colors = self._controller.fix_colors()
        if fix_colors:
            zoomDataRange = plotter._zaxis.linear2data(plotter._zaxis._zoomRange())

        plotter.setData(frame.T, w1, w3, dataRange=dataRange,
                        alpha=alpha)  # Todo: keep contours settings from previous (Traced2DPlot is originaly inteded for one time use, not interactive)
        # ~ self.plotMain.replot()
        plot.zoomToData()

        if fix_colors:
            plotter._zaxis._setZoomRange(*plotter._zaxis.data2linear(zoomDataRange))
