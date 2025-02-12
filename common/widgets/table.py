from PyQt5 import QtWidgets, QtGui


class Table(QtWidgets.QTableView):
    """
    QTableView with multi-selection copy/paste.
    """
    def keyPressEvent(self, event: QtGui.QKeyEvent) -> None:
        if event == QtGui.QKeySequence.Copy:
            indexes = self.selectionModel().selectedIndexes()

            rows = set()
            cols = set()

            for index in indexes:
                rows.add(index.row())
                cols.add(index.column())

            min_col = min(cols)
            max_col = max(cols)
            min_row = min(rows)
            max_row = max(rows)

            columns = max_col - min_col + 1
            rows = max_row - min_row + 1

            text_table = [[""] * columns for i in range(rows)]

            for index in indexes:
                text_table[index.row()-min_row][index.column()-min_col] = "{}".format(self.model().data(index))

            text = "\n".join(("\t".join(i) for i in text_table))

            QtGui.QGuiApplication.clipboard().setText(text)
            event.accept()
        else:
            super().keyPressEvent(event)
