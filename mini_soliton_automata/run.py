import os
import sys

from PyQt5 import QtWidgets

from gui.main_window import MainWindow

if __name__ == "__main__":
    """Initializes object of class `MainWindow`, sets style sheet of the window and then executes the application.
    """
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    with open(os.path.join('mini_soliton_automata','gui', 'styles.css'), 'r') as f:
        style = f.read()
        app.setStyleSheet(style)
    sys.exit(app.exec_())
