"""Starts the application."""
import os
import sys

from PyQt5 import QtWidgets

from soliton_automata.gui.main_window_ba import MainWindow


def main():
    """Initializes object of class `MainWindow`, sets style sheet of the window and then executes the application.
    """
    os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "1"
    sys.argv += ['--style', 'fusion']
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()