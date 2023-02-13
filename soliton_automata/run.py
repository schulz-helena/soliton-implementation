"""Starts the application."""
import os
import sys

from gui.main_window_ba import MainWindow
from PyQt5 import QtWidgets

<<<<<<< HEAD:soliton_automata/run.py
=======
from mini_soliton_automata.gui.main_window import MainWindow

>>>>>>> d860cea084a9aebbbbb7efc487a72e65e1d06619:mini_soliton_automata/run.py

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