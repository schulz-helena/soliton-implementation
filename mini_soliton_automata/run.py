import os
import sys

from PyQt5 import QtWidgets

from mini_soliton_automata.gui.main_window import MainWindow

#from pathlib import Path




def main():
    """Initializes object of class `MainWindow`, sets style sheet of the window and then executes the application.
    """
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    readme_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'styles.css')
    with open(readme_path, "r", encoding="utf-8") as fh:
        style = fh.read()
    #with open(os.path.join('styles.css'), 'r') as f:
        #style = f.read()
    #this_directory = Path(__file__).parent
    #style = (this_directory / 'styles.css').read_text()
    app.setStyleSheet(style)
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
