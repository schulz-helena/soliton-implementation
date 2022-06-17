# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'main_window.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QMessageBox

from soliton_automata import SolitonAutomata
from soliton_graph import SolitonGraph
from visualisation import Visualisation


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(800, 600)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.display_molecule = QtWidgets.QLabel(self.centralwidget)
        self.display_molecule.setGeometry(QtCore.QRect(140, 10, 541, 381))
        self.display_molecule.setAutoFillBackground(False)
        self.display_molecule.setText("")
        self.display_molecule.setScaledContents(True)
        self.display_molecule.setObjectName("display_molecule")
        self.molecule_label = QtWidgets.QLabel(self.centralwidget)
        self.molecule_label.setGeometry(QtCore.QRect(140, 440, 71, 21))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.molecule_label.setFont(font)
        self.molecule_label.setObjectName("molecule_label")
        self.exterior_nodes_label = QtWidgets.QLabel(self.centralwidget)
        self.exterior_nodes_label.setGeometry(QtCore.QRect(140, 470, 111, 21))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.exterior_nodes_label.setFont(font)
        self.exterior_nodes_label.setObjectName("exterior_nodes_label")
        self.molecule_lineedit = QtWidgets.QLineEdit(self.centralwidget)
        self.molecule_lineedit.setGeometry(QtCore.QRect(260, 440, 311, 21))
        self.molecule_lineedit.setObjectName("molecule_lineedit")
        self.submit_molecule = QtWidgets.QPushButton(self.centralwidget)
        self.submit_molecule.setGeometry(QtCore.QRect(600, 440, 91, 31))
        self.submit_molecule.setObjectName("submit_molecule")
        self.node_1 = QtWidgets.QComboBox(self.centralwidget)
        self.node_1.setGeometry(QtCore.QRect(260, 470, 131, 26))
        self.node_1.setObjectName("node_1")
        self.node_2 = QtWidgets.QComboBox(self.centralwidget)
        self.node_2.setGeometry(QtCore.QRect(440, 470, 131, 26))
        self.node_2.setObjectName("node_2")
        self.submit_exterior_nodes = QtWidgets.QPushButton(self.centralwidget)
        self.submit_exterior_nodes.setGeometry(QtCore.QRect(600, 470, 91, 31))
        self.submit_exterior_nodes.setObjectName("submit_exterior_nodes")
        self.exterior_nodes_label_2 = QtWidgets.QLabel(self.centralwidget)
        self.exterior_nodes_label_2.setGeometry(QtCore.QRect(410, 470, 21, 21))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.exterior_nodes_label_2.setFont(font)
        self.exterior_nodes_label_2.setObjectName("exterior_nodes_label_2")
        self.welcome_label_1 = QtWidgets.QLabel(self.centralwidget)
        self.welcome_label_1.setGeometry(QtCore.QRect(130, 120, 551, 41))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.welcome_label_1.setFont(font)
        self.welcome_label_1.setAlignment(QtCore.Qt.AlignCenter)
        self.welcome_label_1.setObjectName("welcome_label_1")
        self.welcome_label_2 = QtWidgets.QLabel(self.centralwidget)
        self.welcome_label_2.setGeometry(QtCore.QRect(130, 150, 551, 41))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.welcome_label_2.setFont(font)
        self.welcome_label_2.setAlignment(QtCore.Qt.AlignCenter)
        self.welcome_label_2.setObjectName("welcome_label_2")
        self.save = QtWidgets.QPushButton(self.centralwidget)
        self.save.setGeometry(QtCore.QRect(610, 360, 71, 31)) # Save Button direkt unter Bild: 620, 391
        self.save.setObjectName("save")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 24))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

        # What I added:
        self.submit_molecule.clicked.connect(self.submit_molecule_clicked)
        self.submit_exterior_nodes.clicked.connect(self.submit_exterior_nodes_clicked)
        self.save.hide()
        self.exterior_nodes_label.hide()
        self.exterior_nodes_label_2.hide()
        self.node_1.hide()
        self.node_2.hide()
        self.submit_exterior_nodes.hide()

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Soliton Automata"))
        self.molecule_label.setText(_translate("MainWindow", "Molecule:"))
        self.exterior_nodes_label.setText(_translate("MainWindow", "Exterior nodes:"))
        self.submit_molecule.setText(_translate("MainWindow", "Submit"))
        self.submit_exterior_nodes.setText(_translate("MainWindow", "Submit"))
        self.exterior_nodes_label_2.setText(_translate("MainWindow", "&"))
        self.welcome_label_1.setText(_translate("MainWindow", "Welcome to the soliton automata software!"))
        self.welcome_label_2.setText(_translate("MainWindow", "Please specify your molecule below:"))
        self.save.setText(_translate("MainWindow", "Save"))

    def submit_molecule_clicked(self):
        self.node_1.clear()
        self.node_2.clear()
        smiles_string = self.molecule_lineedit.text()
        try:
            self.my_graph = SolitonGraph(smiles_string)
            errors = self.my_graph.validate_soliton_graph()
            Visualisation.visualize_soliton_graph(self.my_graph, self.my_graph.bindings, False)
            self.display_molecule.setPixmap(QtGui.QPixmap("database/graph.jpg"))
            self.welcome_label_1.hide()
            self.welcome_label_2.hide()
            if errors != []:
                self.save.hide()
                self.exterior_nodes_label.hide()
                self.exterior_nodes_label_2.hide()
                self.node_1.hide()
                self.node_2.hide()
                self.submit_exterior_nodes.hide()
                msg = QMessageBox()
                msg.setWindowTitle("No soliton graph")
                msg.setText("You specified a molecule that does not fulfill the requirements of a soliton graph")
                msg.setIcon(QMessageBox.Warning)
                msg.setStandardButtons(QMessageBox.Retry)
                msg.setInformativeText("See details for all incorrect parts of your molecule.")
                details = ""
                for error in errors:
                    details = details + error + "\n"
                msg.setDetailedText(details)
                x = msg.exec_() # show messagebox
            else:
                self.save.show()
                self.exterior_nodes_label.show()
                self.exterior_nodes_label_2.show()
                self.node_1.show()
                self.node_2.show()
                self.submit_exterior_nodes.show()
                for key in self.my_graph.exterior_nodes_reverse:
                    self.node_1.addItem(str(key))
                    self.node_2.addItem(str(key))
        except:
            self.save.hide()
            self.exterior_nodes_label.hide()
            self.exterior_nodes_label_2.hide()
            self.node_1.hide()
            self.node_2.hide()
            self.submit_exterior_nodes.hide()
            msg = QMessageBox()
            msg.setWindowTitle("Incorrect input")
            msg.setText("The syntax of your input is not correct")
            msg.setIcon(QMessageBox.Warning)
            msg.setStandardButtons(QMessageBox.Retry)
            msg.setInformativeText("Please try again with another input string.")
            msg.setDetailedText("Some details..") #TODO: Write a more detailed text to help user
            x = msg.exec_()

    def submit_exterior_nodes_clicked(self):
        node1 = int(self.node_1.currentText())
        node2 = int(self.node_2.currentText())
        if node1 != node2:
            automata = SolitonAutomata(self.my_graph, node1, node2)
            if automata.paths == []:
                msg = QMessageBox()
                msg.setWindowTitle("No path found")
                msg.setText("There exists no soliton path between these exterior nodes")
                msg.setIcon(QMessageBox.Information)
                msg.setStandardButtons(QMessageBox.Retry)
                msg.setInformativeText("Please try again with different exterior nodes.")
                x = msg.exec_()
            else:
                print(automata.paths) # TODO: instead of printing open new window with found paths
        else:
            msg = QMessageBox()
            msg.setWindowTitle("Equal nodes")
            msg.setText("You chose the same exterior node twice")
            msg.setIcon(QMessageBox.Warning)
            msg.setStandardButtons(QMessageBox.Retry)
            msg.setInformativeText("Please try again with two differing nodes.")
            x = msg.exec_()


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

