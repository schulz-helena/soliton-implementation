# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'main_window2.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

import re

from PIL import Image
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QDialog, QMessageBox, QScrollArea

from animation import Animation
from soliton_automata import SolitonAutomata
from soliton_graph import SolitonGraph
from soliton_path import SolitonPath
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
        self.molecule_label.setGeometry(QtCore.QRect(130, 420, 71, 21))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.molecule_label.setFont(font)
        self.molecule_label.setObjectName("molecule_label")
        self.exterior_nodes_label = QtWidgets.QLabel(self.centralwidget)
        self.exterior_nodes_label.setGeometry(QtCore.QRect(130, 450, 111, 21))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.exterior_nodes_label.setFont(font)
        self.exterior_nodes_label.setObjectName("exterior_nodes_label")
        self.molecule_lineedit = QtWidgets.QLineEdit(self.centralwidget)
        self.molecule_lineedit.setGeometry(QtCore.QRect(260, 420, 311, 21))
        self.molecule_lineedit.setObjectName("molecule_lineedit")
        self.submit_molecule = QtWidgets.QPushButton(self.centralwidget)
        self.submit_molecule.setGeometry(QtCore.QRect(600, 420, 91, 31))
        self.submit_molecule.setObjectName("submit_molecule")
        self.node_1 = QtWidgets.QComboBox(self.centralwidget)
        self.node_1.setGeometry(QtCore.QRect(260, 450, 131, 26))
        self.node_1.setObjectName("node_1")
        self.node_2 = QtWidgets.QComboBox(self.centralwidget)
        self.node_2.setGeometry(QtCore.QRect(440, 450, 131, 26))
        self.node_2.setObjectName("node_2")
        self.submit_exterior_nodes = QtWidgets.QPushButton(self.centralwidget)
        self.submit_exterior_nodes.setGeometry(QtCore.QRect(600, 450, 91, 31))
        self.submit_exterior_nodes.setObjectName("submit_exterior_nodes")
        self.exterior_nodes_label_2 = QtWidgets.QLabel(self.centralwidget)
        self.exterior_nodes_label_2.setGeometry(QtCore.QRect(410, 450, 21, 21))
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
        self.save.setGeometry(QtCore.QRect(610, 360, 71, 31))
        self.save.setObjectName("save")
        self.soliton_paths_label = QtWidgets.QLabel(self.centralwidget)
        self.soliton_paths_label.setGeometry(QtCore.QRect(130, 480, 131, 21))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.soliton_paths_label.setFont(font)
        self.soliton_paths_label.setObjectName("soliton_paths_label")
        self.paths = QtWidgets.QComboBox(self.centralwidget)
        self.paths.setGeometry(QtCore.QRect(260, 480, 431, 26))
        self.paths.setObjectName("paths")
        self.show_matrices = QtWidgets.QPushButton(self.centralwidget)
        self.show_matrices.setGeometry(QtCore.QRect(260, 510, 131, 31))
        self.show_matrices.setObjectName("show_matrices")
        self.show_end_result = QtWidgets.QPushButton(self.centralwidget)
        self.show_end_result.setGeometry(QtCore.QRect(400, 510, 141, 31))
        self.show_end_result.setObjectName("show_end_result")
        self.show_animation = QtWidgets.QPushButton(self.centralwidget)
        self.show_animation.setGeometry(QtCore.QRect(550, 510, 141, 31))
        self.show_animation.setObjectName("show_animation")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 24))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

        # What I added:
        self.submit_molecule.clicked.connect(self.submit_molecule_clicked)
        self.submit_exterior_nodes.clicked.connect(self.submit_exterior_nodes_clicked)
        self.show_matrices.clicked.connect(self.show_matrices_clicked)
        self.show_end_result.clicked.connect(self.show_end_result_clicked)
        self.show_animation.clicked.connect(self.show_animation_clicked)
        self.save.clicked.connect(self.save_clicked)
        self.save.hide()
        self.exterior_nodes_label.hide()
        self.exterior_nodes_label_2.hide()
        self.node_1.hide()
        self.node_2.hide()
        self.submit_exterior_nodes.hide()
        self.soliton_paths_label.hide()
        self.paths.hide()
        self.show_matrices.hide()
        self.show_end_result.hide()
        self.show_animation.hide()

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
        self.soliton_paths_label.setText(_translate("MainWindow", "Soliton paths:"))
        self.show_matrices.setText(_translate("MainWindow", "Show matrices"))
        self.show_end_result.setText(_translate("MainWindow", "Show end result"))
        self.show_animation.setText(_translate("MainWindow", "Show animation"))

    def submit_molecule_clicked(self):
        self.node_1.clear()
        self.node_2.clear()
        smiles_string = self.molecule_lineedit.text()
        try:
            self.my_graph = SolitonGraph(smiles_string)
            errors = self.my_graph.validate_soliton_graph()
            Visualisation.visualize_soliton_graph(self.my_graph, self.my_graph.bindings, False, "graph")
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
                self.soliton_paths_label.hide()
                self.paths.hide()
                self.show_matrices.hide()
                self.show_end_result.hide()
                self.show_animation.hide()
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
                self.soliton_paths_label.hide()
                self.paths.hide()
                self.show_matrices.hide()
                self.show_end_result.hide()
                self.show_animation.hide()
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
            self.soliton_paths_label.hide()
            self.paths.hide()
            self.show_matrices.hide()
            self.show_end_result.hide()
            self.show_animation.hide()
            msg = QMessageBox()
            msg.setWindowTitle("Incorrect input")
            msg.setText("The syntax of your input is not correct")
            msg.setIcon(QMessageBox.Warning)
            msg.setStandardButtons(QMessageBox.Retry)
            msg.setInformativeText("Please try again with another input string.")
            msg.setDetailedText("Some details..") #TODO: Write a more detailed text to help user
            x = msg.exec_()

    def save_clicked(self):
        option = QtWidgets.QFileDialog.Options()
        name = QtWidgets.QFileDialog.getSaveFileName(self.centralwidget, 'Save File', 'graph.jpg', options = option)
        path = name[0]
        file = open('database/graph.jpg', 'rb')
        data = file.read()
        file.close()
        file = open(path, "wb")
        file.write(data)
        file.close()

    def submit_exterior_nodes_clicked(self):
        self.path_index = None # we need this variable later in show_animation_clicked
        self.paths.clear()
        node1 = int(self.node_1.currentText())
        node2 = int(self.node_2.currentText())
        if node1 != node2:
            self.automata = SolitonAutomata(self.my_graph, node1, node2)
            if self.automata.paths == []:
                self.soliton_paths_label.hide()
                self.paths.hide()
                self.show_matrices.hide()
                self.show_end_result.hide()
                self.show_animation.hide()
                msg = QMessageBox()
                msg.setWindowTitle("No path found")
                msg.setText("There exists no soliton path between these exterior nodes")
                msg.setIcon(QMessageBox.Information)
                msg.setStandardButtons(QMessageBox.Retry)
                msg.setInformativeText("Please try again with different exterior nodes.")
                x = msg.exec_()
            else:
                self.soliton_paths_label.show()
                self.paths.show()
                self.show_matrices.show()
                self.show_end_result.show()
                self.show_animation.show()
                for path in self.automata.paths_for_user:
                    self.paths.addItem(str(path))
                #print(self.automata.paths)
        else:
            self.soliton_paths_label.hide()
            self.paths.hide()
            self.show_matrices.hide()
            self.show_end_result.hide()
            self.show_animation.hide()
            msg = QMessageBox()
            msg.setWindowTitle("Equal nodes")
            msg.setText("You chose the same exterior node twice")
            msg.setIcon(QMessageBox.Warning)
            msg.setStandardButtons(QMessageBox.Retry)
            msg.setInformativeText("Please try again with two differing nodes.")
            x = msg.exec_()
    
    def show_matrices_clicked(self):

        def save_matrices():
            option = QtWidgets.QFileDialog.Options()
            name = QtWidgets.QFileDialog.getSaveFileName(self.centralwidget, 'Save File', 'matrices.txt', options = option)
            path = name[0]
            file = open('database/matrices.txt', 'rb')
            data = file.read()
            file.close()
            file = open(path, "wb")
            file.write(data)
            file.close()

        index = self.paths.currentIndex()
        desired_path = SolitonPath(self.automata.paths_ids[index], self.my_graph)
    
        dlg = QDialog()
        scrollArea = QScrollArea(dlg)
        widget = QtWidgets.QWidget()
        vbox = QtWidgets.QHBoxLayout()
        file = open("database/matrices.txt", 'w')
        for object in desired_path.adjacency_matrices_list:
            matrix = str(object)
            matrix = re.sub(r"[matrix(]", "", matrix)
            matrix = re.sub(r"[)]", "", matrix)
            file.write(matrix)
            file.write('\n')
            label = QtWidgets.QLabel(matrix)
            vbox.addWidget(label)
        file.close()
        widget.setLayout(vbox)
        scrollArea.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        scrollArea.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        #scrollArea.setWidgetResizable(True)
        scrollArea.setWidget(widget)
        scrollArea.setFixedSize(540, 380)
        save_button = QtWidgets.QPushButton("Save", dlg)
        save_button.setGeometry(QtCore.QRect(450, 330, 70, 30))
        save_button.clicked.connect(save_matrices)

        dlg.setWindowTitle("Adjacency Matrices")
        dlg.setFixedSize(540, 380)
        dlg.exec_()
    
    def show_end_result_clicked(self):

        def save_end_result():
            option = QtWidgets.QFileDialog.Options()
            name = QtWidgets.QFileDialog.getSaveFileName(self.centralwidget, 'Save File', 'result.jpg', options = option)
            path = name[0]
            file = open('database/result.jpg', 'rb')
            data = file.read()
            file.close()
            file = open(path, "wb")
            file.write(data)
            file.close()

        index = self.paths.currentIndex()
        desired_path = SolitonPath(self.automata.paths_ids[index], self.my_graph)
        bindings_index = len(desired_path.path) - 1
        Visualisation.visualize_soliton_graph(self.my_graph, desired_path.bindings_list[bindings_index], False, "result")
        #Animation.graph_animation(self.my_graph, desired_path)

        dlg = QDialog()
        label = QtWidgets.QLabel(dlg)
        label.setGeometry(QtCore.QRect(0, 0, 540, 380))
        label.setPixmap(QtGui.QPixmap("database/result.jpg"))
        label.setScaledContents(True)
        save_button = QtWidgets.QPushButton("Save", dlg)
        save_button.setGeometry(QtCore.QRect(470, 350, 70, 30))
        save_button.clicked.connect(save_end_result)

        dlg.setWindowTitle("End result")
        dlg.setFixedSize(540, 380)
        dlg.exec_()

    def show_animation_clicked(self):
        
        def save_animation():
            option = QtWidgets.QFileDialog.Options()
            name = QtWidgets.QFileDialog.getSaveFileName(self.centralwidget, 'Save File', 'animation.gif', options = option)
            path = name[0]
            file = open('database/animation.gif', 'rb')
            data = file.read()
            file.close()
            file = open(path, "wb")
            file.write(data)
            file.close()

        if self.path_index != self.paths.currentIndex():
            self.path_index = self.paths.currentIndex()
            desired_path = SolitonPath(self.automata.paths_ids[self.path_index], self.my_graph)
            Animation.graph_animation(self.my_graph, desired_path)

        dlg = QDialog()
        label = QtWidgets.QLabel(dlg)
        label.setGeometry(QtCore.QRect(0, 0, 540, 380))
        movie = QtGui.QMovie('database/animation.gif')
        movie.setScaledSize(QtCore.QSize(540,380))
        label.setMovie(movie)
        movie.start()
        save_button = QtWidgets.QPushButton("Save", dlg)
        save_button.setGeometry(QtCore.QRect(470, 350, 70, 30))
        save_button.clicked.connect(save_animation)

        dlg.setWindowTitle("Animation")
        dlg.setFixedSize(540, 380)
        dlg.exec_()


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

