# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'main_window3.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

import math
import re

from PIL import Image
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QDialog, QMainWindow, QMessageBox, QScrollArea

from animation import Animation
from soliton_automata import SolitonAutomata
from soliton_graph import SolitonGraph
from soliton_path import SolitonPath
from visualisation import Visualisation


class Ui_MainWindow(QMainWindow):
    def __init__(self):
        super(Ui_MainWindow, self).__init__()
        self.setObjectName("MainWindow")
        self.resize(564, 599)
        self.centralwidget = QtWidgets.QWidget(self)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.show_animation = QtWidgets.QPushButton(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.show_animation.sizePolicy().hasHeightForWidth())
        self.show_animation.setSizePolicy(sizePolicy)
        self.show_animation.setObjectName("show_animation")
        self.gridLayout.addWidget(self.show_animation, 6, 4, 1, 1)
        self.exterior_nodes_label2 = QtWidgets.QLabel(self.centralwidget)
        self.exterior_nodes_label2.setObjectName("exterior_nodes_label2")
        self.gridLayout.addWidget(self.exterior_nodes_label2, 4, 2, 1, 1)
        self.show_end_result = QtWidgets.QPushButton(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.show_end_result.sizePolicy().hasHeightForWidth())
        self.show_end_result.setSizePolicy(sizePolicy)
        self.show_end_result.setObjectName("show_end_result")
        self.gridLayout.addWidget(self.show_end_result, 6, 2, 1, 2)
        self.save = QtWidgets.QPushButton(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.save.sizePolicy().hasHeightForWidth())
        self.save.setSizePolicy(sizePolicy)
        self.save.setMaximumSize(QtCore.QSize(16777215, 32))
        self.save.setObjectName("save")
        self.gridLayout.addWidget(self.save, 2, 4, 1, 1)
        self.molecule_lineedit = QtWidgets.QLineEdit(self.centralwidget)
        self.molecule_lineedit.setObjectName("molecule_lineedit")
        self.gridLayout.addWidget(self.molecule_lineedit, 3, 1, 1, 3)
        self.show_matrices = QtWidgets.QPushButton(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.show_matrices.sizePolicy().hasHeightForWidth())
        self.show_matrices.setSizePolicy(sizePolicy)
        self.show_matrices.setObjectName("show_matrices")
        self.gridLayout.addWidget(self.show_matrices, 6, 1, 1, 1)
        self.molecule_label = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(14)
        self.molecule_label.setFont(font)
        self.molecule_label.setObjectName("molecule_label")
        self.gridLayout.addWidget(self.molecule_label, 3, 0, 1, 1)
        self.node_2 = QtWidgets.QComboBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.node_2.sizePolicy().hasHeightForWidth())
        self.node_2.setSizePolicy(sizePolicy)
        self.node_2.setObjectName("node_2")
        self.gridLayout.addWidget(self.node_2, 4, 3, 1, 1)
        self.soliton_paths_label = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(14)
        self.soliton_paths_label.setFont(font)
        self.soliton_paths_label.setObjectName("soliton_paths_label")
        self.gridLayout.addWidget(self.soliton_paths_label, 5, 0, 1, 1)
        self.paths = QtWidgets.QComboBox(self.centralwidget)
        self.paths.setObjectName("paths")
        self.gridLayout.addWidget(self.paths, 5, 1, 1, 4)
        self.submit_exterior_nodes = QtWidgets.QPushButton(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.submit_exterior_nodes.sizePolicy().hasHeightForWidth())
        self.submit_exterior_nodes.setSizePolicy(sizePolicy)
        self.submit_exterior_nodes.setObjectName("submit_exterior_nodes")
        self.gridLayout.addWidget(self.submit_exterior_nodes, 4, 4, 1, 1)
        self.node_1 = QtWidgets.QComboBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.node_1.sizePolicy().hasHeightForWidth())
        self.node_1.setSizePolicy(sizePolicy)
        self.node_1.setObjectName("node_1")
        self.gridLayout.addWidget(self.node_1, 4, 1, 1, 1)
        self.exterior_nodes_label = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(14)
        self.exterior_nodes_label.setFont(font)
        self.exterior_nodes_label.setObjectName("exterior_nodes_label")
        self.gridLayout.addWidget(self.exterior_nodes_label, 4, 0, 1, 1)
        self.submit_molecule = QtWidgets.QPushButton(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.submit_molecule.sizePolicy().hasHeightForWidth())
        self.submit_molecule.setSizePolicy(sizePolicy)
        self.submit_molecule.setObjectName("submit_molecule")
        self.gridLayout.addWidget(self.submit_molecule, 3, 4, 1, 1)
        self.display_molecule = QtWidgets.QLabel(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Ignored)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.display_molecule.sizePolicy().hasHeightForWidth())
        self.display_molecule.setSizePolicy(sizePolicy)
        self.display_molecule.setMinimumSize(QtCore.QSize(540, 380))
        self.display_molecule.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.display_molecule.setSizeIncrement(QtCore.QSize(0, 0))
        self.display_molecule.setBaseSize(QtCore.QSize(0, 0))
        self.display_molecule.setAutoFillBackground(False)
        self.display_molecule.setText("")
        self.display_molecule.setPixmap(QtGui.QPixmap("database/startscreen.jpg"))
        #self.display_molecule.setPixmap(QtGui.QPixmap("database/startscreen.jpg").scaled(self.display_molecule.size(), QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation)) #.scaled(540 ,380, QtCore.Qt.AspectRatioMode.KeepAspectRatio)
        self.display_molecule.setScaledContents(True)
        self.display_molecule.setAutoFillBackground(True)
        p = self.display_molecule.palette()
        p.setColor(self.display_molecule.backgroundRole(), QtCore.Qt.red)
        self.display_molecule.setPalette(p)
        #self.display_molecule.setFixedSize(0,0)
        self.display_molecule.setObjectName("display_molecule")
        self.gridLayout.addWidget(self.display_molecule, 0, 0, 2, 5, QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter) # !!!
        self.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(self)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 564, 24))
        self.menubar.setObjectName("menubar")
        self.setMenuBar(self.menubar)

        self.retranslateUi()
        QtCore.QMetaObject.connectSlotsByName(self)

        self._heightForWidthFactor = 1.0 * 380 / 540
        #self._resizeImage(self.centralwidget.size())

        # What I added:
        self.submit_molecule.clicked.connect(self.submit_molecule_clicked)
        self.submit_exterior_nodes.clicked.connect(self.submit_exterior_nodes_clicked)
        self.show_matrices.clicked.connect(self.show_matrices_clicked)
        self.show_end_result.clicked.connect(self.show_end_result_clicked)
        self.show_animation.clicked.connect(self.show_animation_clicked)
        self.save.clicked.connect(self.save_clicked)

        self.hide_retain_space(self.save)
        self.save.hide()
        self.hide_retain_space(self.exterior_nodes_label)
        self.exterior_nodes_label.hide()
        self.hide_retain_space(self.exterior_nodes_label2)
        self.exterior_nodes_label2.hide()
        self.hide_retain_space(self.node_1)
        self.node_1.hide()
        self.hide_retain_space(self.node_2)
        self.node_2.hide()
        self.hide_retain_space(self.submit_exterior_nodes)
        self.submit_exterior_nodes.hide()
        self.hide_retain_space(self.soliton_paths_label)
        self.soliton_paths_label.hide()
        self.hide_retain_space(self.paths)
        self.paths.hide()
        self.hide_retain_space(self.show_matrices)
        self.show_matrices.hide()
        self.hide_retain_space(self.show_end_result)
        self.show_end_result.hide()
        self.hide_retain_space(self.show_animation)
        self.show_animation.hide()


    def retranslateUi(self):
        _translate = QtCore.QCoreApplication.translate
        self.setWindowTitle(_translate("MainWindow", "Soliton Automata"))
        self.show_animation.setText(_translate("MainWindow", "Show animation"))
        self.exterior_nodes_label2.setText(_translate("MainWindow", "&"))
        self.show_end_result.setText(_translate("MainWindow", "Show end result"))
        self.save.setText(_translate("MainWindow", "Save"))
        self.show_matrices.setText(_translate("MainWindow", "Show matrices"))
        self.molecule_label.setText(_translate("MainWindow", "Molecule:"))
        self.soliton_paths_label.setText(_translate("MainWindow", "Soliton paths:"))
        self.submit_exterior_nodes.setText(_translate("MainWindow", "Submit"))
        self.exterior_nodes_label.setText(_translate("MainWindow", "Exterior nodes:"))
        self.submit_molecule.setText(_translate("MainWindow", "Submit"))

    def submit_molecule_clicked(self):
        self.node_1.clear()
        self.node_2.clear()
        smiles_string = self.molecule_lineedit.text()
        try:
            self.my_graph = SolitonGraph(smiles_string)
            errors = self.my_graph.validate_soliton_graph()
            Visualisation.visualize_soliton_graph(self.my_graph, self.my_graph.bindings, False, "graph")
            self.display_molecule.setPixmap(QtGui.QPixmap("database/graph.jpg"))
            if errors != []:
                self.save.hide()
                self.exterior_nodes_label.hide()
                self.exterior_nodes_label2.hide()
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
                self.exterior_nodes_label2.show()
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
            self.exterior_nodes_label2.hide()
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

    def hide_retain_space(self, widget):
        retain = widget.sizePolicy()
        retain.setRetainSizeWhenHidden(True)
        widget.setSizePolicy(retain)


    def hasHeightForWidth(self):
        # This tells the layout manager that the banner's height does depend on its width
        return True

    def heightForWidth(self, width):
        # This tells the layout manager what the preferred and minimum height are, for a given width
        return math.ceil(width * self._heightForWidthFactor)

    def resizeEvent(self, event: QtGui.QResizeEvent):
        super(Ui_MainWindow, self).resizeEvent(event)
        # For efficiency, we pass the size from the event to _resizeImage()
        self._resizeImage(event.size())
        #self._resizeImage(self.centralwidget.size())
    
    '''def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        old_size = event.oldSize()
        new_size = QtCore.QSize(self.geometry().width(), self.geometry().height())
        print("old_size = {0}, new_size = {1}".format(old_size, new_size))
        QMainWindow.resizeEvent(self, event)
        width = self.centralwidget.size().width()
        height = self.heightForWidth(width)
        self.display_molecule.setFixedSize(width, height)'''

    def _resizeImage(self, size):
        # Since we're keeping _heightForWidthFactor, we can code a more efficient implementation of this, too
        width = size.width()
        height = self.heightForWidth(width)
        #self.display_molecule.setFixedSize(width, height)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    window = Ui_MainWindow()
    window.show()
    sys.exit(app.exec_())

