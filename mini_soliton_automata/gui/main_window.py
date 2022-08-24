"""GUI of the software.
 """
import io
import math
import re

import mini_soliton_automata.res.resources
from mini_soliton_automata.gui.startscreen import Startscreen
from mini_soliton_automata.soliton_classes.soliton_automata import \
    MiniSolitonAutomata
from mini_soliton_automata.soliton_classes.soliton_graph import SolitonGraph
from mini_soliton_automata.soliton_classes.soliton_path import SolitonPath
from mini_soliton_automata.visualisations.animation import Animation
from mini_soliton_automata.visualisations.visualisation import Visualisation
from PIL.ImageQt import ImageQt
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QDialog, QMainWindow, QMessageBox, QScrollArea


class MainWindow(QMainWindow):
    """Main Window of the GUI. Inherits from class `QMainWindow`.
    """
    def __init__(self):
        """Initializes the main window.
        Displays a welcoming text and all necessary widgets for the user to specify and submit a molecule.
        All other widgets are hidden for now and are revealed step by step, so user is guided through the use of the application.
        """
        super(MainWindow, self).__init__()
        self.setObjectName("MainWindow")
        # move window to the top of the screen + to the center horizontally
        qt_rectangle = self.frameGeometry()
        center_point = QtWidgets.QDesktopWidget().availableGeometry().center()
        qt_rectangle.moveCenter(center_point)
        self.move(qt_rectangle.topLeft().x(), 0)

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
        self.show_animation.setMinimumSize(QtCore.QSize(0, 32))
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
        self.show_end_result.setMinimumSize(QtCore.QSize(0, 32))
        self.show_end_result.setObjectName("show_end_result")
        self.gridLayout.addWidget(self.show_end_result, 6, 2, 1, 2)
        self.save = QtWidgets.QPushButton(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.save.sizePolicy().hasHeightForWidth())
        self.save.setSizePolicy(sizePolicy)
        self.save.setMaximumSize(QtCore.QSize(16777215, 32))
        self.save.setMinimumSize(QtCore.QSize(0, 32))
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
        self.show_matrices.setMinimumSize(QtCore.QSize(0, 32))
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
        self.submit_exterior_nodes.setMinimumSize(QtCore.QSize(0, 32))
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
        self.submit_molecule.setMinimumSize(QtCore.QSize(0, 32))
        self.submit_molecule.setObjectName("submit_molecule")
        self.gridLayout.addWidget(self.submit_molecule, 3, 4, 1, 1)

        self.display_molecule = QtWidgets.QLabel(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Ignored)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setWidthForHeight(self.display_molecule.sizePolicy().hasWidthForHeight())
        self.display_molecule.setSizePolicy(sizePolicy)
        self.display_molecule.setMinimumSize(QtCore.QSize(600, 450))
        self.display_molecule.setSizeIncrement(QtCore.QSize(0, 0))
        self.display_molecule.setBaseSize(QtCore.QSize(0, 0))
        self.display_molecule.setText("")
        startscreen = Startscreen().image
        self.qim = ImageQt(startscreen)
        self.display_molecule.setPixmap(QtGui.QPixmap.fromImage(self.qim).scaled(self.display_molecule.size(), QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation))
        self.display_molecule.setObjectName("display_molecule")
        self.gridLayout.addWidget(self.display_molecule, 0, 0, 2, 5, alignment = QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
        self.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(self)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 564, 24))
        self.menubar.setObjectName("menubar")
        self.setMenuBar(self.menubar)

        self.retranslateUi()
        QtCore.QMetaObject.connectSlotsByName(self)

        self.height_for_width_factor = 1.0 * 519 / 688
        self.width_for_height_factor = 1.0 * 688 / 519

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
        """Implements multi-language suppport. Is generated automatically when using PyQt5 UI code generator.
        """
        _translate = QtCore.QCoreApplication.translate
        self.setWindowTitle(_translate("MainWindow", "Mini Soliton Automata Software"))
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
        """Method that is called when user clicks button to submit the specified molecule.
        Catches errors if user used the wrong syntax or specified a molecule that does not fulfill the requirements of a soliton graph.
        If the user's molecule is valid it displays the graph of the molecule. It then also reveals a save button for the graph visualisation and
        all the necessary widgets for the user to choose a pair of exterior nodes.
        """
        self.node_1.clear()
        self.node_2.clear()
        self.smiles_string = self.molecule_lineedit.text()
        try:
            self.my_graph = SolitonGraph(self.smiles_string)
            errors = self.my_graph.validate_soliton_graph()
            self.graph_pic = Visualisation.visualize_soliton_graph(self.my_graph, self.my_graph.bindings, False, True)
            self.qim = ImageQt(self.graph_pic)
            self.display_molecule.setPixmap(QtGui.QPixmap.fromImage(self.qim).scaled(self.display_molecule.size(), QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation))
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
                msg.setStyleSheet(" QPushButton{ height: 32px; width: 130px;}")
                details = ""
                for i, error in enumerate(errors):
                    if i == len(errors) - 1:
                        details = details + error
                    else:
                        details = details + error + "\n"
                msg.setDetailedText(details)
                x = msg.exec_() # show messagebox
            elif self.my_graph.exterior_nodes_name_collision() == True:
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
                msg.setWindowTitle("Name collision")
                msg.setText("You specified two or more exterior nodes with the same name")
                msg.setStyleSheet(" QPushButton{ height: 32px; width: 130px}")
                msg.setIcon(QMessageBox.Warning)
                msg.setStandardButtons(QMessageBox.Retry)
                x = msg.exec_()
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
                    self.node_1.addItem(key)
                    self.node_2.addItem(key)
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
            msg.setStyleSheet(" QPushButton{ height: 32px; width: 130px}")
            details = f"Reminder - this is how you define a molecule: \n"
            details = details + f"- Carbon atoms are marked with 'C' \n"
            details = details + f"- Single edges are marked with '-' or no character at all \n"
            details = details + f"- Double edges are marked with '=' \n"
            details = details + f"- Branches are embedded in round brackets (e.g. 'C=(CC=C)C')\n"
            details = details + f"- The two connecting atoms of a ring are marked with the same number (e.g. 'C1' and 'C1') \n"
            details = details + "- Exterior nodes are marked with braces and a number (e.g. '{1}')"
            msg.setDetailedText(details)
            x = msg.exec_()


    def save_clicked(self):
        """Method that is called when user clicks button to save the graph visualisation.
        Opens a file dialog in which user can specify a path where the image should be saved.
        Only allows `.jpg`, `.png` and `.jpeg` file suffixes.
        """
        option = QtWidgets.QFileDialog.Options()
        name = QtWidgets.QFileDialog.getSaveFileName(self.centralwidget, 'Save File', 'graph.jpg', 'Images (*.jpg *.png *.jpeg)', options = option)
        if name != ('', ''):
            path = name[0]
            # turn PIL Image into byte array 
            imgByteArr = io.BytesIO()
            self.graph_pic.save(imgByteArr, format='PNG')
            imgByteArr = imgByteArr.getvalue()
            # save image at specified path
            file = open(path, "wb")
            file.write(imgByteArr) # before: file.write(data)
            file.close()


    def submit_exterior_nodes_clicked(self):
        """Method that is called when user clicks button to submit a pair of exterior nodes.
        Computes all possible soliton paths between the two chosen nodes.
        Informs the user if no soliton path exists between them.
        Otherwise all the necessary widgets for the user to choose a computed soliton path and look at further information on it are revealed.
        """
        self.path_index = None # we need this variable later in show_animation_clicked
        self.paths.clear()
        node1 = int(self.node_1.currentText())
        node2 = int(self.node_2.currentText())
        self.automata = MiniSolitonAutomata(self.my_graph, node1, node2)
        if self.automata.soliton_paths == []:
            self.soliton_paths_label.hide()
            self.paths.hide()
            self.show_matrices.hide()
            self.show_end_result.hide()
            self.show_animation.hide()
            msg = QMessageBox()
            msg.setWindowTitle("No path found")
            msg.setText("There exists no soliton path between these exterior nodes.")
            msg.setIcon(QMessageBox.Information)
            msg.setStandardButtons(QMessageBox.Retry)
            msg.setInformativeText("Please try again with different exterior nodes.")
            msg.setStyleSheet(" QPushButton{ height: 32px; width: 130px}")
            x = msg.exec_()
        else:
            self.soliton_paths_label.setText(f"Soliton paths ({len(self.automata.soliton_paths)}):")
            self.soliton_paths_label.show()
            self.paths.show()
            self.show_matrices.show()
            self.show_end_result.show()
            self.show_animation.show()
            for soliton_path in self.automata.soliton_paths:
                self.paths.addItem(str(soliton_path.path_for_user))
    

    def show_matrices_clicked(self):
        """Method that is called when user clicks button to have the adjacency matrices of every timestep displayed.
        Makes a small window pop up that shows the labelled adjacency matrices and provides a save button.
        """

        def save_matrices():
            """Opens a file dialog in which user can specify a path where a text file with all adjacency matrices should be saved.
            Text file also contains the input string representing the molecule and the soliton path.
            Only allows `.txt` file suffix.
            """
            option = QtWidgets.QFileDialog.Options()
            name = QtWidgets.QFileDialog.getSaveFileName(self.centralwidget, 'Save File', 'matrices.txt', 'Text files (*.txt)', options = option)
            if name != ('', ''): # only do the following if user clicked on save button (without this line the application closes with an error if save action is cancelled)
                path = name[0]
                file = open(path, "w")
                file.write(txt_text)
                file.close()

        index = self.paths.currentIndex()
        desired_path = self.automata.soliton_paths[index]
    
        dlg = QDialog()
        scrollArea = QScrollArea(dlg)
        widget = QtWidgets.QWidget()
        vbox = QtWidgets.QVBoxLayout()
        txt_text = f"Molecule: {self.smiles_string} \n"
        txt_text = txt_text + f"Soliton path: {desired_path.path_for_user} \n \n"
        # labelling of matrix depends on wether we have long node labels ("aa", "ab", ...) or short ones ("a", "b", ...)
        if (len(self.my_graph.labels) - len(self.my_graph.exterior_nodes)) > 26:
            matrix_label_horizontal = "    "
            long_labels = True
        else:
            matrix_label_horizontal = "   "
            long_labels = False
        for key in self.my_graph.labels:
            if long_labels:
                if self.my_graph.labels[key] in self.my_graph.exterior_nodes_reverse and int(self.my_graph.labels[key]) < 10:
                    matrix_label_horizontal = matrix_label_horizontal + f"{self.my_graph.labels[key]}  "
                else:
                    matrix_label_horizontal = matrix_label_horizontal + f"{self.my_graph.labels[key]} "
            else:
                if self.my_graph.labels[key] in self.my_graph.exterior_nodes_reverse and int(self.my_graph.labels[key]) >= 10:
                    matrix_label_horizontal = matrix_label_horizontal + f"{self.my_graph.labels[key]} "
                else:
                    matrix_label_horizontal = matrix_label_horizontal + f"{self.my_graph.labels[key]}  "
        # for show-matrices-window the node labels are added to the matrix string and this labelled matrix is added to scroll area
        # for matrices.txt we add everything (horizontal label and every row of matrix) line by line
        for i in range(len(desired_path.adjacency_matrices_list)):
            txt_text = txt_text + f"Timestep {i}: \n"
            txt_text = txt_text + f"{matrix_label_horizontal} \n"
            matrix_labelled = ""
            matrix_labelled = matrix_labelled + f"{matrix_label_horizontal}\n"
            matrix = str(desired_path.adjacency_matrices_list[i])
            matrix = re.sub(r"[matrix(]", "", matrix)
            matrix = re.sub(r"[)]", "", matrix)
            for j in range(len(matrix.splitlines())):
                if long_labels:
                    if self.my_graph.labels[j] in self.my_graph.exterior_nodes_reverse and int(self.my_graph.labels[j]) < 10:
                        matrix_labelled = matrix_labelled + f"{self.my_graph.labels[j]} {matrix.splitlines()[j]} \n"
                        txt_text = txt_text + f"{self.my_graph.labels[j]} {matrix.splitlines()[j]} \n"
                    else:
                        matrix_labelled = matrix_labelled + f"{self.my_graph.labels[j]}{matrix.splitlines()[j]} \n"
                        txt_text = txt_text + f"{self.my_graph.labels[j]}{matrix.splitlines()[j]} \n"
                else:
                    if self.my_graph.labels[j] in self.my_graph.exterior_nodes_reverse and int(self.my_graph.labels[j]) >= 10:
                        matrix_labelled = matrix_labelled + f"{self.my_graph.labels[j]}{matrix.splitlines()[j][1:]} \n"
                        txt_text = txt_text + f"{self.my_graph.labels[j]}{matrix.splitlines()[j][1:]} \n"
                    else:
                        matrix_labelled = matrix_labelled + f"{self.my_graph.labels[j]}{matrix.splitlines()[j]} \n"
                        txt_text = txt_text + f"{self.my_graph.labels[j]}{matrix.splitlines()[j]} \n"
            txt_text = txt_text + f"\n"
            vbox.addWidget(QtWidgets.QLabel(matrix_labelled))
        widget.setLayout(vbox)
        scrollArea.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        scrollArea.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        scrollArea.setWidgetResizable(True)
        scrollArea.setWidget(widget)
        scrollArea.setStyleSheet("font: 13pt Courier;")
        scrollArea.setFixedSize(540, 405)
        save_button = QtWidgets.QPushButton("Save", dlg)
        save_button.setGeometry(QtCore.QRect(454, 359, 70, 30))
        save_button.clicked.connect(save_matrices)

        dlg.setWindowTitle("Adjacency Matrices")
        dlg.setFixedSize(545, 410)
        dlg.exec_()
    

    def show_end_result_clicked(self):
        """Method that is called when user clicks button to have the resulting graph (after soliton path is traversed) displayed.
        Makes a small window pop up that shows the graph visualisation and provides a save button.
        """

        def save_end_result():
            """Opens a file dialog in which user can specify a path where the resulting graph should be saved.
            Only allows `.jpg`, `.png` and `.jpeg` file suffixes.
            """
            option = QtWidgets.QFileDialog.Options()
            name = QtWidgets.QFileDialog.getSaveFileName(self.centralwidget, 'Save File', 'result.jpg', 'Images (*.jpg *.png *.jpeg)', options = option)
            if name != ('', ''):
                path = name[0]
                imgByteArr = io.BytesIO()
                result_pic.save(imgByteArr, format='PNG')
                imgByteArr = imgByteArr.getvalue()
                file = open(path, "wb")
                file.write(imgByteArr)
                file.close()

        index = self.paths.currentIndex()
        desired_path = self.automata.soliton_paths[index]
        bindings_index = len(desired_path.path) - 1
        result_pic = Visualisation.visualize_soliton_graph(self.my_graph, desired_path.bindings_list[bindings_index], False, True)
        qim = ImageQt(result_pic)

        dlg = QDialog()
        label = QtWidgets.QLabel(dlg)
        label.setGeometry(QtCore.QRect(0, 0, 540, 405))
        label.setPixmap(QtGui.QPixmap.fromImage(qim))
        label.setScaledContents(True)
        save_button = QtWidgets.QPushButton("Save", dlg)
        save_button.setGeometry(QtCore.QRect(470, 375, 70, 30))
        save_button.clicked.connect(save_end_result)

        dlg.setWindowTitle("End result")
        dlg.setFixedSize(545, 410)
        dlg.exec_()


    def show_animation_clicked(self):
        """Method that is called when user clicks button to have the animation of the soliton traversing the graph displayed.
        Makes a small window pop up that shows the animation and provides a save button.
        Instead of displaying the `gif` it uses a sequence of `PIL` images and always shows the next image after a certain time.
        """
        
        def save_animation():
            """Opens a file dialog in which user can specify a path where the animation should be saved.
            Only allows `.gif` file suffix.
            """
            option = QtWidgets.QFileDialog.Options()
            name = QtWidgets.QFileDialog.getSaveFileName(self.centralwidget, 'Save File', 'animation.gif', 'Images (*.gif)', options = option)
            if name != ('', ''):
                path = name[0]
                ani = self.my_animation.graph_animation()
                ani.save(path, writer='pillow', dpi = 600)

        dlg = QDialog()
        label = QtWidgets.QLabel(dlg)
        label.setGeometry(QtCore.QRect(0, 0, 540, 405))
        
        if self.path_index != self.paths.currentIndex():
            self.path_index = self.paths.currentIndex()
            self.desired_path = self.automata.soliton_paths[self.path_index]
            self.my_animation = Animation(self.my_graph, self.desired_path)
            self.pil_images = self.my_animation.pil_images

        self.step = 0
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(lambda: self.update_image(label, self.pil_images, self.desired_path))
        self.timer.start(800) # triggers event every 800 milliseconds
        self.update_image(label, self.pil_images, self.desired_path)
        save_button = QtWidgets.QPushButton("Save", dlg)
        save_button.setGeometry(QtCore.QRect(470, 375, 70, 30))
        save_button.clicked.connect(save_animation)

        dlg.setWindowTitle("Animation")
        dlg.setFixedSize(545, 410)
        dlg.closeEvent = self.stop_animation # stop timer when window is closed
        dlg.exec_()
    

    def update_image(self, label: QtWidgets.QLabel, pil_images: list, desired_path: SolitonPath):
        """Displays next image of animation after a certain time (method is triggered by a timer).

        Args:
            label (QtWidgets.QLabel): Displays the animation's images.
            pil_images (list): `PIL` images that should be displayed step by step.
            desired_path (SolitonPath): Soliton path that is traversed in animation.
        """
        im = pil_images[self.step]
        qim = ImageQt(im)
        label.setPixmap(QtGui.QPixmap.fromImage(qim))
        label.setScaledContents(True)
        self.step += 1
        if self.step == len(desired_path.path): # start animation all over again as soon as end is reached (endless loop)
            self.step = 0


    def stop_animation(self, event):
        """Stops the timer that is used for the animation.
        Without this method the application would crash if the animation window would get closed.
        Args:
            event: Close event of the window that shows animation.
        """
        self.timer.disconnect()


    def hide_retain_space(self, widget):
        """Retains the space of a widget even when it's hidden.

        Args:
            widget (QtWidgets.QWidget ): Widget whose space should be retained.
        """
        retain = widget.sizePolicy()
        retain.setRetainSizeWhenHidden(True)
        widget.setSizePolicy(retain)


    def heightForWidth(self, width: float):
        """Computes height for a given width.

        Args:
            width (float): Given width.

        Returns:
            float: Computed height.
        """
        return math.ceil(width * self.height_for_width_factor)


    def widthForHeight(self, height):
        """Computes width for a given height.

        Args:
            height (float): Given height.

        Returns:
            float: Computed width.
        """
        return math.ceil(height * self.width_for_height_factor)


    def resizeEvent(self, event: QtGui.QResizeEvent):
        """Keeps the right aspect ratio of the graph visualisation image (or welcoming screen in the beginning) when window is resized.

        Args:
            event (QtGui.QResizeEvent): Resize event.
        """
        super(MainWindow, self).resizeEvent(event)
        size = event.size()
        # window height minus height of all widgets below the label result in height of pixmap
        height = size.height() - self.save.height() - self.submit_molecule.height() - self.submit_exterior_nodes.height() - self.paths.height() - self.show_animation.height()
        if self.widthForHeight(height) > size.width(): # pixmap width should be as most as big as width of whole window
            width = size.width()
        else:
            width = self.widthForHeight(height)
        self.display_molecule.setPixmap(QtGui.QPixmap.fromImage(self.qim).scaled(width, height, QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation))
