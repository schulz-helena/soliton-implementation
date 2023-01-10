"""GUI of the software.
 """
import copy
import io
import math
import os
import re
import time
from threading import Thread

import networkx as nx
import res.resources
from gui.startscreen import Startscreen
from PIL.ImageQt import ImageQt
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QDialog, QMainWindow, QMessageBox, QScrollArea
from soliton_classes.mini_soliton_automata import MiniSolitonAutomata
from soliton_classes.multiwave_soliton_automata import MultiwaveSolitonAutomata
from soliton_classes.soliton_automata import SolitonAutomata
from soliton_classes.soliton_graph import SolitonGraph
from soliton_classes.soliton_path import SolitonPath
from soliton_classes.traversal import Traversal
from visualisations.animation import Animation
from visualisations.visualisation import Visualisation


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
        # central widget is a stacked layout (so user can switch between the two different windows)
        self.central_wid = QtWidgets.QWidget()
        self.layout_for_wids = QtWidgets.QStackedLayout()
        self.central_wid.setStyleSheet("""background: white;""")
        # move window to the top of the screen + to the center horizontally
        qt_rectangle = self.frameGeometry()
        center_point = QtWidgets.QDesktopWidget().availableGeometry().center()
        qt_rectangle.moveCenter(center_point)
        self.move(qt_rectangle.topLeft().x(), 0)
        self.resize(QtCore.QSize(600, 650))


        # SOLITON AUTOMATA WIDGET
        self.wid_single = QtWidgets.QWidget()
        self.gridLayout = QtWidgets.QGridLayout(self.wid_single)
        # Row 0:
        # Rectangle that displays the molecule
        self.display_molecule = QtWidgets.QLabel(self.wid_single)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Ignored)
        self.display_molecule.setSizePolicy(sizePolicy)
        self.display_molecule.setMinimumSize(QtCore.QSize(50, 37.5))
        startscreen = Startscreen().image
        self.qim = ImageQt(startscreen)
        self.display_molecule.setPixmap(QtGui.QPixmap.fromImage(self.qim).scaled(self.display_molecule.size(), QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation))
        self.gridLayout.addWidget(self.display_molecule, 0, 0, 2, 5, alignment = QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
        # Button that changes between the two windows
        self.change_window_button = QtWidgets.QPushButton("Switch Mode")
        self.change_window_button.setMinimumSize(QtCore.QSize(0, 20))
        self.gridLayout.addWidget(self.change_window_button, 0, 0, 1, 1)
        self.change_window_button.setStyleSheet("QPushButton {border-radius: 10px;}")
        # Row 1: -
        # Row 2:
        # "Traversal Mode" Checkbox
        self.traversal_mode = QtWidgets.QCheckBox("Traversal Mode")
        self.traversal_mode.setChecked(False)
        self.gridLayout.addWidget(self.traversal_mode, 2, 0, 1, 1)
        # Groupbox containg last two elements of row 2
        self.row2 = QtWidgets.QGroupBox()
        self.minigrid2 = QtWidgets.QGridLayout(self.row2)
        self.minigrid2.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.addWidget(self.row2, 2, 3, 1, 1)
        # Info button
        self.mol_info = QtWidgets.QPushButton(self.wid_single)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        self.mol_info.setSizePolicy(sizePolicy)
        self.mol_info.setMinimumSize(QtCore.QSize(0, 32))
        self.minigrid2.addWidget(self.mol_info, 2, 0, 1, 1)
        # Save button for the molecule
        self.save = QtWidgets.QPushButton(self.wid_single)
        self.save.setStyleSheet("QPushButton {background-color: rgb(191, 207, 255); image: url(:/icons/save.svg);} QPushButton::pressed {background-color : rgb(132, 145, 193);}")
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        self.save.setSizePolicy(sizePolicy)
        self.save.setMaximumSize(QtCore.QSize(16777215, 32))
        self.save.setMinimumSize(QtCore.QSize(0, 32))
        self.minigrid2.addWidget(self.save, 2, 1, 1, 1)
        # Row 3:
        # "Molecule" label
        self.molecule_label = QtWidgets.QLabel(self.wid_single)
        self.gridLayout.addWidget(self.molecule_label, 3, 0, 1, 1)
        # Text field for molecule
        self.molecule_lineedit = QtWidgets.QLineEdit(self.wid_single)
        self.gridLayout.addWidget(self.molecule_lineedit, 3, 1, 1, 2)
        # Submit button for molecule
        self.submit_molecule = QtWidgets.QPushButton(self.wid_single)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        self.submit_molecule.setSizePolicy(sizePolicy)
        self.submit_molecule.setMinimumSize(QtCore.QSize(0, 32))
        self.gridLayout.addWidget(self.submit_molecule, 3, 3, 1, 1)
        # Row 4: 
        # "Exterior nodes" label
        self.exterior_nodes_label = QtWidgets.QLabel(self.wid_single)
        self.gridLayout.addWidget(self.exterior_nodes_label, 4, 0, 1, 1)
        # Groupbox containg middle elements of row 4
        self.row4 = QtWidgets.QGroupBox()
        self.minigrid4 = QtWidgets.QGridLayout(self.row4)
        self.minigrid4.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.addWidget(self.row4, 4, 1, 1, 2)
        # "All" Checkbox
        self.all_exterior_nodes = QtWidgets.QCheckBox("All")
        self.all_exterior_nodes.setChecked(False)
        self.minigrid4.addWidget(self.all_exterior_nodes, 4, 0, 1, 1)
        # Combobox to choose first exterior node
        self.node_1 = QtWidgets.QComboBox(self.row4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        self.node_1.setSizePolicy(sizePolicy)
        self.minigrid4.addWidget(self.node_1, 4, 1, 1, 1)
        # "&" label
        self.exterior_nodes_label2 = QtWidgets.QLabel(self.row4)
        self.minigrid4.addWidget(self.exterior_nodes_label2, 4, 2, 1, 1)
        # Combobox to choose second exterior node
        self.node_2 = QtWidgets.QComboBox(self.row4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        self.node_2.setSizePolicy(sizePolicy)
        self.minigrid4.addWidget(self.node_2, 4, 3, 1, 1)
        # Submit button for exterior nodes
        self.submit_exterior_nodes = QtWidgets.QPushButton(self.wid_single)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        self.submit_exterior_nodes.setSizePolicy(sizePolicy)
        self.submit_exterior_nodes.setMinimumSize(QtCore.QSize(0, 32))
        self.gridLayout.addWidget(self.submit_exterior_nodes, 4, 3, 1, 1)
        # Row 5:
        # "Soliton paths" label
        self.soliton_paths_label = QtWidgets.QLabel(self.wid_single)
        self.gridLayout.addWidget(self.soliton_paths_label, 5, 0, 1, 1)
        # Combobox to choose a soliton path
        self.paths = QtWidgets.QComboBox(self.wid_single)
        self.gridLayout.addWidget(self.paths, 5, 1, 1, 3)
        # Row 6:
        # "Show matrices" button
        self.show_matrices = QtWidgets.QPushButton(self.wid_single)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        self.show_matrices.setSizePolicy(sizePolicy)
        self.show_matrices.setMinimumSize(QtCore.QSize(0, 32))
        self.gridLayout.addWidget(self.show_matrices, 6, 1, 1, 1)
        # "Show end result" button
        self.show_end_result = QtWidgets.QPushButton(self.wid_single)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        self.show_end_result.setSizePolicy(sizePolicy)
        self.show_end_result.setMinimumSize(QtCore.QSize(0, 32))
        self.gridLayout.addWidget(self.show_end_result, 6, 2, 1, 1)
        # "Show animation" button
        self.show_animation = QtWidgets.QPushButton(self.wid_single)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        self.show_animation.setSizePolicy(sizePolicy)
        self.show_animation.setMinimumSize(QtCore.QSize(0, 32))
        self.gridLayout.addWidget(self.show_animation, 6, 3, 1, 1)
        # Menubar:
        self.menubar = QtWidgets.QMenuBar(self)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 564, 24))
        self.setMenuBar(self.menubar)

        # Function connections of different widgets:
        self.change_window_button.clicked.connect(self.change_window)
        self.traversal_mode.stateChanged.connect(lambda:self.change_mode(self.traversal_mode))
        self.mol_info.clicked.connect(self.mol_info_clicked)
        self.save.clicked.connect(self.save_clicked)
        self.submit_molecule.clicked.connect(self.submit_molecule_clicked)
        self.all_exterior_nodes.stateChanged.connect(self.all_exterior_nodes_statechanged)
        self.submit_exterior_nodes.clicked.connect(self.submit_exterior_nodes_clicked)
        self.show_matrices.clicked.connect(self.show_matrices_clicked)
        self.show_end_result.clicked.connect(self.show_end_result_clicked)
        self.show_animation.clicked.connect(lambda:self.show_animation_clicked(self.show_animation))
        # Hide most widgets at the beginning (while retaining space) 
        self.hide_retain_space([self.traversal_mode, self.row2, self.exterior_nodes_label, self.row4, self.submit_exterior_nodes, self.soliton_paths_label, self.paths, self.show_matrices, self.show_end_result, self.show_animation])
        self.hide_multiple([self.traversal_mode, self.row2, self.exterior_nodes_label, self.row4, self.submit_exterior_nodes, self.soliton_paths_label, self.paths, self.show_matrices, self.show_end_result, self.show_animation])
        # Stylesheet:
        readme_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../styles.css')
        with open(readme_path, "r", encoding="utf-8") as fh:
            self.style = fh.read()
        self.wid_single.setStyleSheet(self.style)


        # MULTI-WAVE SOLITON AUTOMATA WIDGET
        self.wid_mult = QtWidgets.QWidget()
        self.gridLayout_m = QtWidgets.QGridLayout(self.wid_mult)
        # Row 0:
        # Rectangle that displays the molecule
        self.display_molecule_m = QtWidgets.QLabel(self.wid_mult)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Ignored)
        self.display_molecule_m.setSizePolicy(sizePolicy)
        self.display_molecule_m.setMinimumSize(QtCore.QSize(50, 37.5))
        startscreen_m = Startscreen().image
        self.qim_m = ImageQt(startscreen_m)
        self.display_molecule_m.setPixmap(QtGui.QPixmap.fromImage(self.qim_m).scaled(self.display_molecule_m.size(), QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation))
        self.gridLayout_m.addWidget(self.display_molecule_m, 0, 0, 2, 5, alignment = QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
        # Button that changes between the two windows
        self.change_window_button_m = QtWidgets.QPushButton("Switch Mode")
        self.change_window_button_m.setMinimumSize(QtCore.QSize(0, 20))
        self.gridLayout_m.addWidget(self.change_window_button_m, 0, 0, 1, 1)
        self.change_window_button_m.setStyleSheet("QPushButton {border-radius: 10px;}")
        # Row 1: -
        # Row 2:
        # "Traversal Mode" Checkbox
        self.traversal_mode_m = QtWidgets.QCheckBox("Traversal Mode")
        self.traversal_mode_m.setChecked(False)
        self.gridLayout_m.addWidget(self.traversal_mode_m, 2, 0, 1, 1)
        # Groupbox containg last two elements of row 2
        self.row2_m = QtWidgets.QGroupBox()
        self.minigrid2_m = QtWidgets.QGridLayout(self.row2_m)
        self.minigrid2_m.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_m.addWidget(self.row2_m, 2, 3, 1, 1)
        # Info button
        self.mol_info_m = QtWidgets.QPushButton(self.row2_m)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        self.mol_info_m.setSizePolicy(sizePolicy)
        self.mol_info_m.setMinimumSize(QtCore.QSize(0, 32))
        self.minigrid2_m.addWidget(self.mol_info_m, 2, 0, 1, 1)
        # Save button for the molecule
        self.save_m = QtWidgets.QPushButton(self.row2_m)
        self.save_m.setStyleSheet("QPushButton {background-color: rgb(149, 221, 185); image: url(:/icons/save.svg);} QPushButton::pressed {background-color : rgb(90, 159, 123);}")
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        self.save_m.setSizePolicy(sizePolicy)
        self.save_m.setMaximumSize(QtCore.QSize(16777215, 32))
        self.save_m.setMinimumSize(QtCore.QSize(0, 32))
        self.minigrid2_m.addWidget(self.save_m, 2, 1, 1, 1)
        # Row 3:
        # "Molecule" label
        self.molecule_label_m = QtWidgets.QLabel(self.wid_mult)
        self.gridLayout_m.addWidget(self.molecule_label_m, 3, 0, 1, 1)
        # Text field for molecule
        self.molecule_lineedit_m = QtWidgets.QLineEdit(self.wid_mult)
        self.gridLayout_m.addWidget(self.molecule_lineedit_m, 3, 1, 1, 2)
        # Submit button for molecule
        self.submit_molecule_m = QtWidgets.QPushButton(self.wid_mult)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        self.submit_molecule_m.setSizePolicy(sizePolicy)
        self.submit_molecule_m.setMinimumSize(QtCore.QSize(0, 32))
        self.gridLayout_m.addWidget(self.submit_molecule_m, 3, 3, 1, 1)
        # Row 4: 
        # "Set of bursts" label
        self.set_of_bursts_label = QtWidgets.QLabel(self.wid_mult)
        self.gridLayout_m.addWidget(self.set_of_bursts_label, 4, 0, 1, 1)
        # Text field for set of bursts
        self.set_of_bursts_lineedit = QtWidgets.QLineEdit(self.wid_mult)
        self.gridLayout_m.addWidget(self.set_of_bursts_lineedit, 4, 1, 1, 2)
        # Submit button for set of bursts
        self.submit_set_of_bursts = QtWidgets.QPushButton(self.wid_mult)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        self.submit_set_of_bursts.setSizePolicy(sizePolicy)
        self.submit_set_of_bursts.setMinimumSize(QtCore.QSize(0, 32))
        self.gridLayout_m.addWidget(self.submit_set_of_bursts, 4, 3, 1, 1)
        # Row 5:
        # "Bursts" label
        self.bursts_label = QtWidgets.QLabel(self.wid_mult)
        self.gridLayout_m.addWidget(self.bursts_label, 5, 0, 1, 1)
        # Groupbox containg middle elements of row 5
        self.row5_m = QtWidgets.QGroupBox()
        self.minigrid4_m = QtWidgets.QGridLayout(self.row5_m)
        self.minigrid4_m.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_m.addWidget(self.row5_m, 5, 1, 1, 2)
        # "All" Checkbox
        self.all_bursts = QtWidgets.QCheckBox("All")
        self.all_bursts.setChecked(False)
        self.minigrid4_m.addWidget(self.all_bursts, 5, 0, 1, 1)
        # Combobox to choose burst
        self.burst = QtWidgets.QComboBox(self.row5_m)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        self.burst.setSizePolicy(sizePolicy)
        self.minigrid4_m.addWidget(self.burst, 5, 1, 1, 1)
        # Submit button for burst
        self.submit_burst = QtWidgets.QPushButton(self.wid_mult)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        self.submit_burst.setSizePolicy(sizePolicy)
        self.submit_burst.setMinimumSize(QtCore.QSize(0, 32))
        self.gridLayout_m.addWidget(self.submit_burst, 5, 3, 1, 1)
        # Row 6:
        # "Traversals" label
        self.traversals_label = QtWidgets.QLabel(self.wid_mult)
        self.gridLayout_m.addWidget(self.traversals_label, 6, 0, 1, 1)
        # Combobox to choose a traversal
        self.traversals = QtWidgets.QComboBox(self.wid_mult)
        self.gridLayout_m.addWidget(self.traversals, 6, 1, 1, 3)
        # Row 7:
        # "Show matrices" button
        self.show_matrices_m = QtWidgets.QPushButton(self.wid_mult)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        self.show_matrices_m.setSizePolicy(sizePolicy)
        self.show_matrices_m.setMinimumSize(QtCore.QSize(0, 32))
        self.gridLayout_m.addWidget(self.show_matrices_m, 7, 1, 1, 1)
        # "Show end result" button
        self.show_end_result_m = QtWidgets.QPushButton(self.wid_mult)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        self.show_end_result_m.setSizePolicy(sizePolicy)
        self.show_end_result_m.setMinimumSize(QtCore.QSize(0, 32))
        self.gridLayout_m.addWidget(self.show_end_result_m, 7, 2, 1, 1)
        # "Show animation" button
        self.show_animation_m = QtWidgets.QPushButton(self.wid_mult)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        self.show_animation_m.setSizePolicy(sizePolicy)
        self.show_animation_m.setMinimumSize(QtCore.QSize(0, 32))
        self.gridLayout_m.addWidget(self.show_animation_m, 7, 3, 1, 1)

        # Function connections of different widgets:
        self.change_window_button_m.clicked.connect(self.change_window)
        self.traversal_mode_m.stateChanged.connect(lambda:self.change_mode_m(self.traversal_mode_m))
        self.mol_info_m.clicked.connect(self.mol_info_clicked_m)
        self.save_m.clicked.connect(self.save_clicked_m)
        self.all_bursts.stateChanged.connect(self.all_bursts_statechanged)
        self.submit_molecule_m.clicked.connect(self.submit_molecule_clicked_m)
        self.submit_set_of_bursts.clicked.connect(self.submit_set_of_bursts_clicked)
        self.submit_burst.clicked.connect(self.submit_burst_clicked)
        self.traversals.currentIndexChanged.connect(lambda: self.endless_loop_picked(self.traversals))
        self.show_matrices_m.clicked.connect(self.show_matrices_clicked_m)
        self.show_end_result_m.clicked.connect(self.show_end_result_clicked_m)
        self.show_animation_m.clicked.connect(lambda:self.show_animation_clicked(self.show_animation_m))
        # Hide most widgets at the beginning (while retaining space)
        self.hide_retain_space([self.traversal_mode_m, self.save_m, self.mol_info_m, self.set_of_bursts_label, self.set_of_bursts_lineedit, self.submit_set_of_bursts, self.bursts_label, self.row5_m, self.submit_burst, self.traversals_label, self.traversals, self.show_matrices_m, self.show_end_result_m, self.show_animation_m])
        self.hide_multiple([self.traversal_mode_m, self.save_m, self.mol_info_m, self.set_of_bursts_label, self.set_of_bursts_lineedit, self.submit_set_of_bursts, self.bursts_label, self.row5_m, self.submit_burst, self.traversals_label, self.traversals, self.show_matrices_m, self.show_end_result_m, self.show_animation_m])
        # Stylesheet:
        readme_path_m = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../styles_m.css')
        with open(readme_path_m, "r", encoding="utf-8") as fh_m:
            self.style_m = fh_m.read()
        self.wid_mult.setStyleSheet(self.style_m)


        # Add both widgets to stacked layout
        self.layout_for_wids.addWidget(self.wid_single)
        self.layout_for_wids.addWidget(self.wid_mult)
        self.central_wid.setLayout(self.layout_for_wids)
        self.setCentralWidget(self.central_wid)
        self.front_wid = 1
        self.retranslateUi()
        QtCore.QMetaObject.connectSlotsByName(self)


    def retranslateUi(self):
        """Implements multi-language suppport. Is generated automatically when using PyQt5 UI code generator.
        """
        _translate = QtCore.QCoreApplication.translate
        self.setWindowTitle(_translate("MainWindow", "Soliton Automata Software"))
        self.mol_info.setText(_translate("MainWindow", "Info"))
        self.molecule_label.setText(_translate("MainWindow", "Molecule:"))
        self.submit_molecule.setText(_translate("MainWindow", "Submit"))
        self.exterior_nodes_label.setText(_translate("MainWindow", "Exterior nodes:"))
        self.exterior_nodes_label2.setText(_translate("MainWindow", "&"))
        self.submit_exterior_nodes.setText(_translate("MainWindow", "Submit"))
        self.soliton_paths_label.setText(_translate("MainWindow", "Soliton paths:"))
        self.show_matrices.setText(_translate("MainWindow", "Show matrices"))
        self.show_end_result.setText(_translate("MainWindow", "Show end result"))
        self.show_animation.setText(_translate("MainWindow", "Show animation"))

        self.mol_info_m.setText(_translate("MainWindow", "Info"))
        self.molecule_label_m.setText(_translate("MainWindow", "Molecule:"))
        self.submit_molecule_m.setText(_translate("MainWindow", "Submit"))
        self.set_of_bursts_label.setText(_translate("MainWindow", "Set of bursts:"))
        self.submit_set_of_bursts.setText(_translate("MainWindow", "Submit"))
        self.bursts_label.setText(_translate("MainWindow", "Bursts:"))
        self.submit_burst.setText(_translate("MainWindow", "Submit"))
        self.traversals_label.setText(_translate("MainWindow", "Sets of paths:"))
        self.show_matrices_m.setText(_translate("MainWindow", "Show matrices"))
        self.show_end_result_m.setText(_translate("MainWindow", "Show end result"))
        self.show_animation_m.setText(_translate("MainWindow", "Show animation"))


    def change_window(self):
        if self.front_wid == 1:
            self.layout_for_wids.setCurrentIndex(1)
            self.front_wid = 2
        else:
            self.layout_for_wids.setCurrentIndex(0)
            self.front_wid = 1

    
    def mol_info_clicked(self):
        if self.automata is None:
            self.automata = SolitonAutomata(self.my_graph)
        
        dlg = QDialog(self.wid_single)
        grid = QtWidgets.QGridLayout(dlg)
        label_det = QtWidgets.QLabel(dlg)
        label_det.setText("Deterministic:")
        grid.addWidget(label_det, 0, 0, 1, 1)
        det_bool = QtWidgets.QLabel(dlg)
        if self.automata.deterministic:
            det_bool.setText("Yes")
        else: det_bool.setText("No")
        grid.addWidget(det_bool, 0, 1, 1, 1)
        label_strong_det = QtWidgets.QLabel(dlg)
        label_strong_det.setText("Strongly deterministic:")
        grid.addWidget(label_strong_det, 1, 0, 1, 1)
        strong_det_bool = QtWidgets.QLabel(dlg)
        if self.automata.strongly_deterministic:
            strong_det_bool.setText("Yes")
        else: strong_det_bool.setText("No")
        grid.addWidget(strong_det_bool, 1, 1, 1, 1)
        label_imp_paths = QtWidgets.QLabel(dlg)
        label_imp_paths.setText("Impervious path(s):")
        grid.addWidget(label_imp_paths, 2, 0, 1, 1, alignment = QtCore.Qt.AlignTop)

        group = QtWidgets.QGroupBox()
        layout = QtWidgets.QGridLayout(group)
        layout.setContentsMargins(0, 0, 0, 0)
        impervs = self.automata.find_impervious_paths()
        if impervs == []:
            layout.addWidget(QtWidgets.QLabel("-"))
        for i, path in enumerate(impervs):
            layout.addWidget(QtWidgets.QLabel(path), i, 0, 1, 1, alignment = QtCore.Qt.AlignTop)
        grid.addWidget(group, 2, 1, 1, 1, alignment = QtCore.Qt.AlignTop)

        dlg.setWindowTitle("Info")
        dlg.exec_()

    
    def mol_info_clicked_m(self):
        
        dlg = QDialog(self.wid_single)
        grid = QtWidgets.QGridLayout(dlg)
        label_det = QtWidgets.QLabel(dlg)
        label_det.setText("Deterministic:")
        grid.addWidget(label_det, 0, 0, 1, 1)
        det_bool = QtWidgets.QLabel(dlg)
        if self.multi_automata.deterministic:
            det_bool.setText("Yes")
        else: det_bool.setText("No")
        grid.addWidget(det_bool, 0, 1, 1, 1)
        label_strong_det = QtWidgets.QLabel(dlg)
        label_strong_det.setText("Strongly deterministic:")
        grid.addWidget(label_strong_det, 1, 0, 1, 1)
        strong_det_bool = QtWidgets.QLabel(dlg)
        if self.multi_automata.strongly_deterministic:
            strong_det_bool.setText("Yes")
        else: strong_det_bool.setText("No")
        grid.addWidget(strong_det_bool, 1, 1, 1, 1)

        dlg.setWindowTitle("Info")
        dlg.exec_()

    
    def endless_loop_picked(self, combobox: QtWidgets.QComboBox):
        if not isinstance(self.found_traversals[combobox.currentIndex()], Traversal):
            try:
                self.show_matrices_m.clicked.disconnect()
                self.show_matrices_m.setStyleSheet("QPushButton {background-color: rgb(230, 230, 230);}")
                self.show_end_result_m.clicked.disconnect()
                self.show_end_result_m.setStyleSheet("QPushButton {background-color: rgb(230, 230, 230);}")
            except:
                pass
        else:
            try:
                self.show_matrices_m.clicked.disconnect()
                self.show_matrices_m.clicked.connect(self.show_matrices_clicked_m)
                self.show_matrices_m.setStyleSheet("QPushButton {background-color: rgb(149, 221, 185);} QPushButton::pressed {background-color : rgb(90, 159, 123);}")
                self.show_end_result_m.clicked.disconnect()
                self.show_end_result_m.clicked.connect(self.show_end_result_clicked_m)
                self.show_end_result_m.setStyleSheet("QPushButton {background-color: rgb(149, 221, 185);} QPushButton::pressed {background-color : rgb(90, 159, 123);}")
            except:
                self.show_matrices_m.clicked.connect(self.show_matrices_clicked_m)
                self.show_matrices_m.setStyleSheet("QPushButton {background-color: rgb(149, 221, 185);} QPushButton::pressed {background-color : rgb(90, 159, 123);}")
                self.show_end_result_m.clicked.connect(self.show_end_result_clicked_m)
                self.show_end_result_m.setStyleSheet("QPushButton {background-color: rgb(149, 221, 185);} QPushButton::pressed {background-color : rgb(90, 159, 123);}")


    def change_mode(self, checkbox: QtWidgets.QCheckBox):
        """Realizes the change between being in traversal mode and not being in traversal mode.
        In traversal mode, a soliton automata can be traversed by using an end result as the new soliton graph.
        When in traversal mode, the input molecule can't be edited.

        Args:
            checkbox (QtWidgets.QCheckBox): The checkbox that changes the mode.
        """
        if checkbox.isChecked() == True:
            # make text field uneditable and button unclickable, turn both objects grey to make those properties visually recognizable
            self.submit_molecule.clicked.disconnect()
            self.submit_molecule.setStyleSheet("QPushButton {background-color: rgb(230, 230, 230);}")
            self.molecule_lineedit.setReadOnly(True)
            self.molecule_lineedit.setStyleSheet("QLineEdit {border: 2px solid rgb(230, 230, 230);border-radius: 10px;padding: 0 8px;}")
            self.molecule_label.setText("Original molecule:")
        else:
            self.submit_molecule.clicked.connect(self.submit_molecule_clicked)
            self.submit_molecule.setStyleSheet("QPushButton {background-color: rgb(191, 207, 255);} QPushButton::pressed {background-color : rgb(132, 145, 193);}")
            self.molecule_lineedit.setReadOnly(False)
            self.molecule_lineedit.setStyleSheet("QLineEdit {border: 2px solid rgb(191, 207, 255);border-radius: 10px;padding: 0 8px;}")
            self.molecule_label.setText("Molecule:")
            self.submit_molecule_clicked()

    
    def change_mode_m(self, checkbox: QtWidgets.QCheckBox):
        """Realizes the change between being in traversal mode and not being in traversal mode.
        In traversal mode, a soliton automata can be traversed by using an end result as the new soliton graph.
        When in traversal mode, the input molecule can't be edited.

        Args:
            checkbox (QtWidgets.QCheckBox): The checkbox that changes the mode.
        """
        if checkbox.isChecked() == True:
            # make text field uneditable and button unclickable, turn both objects grey to make those properties visually recognizable
            self.submit_molecule_m.clicked.disconnect()
            self.submit_molecule_m.setStyleSheet("QPushButton {background-color: rgb(230, 230, 230);}")
            self.molecule_lineedit_m.setReadOnly(True)
            self.molecule_lineedit_m.setStyleSheet("QLineEdit {border: 2px solid rgb(230, 230, 230);border-radius: 10px;padding: 0 8px;}")
            self.molecule_label_m.setText("Original molecule:")
            self.submit_set_of_bursts.clicked.disconnect()
            self.submit_set_of_bursts.setStyleSheet("QPushButton {background-color: rgb(230, 230, 230);}")
            self.set_of_bursts_lineedit.setReadOnly(True)
            self.set_of_bursts_lineedit.setStyleSheet("QLineEdit {border: 2px solid rgb(230, 230, 230);border-radius: 10px;padding: 0 8px;}")
        else:
            self.submit_molecule_m.clicked.connect(self.submit_molecule_clicked_m)
            self.submit_molecule_m.setStyleSheet("QPushButton {background-color: rgb(149, 221, 185);} QPushButton::pressed {background-color : rgb(90, 159, 123);}")
            self.molecule_lineedit_m.setReadOnly(False)
            self.molecule_lineedit_m.setStyleSheet("QLineEdit {border: 2px solid rgb(149, 221, 185);border-radius: 10px;padding: 0 8px;}")
            self.molecule_label_m.setText("Molecule:")
            self.submit_molecule_clicked_m()
            self.submit_set_of_bursts.clicked.connect(self.submit_set_of_bursts_clicked)
            self.submit_set_of_bursts.setStyleSheet("QPushButton {background-color: rgb(149, 221, 185);} QPushButton::pressed {background-color : rgb(90, 159, 123);}")
            self.set_of_bursts_lineedit.setReadOnly(False)
            self.set_of_bursts_lineedit.setStyleSheet("QLineEdit {border: 2px solid rgb(149, 221, 185);border-radius: 10px;padding: 0 8px;}")


    def submit_molecule_clicked(self):
        """Method that is called when user clicks button to submit the specified molecule.
        Catches errors if user used the wrong syntax or specified a molecule that does not fulfill the requirements of a soliton graph.
        If the user's molecule is valid it displays the graph of the molecule. It then also reveals a save button for the graph visualisation and
        all the necessary widgets for the user to choose a pair of exterior nodes.
        """
        self.node_1.clear()
        self.node_2.clear()
        smiles_string = self.molecule_lineedit.text()
        self.automata = None
        try:
            self.my_graph = SolitonGraph(smiles_string)
            errors = self.my_graph.validate_soliton_graph()
            self.graph_pic = Visualisation.visualize_soliton_graph(self.my_graph, self.my_graph.bindings, False, True)
            self.qim = ImageQt(self.graph_pic)
            self.display_molecule.setPixmap(QtGui.QPixmap.fromImage(self.qim).scaled(self.display_molecule.size(), QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation))
            if errors != []:
                self.hide_multiple([self.traversal_mode, self.row2, self.exterior_nodes_label, self.row4, self.submit_exterior_nodes, self.soliton_paths_label, self.paths, self.show_matrices, self.show_end_result, self.show_animation])
                msg = QMessageBox(self.wid_single)
                msg.setStyleSheet(" QPushButton{ height: 32px; width: 130px;}")
                msg.setWindowTitle("No soliton graph")
                msg.setText("You specified a molecule that does not fulfill the requirements of a soliton graph.")
                msg.setIcon(QMessageBox.Warning)
                msg.setStandardButtons(QMessageBox.Retry)
                msg.setInformativeText("See details for all incorrect parts of your molecule.")
                details = ""
                for i, error in enumerate(errors):
                    if i == len(errors) - 1:
                        details = details + f"- {error}"
                    else:
                        details = details + f"- {error}" + "\n"
                msg.setDetailedText(details)
                x = msg.exec_() # show messagebox
            elif self.my_graph.exterior_nodes_name_collision() == True:
                self.hide_multiple([self.traversal_mode, self.row2, self.exterior_nodes_label, self.row4, self.submit_exterior_nodes, self.soliton_paths_label, self.paths, self.show_matrices, self.show_end_result, self.show_animation])
                msg = QMessageBox(self.wid_single)
                msg.setStyleSheet(" QPushButton{ height: 32px; width: 130px}")
                msg.setWindowTitle("Name collision")
                msg.setText("You specified two or more exterior nodes with the same name.")
                msg.setIcon(QMessageBox.Warning)
                msg.setStandardButtons(QMessageBox.Retry)
                x = msg.exec_()
            else:
                self.show_multiple([self.traversal_mode, self.row2, self.exterior_nodes_label, self.row4, self.submit_exterior_nodes])
                self.hide_multiple([self.soliton_paths_label, self.paths, self.show_matrices, self.show_end_result, self.show_animation])
                for key in self.my_graph.exterior_nodes_reverse:
                    self.node_1.addItem(key)
                    self.node_2.addItem(key)
        except:
            self.hide_multiple([self.traversal_mode, self.row2, self.exterior_nodes_label, self.row4, self.submit_exterior_nodes, self.soliton_paths_label, self.paths, self.show_matrices, self.show_end_result, self.show_animation])
            msg = QMessageBox(self.wid_single)
            msg.setStyleSheet(" QPushButton{ height: 32px; width: 130px}")
            msg.setWindowTitle("Incorrect input")
            msg.setText("The syntax of your input is not correct.")
            msg.setIcon(QMessageBox.Warning)
            msg.setStandardButtons(QMessageBox.Retry)
            msg.setInformativeText("Please try again with another input string.")
            details = f"Reminder - this is how you define a molecule: \n"
            details = details + f"- Carbon atoms are marked with 'C' \n"
            details = details + f"- Single bonds are marked with '-' or no character at all \n"
            details = details + f"- Double bonds are marked with '=' \n"
            details = details + f"- Branches are embedded in round brackets (e.g. 'C(=CC=C)C')\n"
            details = details + f"- The two connecting atoms of a ring are marked with the same number (e.g. 'C1' and 'C1') \n"
            details = details + "- Exterior nodes are marked with braces and a number (e.g. '{=1}')"
            msg.setDetailedText(details)
            x = msg.exec_()


    def submit_molecule_clicked_m(self):
    
        smiles_string = self.molecule_lineedit_m.text()
        try:
            self.my_graph_m = SolitonGraph(smiles_string)
            errors = self.my_graph_m.validate_soliton_graph()
            self.graph_pic_m = Visualisation.visualize_soliton_graph(self.my_graph_m, self.my_graph_m.bindings, False, True)
            self.qim_m = ImageQt(self.graph_pic_m)
            self.display_molecule_m.setPixmap(QtGui.QPixmap.fromImage(self.qim_m).scaled(self.display_molecule_m.size(), QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation))
            self.set_of_bursts_lineedit.clear()
            if errors != []:
                self.hide_multiple([self.traversal_mode_m, self.save_m, self.mol_info_m, self.set_of_bursts_label, self.set_of_bursts_lineedit, self.submit_set_of_bursts, self.bursts_label, self.row5_m, self.submit_burst, self.traversals_label, self.traversals, self.show_matrices_m, self.show_end_result_m, self.show_animation_m])
                msg = QMessageBox(self.wid_mult)
                msg.setStyleSheet(" QPushButton{ height: 32px; width: 130px}")
                msg.setWindowTitle("No soliton graph")
                msg.setText("You specified a molecule that does not fulfill the requirements of a soliton graph.")
                msg.setIcon(QMessageBox.Warning)
                msg.setStandardButtons(QMessageBox.Retry)
                msg.setInformativeText("See details for all incorrect parts of your molecule.")
                details = ""
                for i, error in enumerate(errors):
                    if i == len(errors) - 1:
                        details = details + f"- {error}"
                    else:
                        details = details + f"- {error}" + "\n"
                msg.setDetailedText(details)
                x = msg.exec_() # show messagebox
            elif self.my_graph_m.exterior_nodes_name_collision() == True:
                self.hide_multiple([self.traversal_mode_m, self.save_m, self.mol_info_m, self.set_of_bursts_label, self.set_of_bursts_lineedit, self.submit_set_of_bursts, self.bursts_label, self.row5_m, self.submit_burst, self.traversals_label, self.traversals, self.show_matrices_m, self.show_end_result_m, self.show_animation_m])
                msg = QMessageBox(self.wid_mult)
                msg.setStyleSheet(" QPushButton{ height: 32px; width: 130px}")
                msg.setWindowTitle("Name collision")
                msg.setText("You specified two or more exterior nodes with the same name.")
                msg.setIcon(QMessageBox.Warning)
                msg.setStandardButtons(QMessageBox.Retry)
                x = msg.exec_()
            else:
                self.show_multiple([self.traversal_mode_m, self.save_m, self.set_of_bursts_label, self.set_of_bursts_lineedit, self.submit_set_of_bursts])
                self.hide_multiple([self.mol_info_m, self.bursts_label, self.row5_m, self.submit_burst, self.traversals_label, self.traversals, self.show_matrices_m, self.show_end_result_m, self.show_animation_m])
        except:
            self.hide_multiple([self.traversal_mode_m, self.save_m, self.mol_info_m, self.set_of_bursts_label, self.set_of_bursts_lineedit, self.submit_set_of_bursts, self.mol_info_m, self.bursts_label, self.row5_m, self.submit_burst, self.traversals_label, self.traversals, self.show_matrices_m, self.show_end_result_m, self.show_animation_m])
            msg = QMessageBox(self.wid_mult)
            msg.setStyleSheet(" QPushButton{ height: 32px; width: 130px}")
            msg.setWindowTitle("Incorrect input")
            msg.setText("The syntax of your input is not correct.")
            msg.setIcon(QMessageBox.Warning)
            msg.setStandardButtons(QMessageBox.Retry)
            msg.setInformativeText("Please try again with another input string.")
            details = f"Reminder - this is how you define a molecule: \n"
            details = details + f"- Carbon atoms are marked with 'C' \n"
            details = details + f"- Single bonds are marked with '-' or no character at all \n"
            details = details + f"- Double bonds are marked with '=' \n"
            details = details + f"- Branches are embedded in round brackets (e.g. 'C(=CC=C)C')\n"
            details = details + f"- The two connecting atoms of a ring are marked with the same number (e.g. 'C1' and 'C1') \n"
            details = details + "- Exterior nodes are marked with braces and a number (e.g. '{=1}')"
            msg.setDetailedText(details)
            x = msg.exec_()

    
    def submit_set_of_bursts_clicked(self):
    
        bursts = self.set_of_bursts_lineedit.text()
        self.burst.clear()
        try:
            self.multi_automata = MultiwaveSolitonAutomata(self.my_graph_m, bursts)
            bursts = bursts.split(";")
            for burst in bursts:
                burst = re.sub(r"[{}]+", "", burst)
                self.burst.addItem(burst)
            self.show_multiple([self.mol_info_m, self.bursts_label, self.row5_m, self.submit_burst])
            self.hide_multiple([self.traversals_label, self.traversals, self.show_matrices_m, self.show_end_result_m, self.show_animation_m])
        except:
            self.hide_multiple([self.mol_info_m, self.bursts_label, self.row5_m, self.submit_burst, self.traversals_label, self.traversals, self.show_matrices_m, self.show_end_result_m, self.show_animation_m])
            msg = QMessageBox(self.wid_mult)
            msg.setStyleSheet(" QPushButton{ height: 32px; width: 130px;}")
            msg.setWindowTitle("Incorrect input")
            msg.setText("The syntax of your input is not correct.")
            msg.setIcon(QMessageBox.Warning)
            msg.setStandardButtons(QMessageBox.Retry)
            msg.setInformativeText("Please try again with another input string.")
            details = f"Reminder - this is how you define a set of bursts: \n"
            details = details + f"- Starting and end node of a soliton are embedded in round brackets (e.g. (3, 2))\n"
            details = details + f"- Solitons are seperated by '||'\n"
            details = details + f"- The number in front of a soliton's pair of exterior nodes defines how many timesteps later than the previous soliton this soliton enters the graph\n"
            details = details + f"- A set of bursts is embedded in braces, where individual bursts are seperated by ';'\n"
            details = details + "- Example: " + "{(3,1)||1(1,2); (3,2)||4(2,1)}\n"
            msg.setDetailedText(details)
            x = msg.exec_()


    def all_exterior_nodes_statechanged(self):
        if self.all_exterior_nodes.isChecked():
            self.hide_multiple([self.node_1, self.exterior_nodes_label2, self.node_2])
        else: self.show_multiple([self.node_1, self.exterior_nodes_label2, self.node_2])


    def all_bursts_statechanged(self):
        if self.all_bursts.isChecked():
            #self.hide_retain_space([self.burst])
            self.burst.hide()
        else: self.burst.show()


    def save_clicked(self):
        """Method that is called when user clicks button to save the graph visualisation.
        Opens a file dialog in which user can specify a path where the image should be saved.
        Only allows `.jpg`, `.png` and `.jpeg` file suffixes.
        """
        option = QtWidgets.QFileDialog.Options()
        name = QtWidgets.QFileDialog.getSaveFileName(self.wid_single, 'Save File', 'graph.jpg', 'Images (*.jpg *.png *.jpeg)', options = option)
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

    def save_clicked_m(self):
        option = QtWidgets.QFileDialog.Options()
        name = QtWidgets.QFileDialog.getSaveFileName(self.wid_mult, 'Save File', 'graph.jpg', 'Images (*.jpg *.png *.jpeg)', options = option)
        if name != ('', ''):
            path = name[0]
            # turn PIL Image into byte array 
            imgByteArr = io.BytesIO()
            self.graph_pic_m.save(imgByteArr, format='PNG')
            imgByteArr = imgByteArr.getvalue()
            # save image at specified path
            file = open(path, "wb")
            file.write(imgByteArr)
            file.close()


    def submit_exterior_nodes_clicked(self):
        """Method that is called when user clicks button to submit a pair of exterior nodes.
        Computes all possible soliton paths between the two chosen nodes.
        Informs the user if no soliton path exists between them.
        Otherwise all the necessary widgets for the user to choose a computed soliton path and look at further information on it are revealed.
        """
        self.path_index = None # we need this variable later in show_animation_clicked
        self.paths.clear()
        if self.all_exterior_nodes.isChecked():
            if self.automata is None:
                self.automata = SolitonAutomata(self.my_graph)
            key = self.automata.matrix_to_string(nx.to_numpy_array(self.my_graph.graph))
            self.found_paths = self.automata.states_plus_soliton_paths[key][1]
        else:
            node1 = int(self.node_1.currentText())
            node2 = int(self.node_2.currentText())
            miniautomata = MiniSolitonAutomata(self.my_graph, node1, node2)
            self.found_paths = miniautomata.soliton_paths
        if self.found_paths == []:
            self.hide_multiple([self.soliton_paths_label, self.paths, self.show_matrices, self.show_end_result, self.show_animation])
            msg = QMessageBox(self.wid_single)
            msg.setStyleSheet(" QPushButton{ height: 32px; width: 130px;}")
            msg.setWindowTitle("No path found")
            msg.setText("There exists no soliton path between these exterior nodes.")
            msg.setIcon(QMessageBox.Information)
            msg.setStandardButtons(QMessageBox.Retry)
            msg.setInformativeText("Please try again with different exterior nodes.")
            x = msg.exec_()
        else:
            self.soliton_paths_label.setText(f"Soliton paths ({len(self.found_paths)}):")
            self.show_multiple([self.soliton_paths_label, self.paths, self.show_matrices, self.show_end_result, self.show_animation])
            for soliton_path in self.found_paths:
                self.paths.addItem(str(soliton_path.path_for_user))

    
    def submit_burst_clicked(self):

        self.traversal_index = None # we need this variable later in show_animation_clicked
        self.traversals.clear()
        if self.all_bursts.isChecked():
            key = self.multi_automata.matrix_to_string(nx.to_numpy_array(self.my_graph_m.graph))
            self.found_traversals = self.multi_automata.states_plus_traversals[key][1]
            self.num_traversals_per_burst = self.multi_automata.states_plus_traversals[key][2]
        else:
            burst_index = int(self.burst.currentIndex())
            self.found_traversals = self.multi_automata.call_find_all_travs_given_burst(self.multi_automata.bursts_dicts[burst_index], self.my_graph_m)
            self.num_traversals_per_burst = None
        if self.found_traversals == []:
            self.hide_multiple([self.traversals_label, self.traversals, self.show_matrices_m, self.show_end_result_m, self.show_animation_m])
            msg = QMessageBox(self.wid_mult)
            msg.setStyleSheet(" QPushButton{ height: 32px; width: 130px;}")
            msg.setWindowTitle("No soliton paths found")
            msg.setText("With this burst(s), not all solitons can traverse the soliton graph successfully.")
            msg.setIcon(QMessageBox.Information)
            msg.setStandardButtons(QMessageBox.Retry)
            msg.setInformativeText("Please try again with different bursts.")
            x = msg.exec_()
        else:
            endless_loops = 0
            for traversal in self.found_traversals:
                if isinstance(traversal, Traversal): # if traversal is a real traversal and no endless loop
                    this_traversal = ""
                    for i, path in enumerate(traversal.traversal_for_user):
                        this_traversal = this_traversal + path
                        if i != len(traversal.traversal_for_user)-1:
                            this_traversal = this_traversal + ", "
                    self.traversals.addItem(this_traversal)
                else:
                    this_traversal = "LOOP! "
                    for i, path in enumerate(traversal[0].traversal_for_user):
                        this_traversal = this_traversal + path + " ..."
                        if i != len(traversal[0].traversal_for_user)-1:
                            this_traversal = this_traversal + ", "
                    self.traversals.addItem(this_traversal)
                    endless_loops += 1
            self.traversals_label.setText(f"Sets of paths ({len(self.found_traversals) - endless_loops}):")
            self.show_multiple([self.traversals_label, self.traversals, self.show_matrices_m, self.show_end_result_m, self.show_animation_m])
            sep_index = 0
            if self.num_traversals_per_burst:
                for n, num in enumerate(self.num_traversals_per_burst): # add seperators to distinct which traversals resulted from which burst
                    if num != 0 and n != len(self.num_traversals_per_burst)-1: # don't add if there is no traversal for this burst of if it's the last burst
                        sep_index += num
                        self.traversals.insertSeparator(sep_index)
                        self.found_traversals.insert(sep_index, None) # so the indices in found_traversals correspond to the indices of the combobox items
                        sep_index += 1 # indices got shifted by one because seperator was added


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
            name = QtWidgets.QFileDialog.getSaveFileName(self.wid_single, 'Save File', 'matrices.txt', 'Text files (*.txt)', options = option)
            if name != ('', ''): # only do the following if user clicked on save button (without this line the application closes with an error if save action is cancelled)
                path = name[0]
                file = open(path, "w")
                file.write(txt_text)
                file.close()

        index = self.paths.currentIndex()
        desired_path = self.found_paths[index]
    
        dlg = QDialog(self.wid_single)
        scrollArea = QScrollArea(dlg)
        scrollArea.setStyleSheet("font: 13pt Courier;")
        widget = QtWidgets.QWidget()
        vbox = QtWidgets.QVBoxLayout()
        if re.search(',', self.my_graph.way) is not None:
            txt_text = f"Way to soliton graph: {self.my_graph.way} \n"
        else:
            txt_text = f"Soliton graph: {self.my_graph.way} \n"
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
        scrollArea.setFixedSize(540, 405)
        save_button = QtWidgets.QPushButton(dlg)
        save_button.setStyleSheet("QPushButton {background-color: rgb(191, 207, 255); image: url(:/icons/save.svg);} QPushButton::pressed {background-color : rgb(132, 145, 193);}")
        save_button.setGeometry(QtCore.QRect(454, 359, 70, 30))
        save_button.clicked.connect(save_matrices)

        dlg.setWindowTitle("Adjacency Matrices")
        dlg.setFixedSize(545, 410)
        dlg.exec_()


    def show_matrices_clicked_m(self):

        def save_matrices():
            option = QtWidgets.QFileDialog.Options()
            name = QtWidgets.QFileDialog.getSaveFileName(self.wid_mult, 'Save File', 'matrices.txt', 'Text files (*.txt)', options = option)
            if name != ('', ''): # only do the following if user clicked on save button (without this line the application closes with an error if save action is cancelled)
                path = name[0]
                file = open(path, "w")
                file.write(txt_text)
                file.close()

        index = self.traversals.currentIndex()
        desired_traversal = self.found_traversals[index]
    
        dlg = QDialog(self.wid_mult)
        scrollArea = QScrollArea(dlg)
        scrollArea.setStyleSheet("font: 13pt Courier;")
        widget = QtWidgets.QWidget()
        vbox = QtWidgets.QVBoxLayout()
        if re.search(',', self.my_graph_m.way) is not None:
            txt_text = f"Way to soliton graph: {self.my_graph_m.way} \n"
        else:
            txt_text = f"Soliton graph: {self.my_graph_m.way} \n"
        txt_text = txt_text + f"Traversal: {desired_traversal.traversal_for_user} \n \n"
        # labelling of matrix depends on wether we have long node labels ("aa", "ab", ...) or short ones ("a", "b", ...)
        if (len(self.my_graph_m.labels) - len(self.my_graph_m.exterior_nodes)) > 26:
            matrix_label_horizontal = "    "
            long_labels = True
        else:
            matrix_label_horizontal = "   "
            long_labels = False
        for key in self.my_graph_m.labels:
            if long_labels:
                if self.my_graph_m.labels[key] in self.my_graph_m.exterior_nodes_reverse and int(self.my_graph_m.labels[key]) < 10:
                    matrix_label_horizontal = matrix_label_horizontal + f"{self.my_graph_m.labels[key]}  "
                else:
                    matrix_label_horizontal = matrix_label_horizontal + f"{self.my_graph_m.labels[key]} "
            else:
                if self.my_graph_m.labels[key] in self.my_graph_m.exterior_nodes_reverse and int(self.my_graph_m.labels[key]) >= 10:
                    matrix_label_horizontal = matrix_label_horizontal + f"{self.my_graph_m.labels[key]} "
                else:
                    matrix_label_horizontal = matrix_label_horizontal + f"{self.my_graph_m.labels[key]}  "
        # for show-matrices-window the node labels are added to the matrix string and this labelled matrix is added to scroll area
        # for matrices.txt we add everything (horizontal label and every row of matrix) line by line
        for i in range(len(desired_traversal.adjacency_matrices_list)):
            txt_text = txt_text + f"Timestep {i}: \n"
            txt_text = txt_text + f"{matrix_label_horizontal} \n"
            matrix_labelled = ""
            matrix_labelled = matrix_labelled + f"{matrix_label_horizontal}\n"
            matrix = str(desired_traversal.adjacency_matrices_list[i])
            matrix = re.sub(r"[matrix(]", "", matrix)
            matrix = re.sub(r"[)]", "", matrix)
            for j in range(len(matrix.splitlines())):
                if long_labels:
                    if self.my_graph_m.labels[j] in self.my_graph_m.exterior_nodes_reverse and int(self.my_graph_m.labels[j]) < 10:
                        matrix_labelled = matrix_labelled + f"{self.my_graph_m.labels[j]} {matrix.splitlines()[j]} \n"
                        txt_text = txt_text + f"{self.my_graph_m.labels[j]} {matrix.splitlines()[j]} \n"
                    else:
                        matrix_labelled = matrix_labelled + f"{self.my_graph_m.labels[j]}{matrix.splitlines()[j]} \n"
                        txt_text = txt_text + f"{self.my_graph_m.labels[j]}{matrix.splitlines()[j]} \n"
                else:
                    if self.my_graph_m.labels[j] in self.my_graph_m.exterior_nodes_reverse and int(self.my_graph_m.labels[j]) >= 10:
                        matrix_labelled = matrix_labelled + f"{self.my_graph_m.labels[j]}{matrix.splitlines()[j][1:]} \n"
                        txt_text = txt_text + f"{self.my_graph_m.labels[j]}{matrix.splitlines()[j][1:]} \n"
                    else:
                        matrix_labelled = matrix_labelled + f"{self.my_graph_m.labels[j]}{matrix.splitlines()[j]} \n"
                        txt_text = txt_text + f"{self.my_graph_m.labels[j]}{matrix.splitlines()[j]} \n"
            txt_text = txt_text + f"\n"
            vbox.addWidget(QtWidgets.QLabel(matrix_labelled))
        widget.setLayout(vbox)
        scrollArea.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        scrollArea.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        scrollArea.setWidgetResizable(True)
        scrollArea.setWidget(widget)
        scrollArea.setFixedSize(540, 405)
        save_button = QtWidgets.QPushButton(dlg)
        save_button.setStyleSheet("QPushButton {background-color: rgb(149, 221, 185); image: url(:/icons/save.svg);} QPushButton::pressed {background-color : rgb(90, 159, 123);}")
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
            name = QtWidgets.QFileDialog.getSaveFileName(self.wid_single, 'Save File', 'result.jpg', 'Images (*.jpg *.png *.jpeg)', options = option)
            if name != ('', ''):
                path = name[0]
                imgByteArr = io.BytesIO()
                result_pic.save(imgByteArr, format='PNG')
                imgByteArr = imgByteArr.getvalue()
                file = open(path, "wb")
                file.write(imgByteArr)
                file.close()

        def use_as_new_soliton_graph(dlg: QDialog):
            """Uses the end result of the selected path as the new soliton graph. Puts the software in traversal mode.

            Args:
                dlg (QDialog): The dialog window that shows the end result.
            """
            current_way = self.my_graph.way
            new_soliton_graph = copy.deepcopy(self.desired_path.resulting_soliton_graph)
            self.my_graph = new_soliton_graph
            self.my_graph.way = f"{current_way}, {self.desired_path.path_for_user}"
            self.graph_pic = Visualisation.visualize_soliton_graph(self.my_graph, self.my_graph.bindings, False, True)
            self.qim = ImageQt(self.graph_pic)
            self.display_molecule.setPixmap(QtGui.QPixmap.fromImage(self.qim).scaled(self.display_molecule.size(), QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation))

            self.traversal_mode.setChecked(True)
            self.hide_multiple([self.soliton_paths_label, self.paths, self.show_matrices, self.show_end_result, self.show_animation])
            self.node_1.setCurrentIndex(0)
            self.node_2.setCurrentIndex(0)
            dlg.close()

        index = self.paths.currentIndex()
        self.desired_path = self.found_paths[index]
        bindings_index = len(self.desired_path.path) - 1
        result_pic = Visualisation.visualize_soliton_graph(self.my_graph, self.desired_path.bindings_list[bindings_index], False, True)
        qim = ImageQt(result_pic)

        dlg = QDialog(self.wid_single)
        label = QtWidgets.QLabel(dlg)
        label.setGeometry(QtCore.QRect(0, 0, 540, 405))
        label.setPixmap(QtGui.QPixmap.fromImage(qim))
        label.setScaledContents(True)
        save_button = QtWidgets.QPushButton(dlg)
        save_button.setStyleSheet("QPushButton {background-color: rgb(191, 207, 255); image: url(:/icons/save.svg);} QPushButton::pressed {background-color : rgb(132, 145, 193);}")
        save_button.setGeometry(QtCore.QRect(470, 375, 70, 30))
        save_button.clicked.connect(save_end_result)
        use_button = QtWidgets.QPushButton("Use", dlg)
        use_button.clicked.connect(lambda: use_as_new_soliton_graph(dlg))
        use_button.setGeometry(QtCore.QRect(395, 375, 70, 30))

        dlg.setWindowTitle("End result")
        dlg.setFixedSize(545, 410)
        dlg.exec_()

    
    def show_end_result_clicked_m(self):

        def save_end_result():

            option = QtWidgets.QFileDialog.Options()
            name = QtWidgets.QFileDialog.getSaveFileName(self.wid_mult, 'Save File', 'result.jpg', 'Images (*.jpg *.png *.jpeg)', options = option)
            if name != ('', ''):
                path = name[0]
                imgByteArr = io.BytesIO()
                result_pic.save(imgByteArr, format='PNG')
                imgByteArr = imgByteArr.getvalue()
                file = open(path, "wb")
                file.write(imgByteArr)
                file.close()

        def use_as_new_soliton_graph(dlg: QDialog):
            
            current_way = self.my_graph_m.way
            new_soliton_graph = copy.deepcopy(self.desired_traversal.resulting_soliton_graph)
            self.my_graph_m = new_soliton_graph
            self.my_graph_m.way = f"{current_way}, {self.desired_traversal.traversal_for_user}"
            self.graph_pic = Visualisation.visualize_soliton_graph(self.my_graph_m, self.my_graph_m.bindings, False, True)
            self.qim_m = ImageQt(self.graph_pic)
            self.display_molecule_m.setPixmap(QtGui.QPixmap.fromImage(self.qim_m).scaled(self.display_molecule_m.size(), QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation))

            self.traversal_mode_m.setChecked(True)
            self.hide_multiple([self.traversals_label, self.traversals, self.show_matrices_m, self.show_end_result_m, self.show_animation_m])
            self.burst.setCurrentIndex(0)
            dlg.close()

        index = self.traversals.currentIndex()
        self.desired_traversal = self.found_traversals[index]
        result_pic = Visualisation.visualize_soliton_graph(self.desired_traversal.resulting_soliton_graph, self.desired_traversal.resulting_soliton_graph.bindings, False, True)
        qim = ImageQt(result_pic)

        dlg = QDialog(self.wid_mult)
        label = QtWidgets.QLabel(dlg)
        label.setGeometry(QtCore.QRect(0, 0, 540, 405))
        label.setPixmap(QtGui.QPixmap.fromImage(qim))
        label.setScaledContents(True)
        save_button = QtWidgets.QPushButton(dlg)
        save_button.setStyleSheet("QPushButton {background-color: rgb(149, 221, 185); image: url(:/icons/save.svg);} QPushButton::pressed {background-color : rgb(90, 159, 123);}")
        save_button.setGeometry(QtCore.QRect(470, 375, 70, 30))
        save_button.clicked.connect(save_end_result)
        use_button = QtWidgets.QPushButton("Use", dlg)
        use_button.clicked.connect(lambda: use_as_new_soliton_graph(dlg))
        use_button.setGeometry(QtCore.QRect(395, 375, 70, 30))

        dlg.setWindowTitle("End result")
        dlg.setFixedSize(545, 410)
        dlg.exec_()


    def show_animation_clicked(self, button: QtWidgets.QPushButton):
        """Method that is called when user clicks button to have the animation of the soliton traversing the graph displayed.
        Makes a small window pop up that shows the animation and provides a save button.
        Instead of displaying the `gif` it uses a sequence of `PIL` images and always shows the next image after a certain time.
        """

        def save_animation():
            """Opens a file dialog in which user can specify a path where the animation should be saved.
            Only allows `.gif` file suffix.
            """
            option = QtWidgets.QFileDialog.Options()
            name = QtWidgets.QFileDialog.getSaveFileName(wid, 'Save File', 'animation.gif', 'Images (*.gif)', options = option)
            if name != ('', ''):
                path = name[0]
                if button == self.show_animation:
                    ani = Animation.graph_animation(self.my_graph, self.desired_path)
                else: ani = Animation.graph_animation_multiwave(self.my_graph_m, self.desired_traversal)
                ani.save(path, writer='pillow', dpi = 600)

        def update_image():
            """Displays next image of animation after a certain time (method is triggered by a timer).
            """
            self.step += 1
            if button == self.show_animation:
                if self.step == len(self.desired_path.path): # start animation all over again as soon as end is reached (endless loop)
                    self.step = 0
            else:
                if isinstance(self.desired_traversal, Traversal):
                    if self.step == len(self.desired_traversal.pos):
                        self.step = 0
                else:
                    if self.step == len(self.desired_traversal[0].pos):
                        self.over_looppoint_count += 1 # another "round" of loop made
                        self.step = self.desired_traversal[1] # we reached the "end" (which in reality is the point where we noticed we are stuck in a loop) so we have to continue at the loop point now
            im = self.pil_images[self.step]
            qim = ImageQt(im)
            self.label.setPixmap(QtGui.QPixmap.fromImage(qim.copy()))
            self.label.setScaledContents(True)

        def update_image_reverse():
            self.step -= 1
            if button == self.show_animation:
                if self.step == -1: # if start of animation was reached in the step before then continue at end of animation
                    self.step = len(self.desired_path.path) - 1
            else:
                if isinstance(self.desired_traversal, Traversal):
                    if self.step == -1:
                        self.step = len(self.desired_traversal.pos) - 1
                else:
                    if self.step == -1:
                        self.step = 0 # we can't go back to the end of the animation because it has no end, instead we stay at the beginning
                    elif self.step == self.desired_traversal[1] - 1 and self.over_looppoint_count != 0: # before we go back to before the loop point we have to dismantle all "rounds" of loop that we made
                        self.step = len(self.desired_traversal[0].pos) - 1 # continue again at the "end"
                        self.over_looppoint_count -= 1 # dismantle
            im = self.pil_images[self.step]
            qim = ImageQt(im)
            self.label.setPixmap(QtGui.QPixmap.fromImage(qim.copy()))
            self.label.setScaledContents(True)
            

        def pause_animation(button: QtWidgets.QPushButton):
            if self.timer.isActive():
                self.timer.stop()
                if button.parentWidget().parentWidget() == self.wid_single:
                    button.setStyleSheet("QPushButton {background-color: rgb(191, 207, 255); image: url(:/icons/play.svg);} QPushButton::pressed {background-color : rgb(132, 145, 193);}")
                else: button.setStyleSheet("QPushButton {background-color: rgb(149, 221, 185); image: url(:/icons/play.svg);} QPushButton::pressed {background-color : rgb(90, 159, 123);}")
            else:
                update_image()
                self.timer.start(800)
                if button.parentWidget().parentWidget() == self.wid_single:
                    button.setStyleSheet("QPushButton {background-color: rgb(191, 207, 255); image: url(:/icons/pause.svg);} QPushButton::pressed {background-color : rgb(132, 145, 193);}")
                else: button.setStyleSheet("QPushButton {background-color: rgb(149, 221, 185); image: url(:/icons/pause.svg);} QPushButton::pressed {background-color : rgb(90, 159, 123);}")

        def next_img(button: QtWidgets.QPushButton):
            if self.timer.isActive():
                self.timer.stop()
                if button.parentWidget().parentWidget() == self.wid_single:
                    button.setStyleSheet("QPushButton {background-color: rgb(191, 207, 255); image: url(:/icons/play.svg);} QPushButton::pressed {background-color : rgb(132, 145, 193);}")
                else: button.setStyleSheet("QPushButton {background-color: rgb(149, 221, 185); image: url(:/icons/play.svg);} QPushButton::pressed {background-color : rgb(90, 159, 123);}")
            update_image()

        def prev_img(button: QtWidgets.QPushButton):
            if self.timer.isActive():
                self.timer.stop()
                if button.parentWidget().parentWidget() == self.wid_single:
                    button.setStyleSheet("QPushButton {background-color: rgb(191, 207, 255); image: url(:/icons/play.svg);} QPushButton::pressed {background-color : rgb(132, 145, 193);}")
                else: button.setStyleSheet("QPushButton {background-color: rgb(149, 221, 185); image: url(:/icons/play.svg);} QPushButton::pressed {background-color : rgb(90, 159, 123);}")
            update_image_reverse()
            
        if button == self.show_animation:
            dlg = QDialog(self.wid_single)
        else: dlg = QDialog(self.wid_mult)
        self.label = QtWidgets.QLabel(dlg)
        self.label.setGeometry(QtCore.QRect(0, 0, 540, 405))
        
        if button == self.show_animation:
            wid = self.wid_single
        else:
            wid = self.wid_mult
        
        if button == self.show_animation:
            if self.path_index != self.paths.currentIndex():
                self.path_index = self.paths.currentIndex()
                self.desired_path = self.found_paths[self.path_index]
                plots_and_arrays = Animation.list_of_plots_and_arrays(self.my_graph, self.desired_path)
                self.pil_images = Animation.list_of_pil_images(plots_and_arrays)
        else:
            if self.traversal_index != self.traversals.currentIndex():
                self.traversal_index = self.traversals.currentIndex()
                self.desired_traversal = self.found_traversals[self.traversal_index]
                if isinstance(self.desired_traversal, Traversal):
                    plots_and_arrays = Animation.list_of_plots_and_arrays_multiwave(self.my_graph_m, self.desired_traversal)
                else: 
                    plots_and_arrays = Animation.list_of_plots_and_arrays_multiwave(self.my_graph_m, self.desired_traversal[0])
                self.pil_images = Animation.list_of_pil_images(plots_and_arrays)

        save_button = QtWidgets.QPushButton(dlg)
        save_button.setGeometry(QtCore.QRect(470, 375, 70, 30))
        if isinstance(self.desired_traversal, Traversal):
            save_button.clicked.connect(save_animation)
        else:
            self.over_looppoint_count = 0 # how many times we passed the loop point
        pause_button = QtWidgets.QPushButton(dlg)
        pause_button.setGeometry(QtCore.QRect(240, 375, 30, 30))
        pause_button.clicked.connect(lambda: pause_animation(pause_button))
        next_button = QtWidgets.QPushButton(dlg)
        next_button.setGeometry(QtCore.QRect(275, 375, 30, 30))
        next_button.clicked.connect(lambda: next_img(pause_button))
        prev_button = QtWidgets.QPushButton(dlg)
        prev_button.setGeometry(QtCore.QRect(205, 375, 30, 30))
        prev_button.clicked.connect(lambda: prev_img(pause_button))
        if button == self.show_animation:
            save_button.setStyleSheet("QPushButton {background-color: rgb(191, 207, 255); image: url(:/icons/save.svg);} QPushButton::pressed {background-color : rgb(132, 145, 193);}")
            pause_button.setStyleSheet("QPushButton {background-color: rgb(191, 207, 255); image: url(:/icons/pause.svg);} QPushButton::pressed {background-color : rgb(132, 145, 193);}")
            next_button.setStyleSheet("QPushButton {background-color: rgb(191, 207, 255); image: url(:/icons/right-arrow.svg);} QPushButton::pressed {background-color : rgb(132, 145, 193);}")
            prev_button.setStyleSheet("QPushButton {background-color: rgb(191, 207, 255); image: url(:/icons/left-arrow.svg);} QPushButton::pressed {background-color : rgb(132, 145, 193)}")
        else:
            if isinstance(self.desired_traversal, Traversal):
                save_button.setStyleSheet("QPushButton {background-color: rgb(149, 221, 185); image: url(:/icons/save.svg);} QPushButton::pressed {background-color : rgb(90, 159, 123);}")
            else:
                save_button.setStyleSheet("QPushButton {background-color: rgb(230, 230, 230); image: url(:/icons/save.svg);}")
            pause_button.setStyleSheet("QPushButton {background-color: rgb(149, 221, 185); image: url(:/icons/pause.svg);} QPushButton::pressed {background-color : rgb(90, 159, 123);}")
            next_button.setStyleSheet("QPushButton {background-color: rgb(149, 221, 185); image: url(:/icons/right-arrow.svg);} QPushButton::pressed {background-color : rgb(90, 159, 123);}")
            prev_button.setStyleSheet("QPushButton {background-color: rgb(149, 221, 185); image: url(:/icons/left-arrow.svg);} QPushButton::pressed {background-color : rgb(90, 159, 123)}")

        self.step = -1
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(update_image)
        update_image()
        self.timer.start(800) # triggers event every 800 millisecond

        dlg.setWindowTitle("Animation")
        dlg.setFixedSize(545, 410)
        dlg.closeEvent = self.stop_animation # stop timer when window is closed
        dlg.exec_()


    def stop_animation(self, event):
        """Stops the timer that is used for the animation.
        Without this method the application would crash if the animation window would get closed.
        Args:
            event: Close event of the window that shows animation.
        """
        self.timer.disconnect()


    def hide_retain_space(self, widgets: list):
        """Retains the space of a widget even when it's hidden.

        Args:
            widget (QtWidgets.QWidget ): Widget whose space should be retained.
        """
        for widget in widgets:
            retain = widget.sizePolicy()
            retain.setRetainSizeWhenHidden(True)
            widget.setSizePolicy(retain)


    def hide_multiple(self, widgets: list):
        for widget in widgets:
            widget.hide()


    def show_multiple(self, widgets: list):
        for widget in widgets:
            widget.show()


    def heightForWidth(self, width: float):
        """Computes height for a given width.

        Args:
            width (float): Given width.

        Returns:
            float: Computed height.
        """
        height_for_width_factor = 1.0 * 519 / 688
        return math.ceil(width * height_for_width_factor)


    def widthForHeight(self, height):
        """Computes width for a given height.

        Args:
            height (float): Given height.

        Returns:
            float: Computed width.
        """
        width_for_height_factor = 1.0 * 688 / 519
        return math.ceil(height * width_for_height_factor)


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
        self.display_molecule_m.setPixmap(QtGui.QPixmap.fromImage(self.qim_m).scaled(width, height, QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation))
