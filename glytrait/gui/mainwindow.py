# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'gui.ui'
##
## Created by: Qt User Interface Compiler version 6.5.1
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import QCoreApplication, QMetaObject, QRect, Qt
from PySide6.QtGui import QFont, QPixmap
from PySide6.QtWidgets import (
    QCheckBox,
    QFrame,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMenuBar,
    QPushButton,
    QStatusBar,
    QVBoxLayout,
    QWidget,
)


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName("MainWindow")
        MainWindow.resize(650, 600)
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.groupBox = QGroupBox(self.centralwidget)
        self.groupBox.setObjectName("groupBox")
        self.groupBox.setGeometry(QRect(20, 10, 601, 80))
        self.glytrait_label = QLabel(self.groupBox)
        self.glytrait_label.setObjectName("glytrait_label")
        self.glytrait_label.setGeometry(QRect(10, 10, 71, 31))
        font = QFont()
        font.setPointSize(18)
        font.setBold(True)
        self.glytrait_label.setFont(font)
        self.introduction = QLabel(self.groupBox)
        self.introduction.setObjectName("introduction")
        self.introduction.setGeometry(QRect(10, 50, 271, 16))
        self.logo = QLabel(self.groupBox)
        self.logo.setObjectName("logo")
        self.logo.setGeometry(QRect(270, 10, 321, 61))
        self.gridLayoutWidget = QWidget(self.centralwidget)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayoutWidget.setGeometry(QRect(20, 140, 601, 204))
        self.gridLayout = QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setObjectName("gridLayout")
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.input_label = QLabel(self.gridLayoutWidget)
        self.input_label.setObjectName("input_label")
        font1 = QFont()
        font1.setPointSize(14)
        font1.setBold(True)
        self.input_label.setFont(font1)

        self.gridLayout.addWidget(self.input_label, 0, 1, 1, 1)

        self.sia_checkbox = QCheckBox(self.gridLayoutWidget)
        self.sia_checkbox.setObjectName("sia_checkbox")

        self.gridLayout.addWidget(self.sia_checkbox, 6, 1, 1, 1)

        self.output_label = QLabel(self.gridLayoutWidget)
        self.output_label.setObjectName("output_label")
        self.output_label.setFont(font1)

        self.gridLayout.addWidget(self.output_label, 3, 1, 1, 1)

        self.output_line_edit = QLineEdit(self.gridLayoutWidget)
        self.output_line_edit.setObjectName("output_line_edit")
        font2 = QFont()
        font2.setPointSize(14)
        font2.setBold(False)
        self.output_line_edit.setFont(font2)

        self.gridLayout.addWidget(self.output_line_edit, 4, 1, 1, 1)

        self.input_line_edit = QLineEdit(self.gridLayoutWidget)
        self.input_line_edit.setObjectName("input_line_edit")

        self.gridLayout.addWidget(self.input_line_edit, 1, 1, 1, 1)

        self.sia_label = QLabel(self.gridLayoutWidget)
        self.sia_label.setObjectName("sia_label")
        self.sia_label.setFont(font1)

        self.gridLayout.addWidget(self.sia_label, 5, 1, 1, 1)

        self.input_button = QPushButton(self.gridLayoutWidget)
        self.input_button.setObjectName("input_button")

        self.gridLayout.addWidget(self.input_button, 1, 2, 1, 1)

        self.output_button = QPushButton(self.gridLayoutWidget)
        self.output_button.setObjectName("output_button")

        self.gridLayout.addWidget(self.output_button, 4, 2, 1, 1)

        self.basic_label = QLabel(self.centralwidget)
        self.basic_label.setObjectName("basic_label")
        self.basic_label.setGeometry(QRect(20, 100, 101, 21))
        font3 = QFont()
        font3.setPointSize(15)
        font3.setBold(True)
        self.basic_label.setFont(font3)
        self.line = QFrame(self.centralwidget)
        self.line.setObjectName("line")
        self.line.setGeometry(QRect(20, 120, 601, 20))
        self.line.setFrameShape(QFrame.HLine)
        self.line.setFrameShadow(QFrame.Sunken)
        self.advance_label = QLabel(self.centralwidget)
        self.advance_label.setObjectName("advance_label")
        self.advance_label.setGeometry(QRect(20, 360, 141, 21))
        self.advance_label.setFont(font3)
        self.line_2 = QFrame(self.centralwidget)
        self.line_2.setObjectName("line_2")
        self.line_2.setGeometry(QRect(20, 380, 601, 20))
        self.line_2.setFrameShape(QFrame.HLine)
        self.line_2.setFrameShadow(QFrame.Sunken)
        self.horizontalLayoutWidget = QWidget(self.centralwidget)
        self.horizontalLayoutWidget.setObjectName("horizontalLayoutWidget")
        self.horizontalLayoutWidget.setGeometry(QRect(20, 400, 601, 80))
        self.horizontalLayout = QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_2 = QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.save_template_label = QLabel(self.horizontalLayoutWidget)
        self.save_template_label.setObjectName("save_template_label")
        self.save_template_label.setFont(font1)

        self.verticalLayout_2.addWidget(self.save_template_label)

        self.save_template_button = QPushButton(self.horizontalLayoutWidget)
        self.save_template_button.setObjectName("save_template_button")

        self.verticalLayout_2.addWidget(self.save_template_button)

        self.horizontalLayout.addLayout(self.verticalLayout_2)

        self.line_4 = QFrame(self.horizontalLayoutWidget)
        self.line_4.setObjectName("line_4")
        self.line_4.setFrameShape(QFrame.VLine)
        self.line_4.setFrameShadow(QFrame.Sunken)

        self.horizontalLayout.addWidget(self.line_4)

        self.gridLayout_2 = QGridLayout()
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.use_template_label = QLabel(self.horizontalLayoutWidget)
        self.use_template_label.setObjectName("use_template_label")
        self.use_template_label.setFont(font1)

        self.gridLayout_2.addWidget(self.use_template_label, 0, 0, 1, 1)

        self.formula_edit_line = QLineEdit(self.horizontalLayoutWidget)
        self.formula_edit_line.setObjectName("formula_edit_line")

        self.gridLayout_2.addWidget(self.formula_edit_line, 1, 0, 1, 1)

        self.formula_button = QPushButton(self.horizontalLayoutWidget)
        self.formula_button.setObjectName("formula_button")

        self.gridLayout_2.addWidget(self.formula_button, 1, 1, 1, 1)

        self.horizontalLayout.addLayout(self.gridLayout_2)

        self.groupBox_2 = QGroupBox(self.centralwidget)
        self.groupBox_2.setObjectName("groupBox_2")
        self.groupBox_2.setGeometry(QRect(20, 525, 625, 31))
        self.run_button = QPushButton(self.groupBox_2)
        self.run_button.setObjectName("run_button")
        self.run_button.setGeometry(QRect(200, 0, 201, 32))
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QMenuBar(MainWindow)
        self.menubar.setObjectName("menubar")
        self.menubar.setGeometry(QRect(0, 0, 650, 24))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)

        QMetaObject.connectSlotsByName(MainWindow)

    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(
            QCoreApplication.translate("MainWindow", "MainWindow", None)
        )
        self.groupBox.setTitle("")
        self.glytrait_label.setText(
            QCoreApplication.translate("MainWindow", "GlyTrait", None)
        )
        self.introduction.setText(
            QCoreApplication.translate(
                "MainWindow", "A tool to calculate glycan traits.", None
            )
        )
        self.logo.setText("")
        self.logo.setAlignment(Qt.AlignRight)
        pixmap = QPixmap("logo.png").scaled(
            self.logo.size(), aspectMode=Qt.KeepAspectRatio
        )
        self.logo.setPixmap(pixmap)
        self.logo.repaint()

        self.input_label.setText(
            QCoreApplication.translate("MainWindow", "Input file", None)
        )
        self.sia_checkbox.setText(
            QCoreApplication.translate(
                "MainWindow", "Include the sialic-acid-linkage traits", None
            )
        )
        self.output_label.setText(
            QCoreApplication.translate("MainWindow", "Output file", None)
        )
        self.output_line_edit.setText(
            QCoreApplication.translate(
                "MainWindow", "Default (the same directory as the input file)", None
            )
        )
        self.sia_label.setText(
            QCoreApplication.translate("MainWindow", "Sialic acid linkage", None)
        )
        self.input_button.setText(
            QCoreApplication.translate("MainWindow", "select", None)
        )
        self.output_button.setText(
            QCoreApplication.translate("MainWindow", "select", None)
        )
        self.basic_label.setText(
            QCoreApplication.translate("MainWindow", "Basic options", None)
        )
        self.advance_label.setText(
            QCoreApplication.translate("MainWindow", "Advanced options", None)
        )
        self.save_template_label.setText(
            QCoreApplication.translate("MainWindow", "Save template file", None)
        )
        self.save_template_button.setText(
            QCoreApplication.translate("MainWindow", "Save to ...", None)
        )
        self.use_template_label.setText(
            QCoreApplication.translate("MainWindow", "Use formula file", None)
        )
        self.formula_button.setText(
            QCoreApplication.translate("MainWindow", "select", None)
        )
        self.groupBox_2.setTitle("")
        self.run_button.setText(QCoreApplication.translate("MainWindow", "Run", None))

    # retranslateUi
