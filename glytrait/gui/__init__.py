import sys
from pathlib import Path
import platform

from PySide6.QtCore import Slot
from PySide6.QtWidgets import QApplication, QMainWindow, QFileDialog

from glytrait.trait import save_trait_formula_template
from glytrait.core import run_workflow
from glytrait.gui.mainwindow_mac import Ui_MainWindow as Ui_MainWindow_mac
from glytrait.gui.mainwindow_win import Ui_MainWindow as Ui_MainWindow_win


def get_os():
    os_name = platform.system()
    if os_name == "Darwin":
        return "MacOS"
    elif os_name == "Windows":
        return "Windows"
    else:
        return "Unknown"


class MainWindow(QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()
        if get_os() == "Windows":
            self.ui = Ui_MainWindow_win()
        elif get_os() == "MacOS":
            self.ui = Ui_MainWindow_mac()
        else:
            raise OSError("Unsupported OS.")
        self.ui.setupUi(self)
        self.setWindowTitle("GlyTrait")
        self.setFixedSize(650, 600)

        self.ui.input_button.clicked.connect(self.input_button_clicked)
        self.ui.output_button.clicked.connect(self.output_button_clicked)
        self.ui.save_template_button.clicked.connect(self.save_template_button_clicked)
        self.ui.formula_button.clicked.connect(self.formula_button_clicked)
        self.ui.run_button.clicked.connect(self.run)

    @Slot()
    def input_button_clicked(self):
        filename, _ = QFileDialog.getOpenFileName(
            self, "Open File", "", "CSV Files (*.csv)"
        )
        self.ui.input_line_edit.setText(filename)

    @Slot()
    def output_button_clicked(self):
        filename, _ = QFileDialog.getSaveFileName(
            self, "Save File", "", "XLSX Files (*.xlsx)"
        )
        self.ui.output_line_edit.setText(filename)

    @Slot()
    def save_template_button_clicked(self):
        # Choose a directory to save the template
        dirname = QFileDialog.getExistingDirectory(
            self, "Choose a directory to save the template"
        )
        save_trait_formula_template(dirname)
        msg = (f"Template saved to {Path(dirname) / 'trait_formula.txt'}.")
        self.statusBar().showMessage(msg)

    @Slot()
    def formula_button_clicked(self):
        filename, _ = QFileDialog.getOpenFileName(
            self, "Open File", "", "TXT Files (*.txt)"
        )
        self.ui.formula_edit_line.setText(filename)

    @Slot()
    def run(self):
        input_file = self.ui.input_line_edit.text()
        output_file = self.ui.output_line_edit.text()
        if "Default" in output_file:
            output_filename = Path(input_file).stem + "_glytrait.xlsx"
            output_file = str(Path(input_file).parent / output_filename)
        formula_file = self.ui.formula_edit_line.text()
        if formula_file == "":
            formula_file = None
        sia_linkage = self.ui.sia_checkbox.isChecked()
        run_workflow(input_file, output_file, sia_linkage, formula_file)
        msg = f"Done! Results saved to {output_file}."
        self.statusBar().showMessage(msg)


if __name__ == "__main__":
    app = QApplication(sys.argv)

    window = MainWindow()
    window.show()

    sys.exit(app.exec())
