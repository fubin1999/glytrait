import platform
import sys
from pathlib import Path

from PySide6.QtCore import Slot
from PySide6.QtWidgets import QApplication, QMainWindow, QFileDialog, QMessageBox

from glytrait.core import run_workflow
from glytrait.exception import GlyTraitError
from glytrait.gui.mainwindow_mac import Ui_MainWindow as Ui_MainWindow_mac
from glytrait.gui.mainwindow_win import Ui_MainWindow as Ui_MainWindow_win
from glytrait.trait import save_trait_formula_template


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
        if dirname == "":
            return
        save_trait_formula_template(dirname)
        msg = f"Template saved to {Path(dirname) / 'trait_formula.txt'}."
        self.pop_message_box(text="Done!", infor_text=msg)

    @Slot()
    def formula_button_clicked(self):
        filename, _ = QFileDialog.getOpenFileName(
            self, "Open File", "", "TXT Files (*.txt)"
        )
        self.ui.formula_edit_line.setText(filename)

    @Slot()
    def run(self):
        input_file = self.get_input_filename()
        output_file = self.get_output_filename()
        formula_file = self.get_formula_filename()
        sia_linkage = self.get_sia_linkage()

        if not self.valid_input_filename(input_file):
            self.pop_message_box(
                text="Please select a valid input file!",
                infor_text="A .csv file is needed.",
            )
            return

        if not self.valid_output_filename(output_file):
            self.pop_message_box(
                text="Please select a valid output file!",
                infor_text="A .xlsx file is needed.",
            )
            return

        try:
            run_workflow(input_file, output_file, sia_linkage, formula_file)
        except GlyTraitError as e:
            self.pop_message_box(
                text="Error!", infor_text=str(e), icon=QMessageBox.Icon.Critical
            )
        else:
            msg = f"Results saved to {output_file}."
            self.pop_message_box(text="Done!", infor_text=msg)

    def get_input_filename(self):
        return self.ui.input_line_edit.text()

    def get_output_filename(self):
        input_file = self.ui.input_line_edit.text()
        output_file = self.ui.output_line_edit.text()
        if "Default" in output_file:
            output_filename = Path(input_file).stem + "_glytrait.xlsx"
            output_file = str(Path(input_file).parent / output_filename)
        return output_file

    def get_formula_filename(self):
        formula_file = self.ui.formula_edit_line.text()
        return formula_file if formula_file != "" else None

    def get_sia_linkage(self):
        return self.ui.sia_checkbox.isChecked()

    def valid_filename(self, filename: str | None, suffix: str):
        if filename == "":
            return False
        if not Path(filename).exists():
            return False
        if not Path(filename).is_file():
            return False
        if not Path(filename).suffix == suffix:
            return False
        return True

    def valid_input_filename(self, filename):
        return self.valid_filename(filename, ".csv")

    def valid_output_filename(self, filename):
        return self.valid_filename(filename, ".xlsx")

    def valid_formula_filename(self, filename):
        return self.valid_filename(filename, ".txt")

    def pop_message_box(self, text, infor_text, icon=QMessageBox.Icon.NoIcon):
        msgBox = QMessageBox()
        msgBox.setIcon(icon)
        msgBox.setWindowTitle("GlyTrait")
        msgBox.setText(text)
        msgBox.setInformativeText(infor_text)
        msgBox.exec()


if __name__ == "__main__":
    app = QApplication(sys.argv)

    window = MainWindow()
    window.show()

    sys.exit(app.exec())
