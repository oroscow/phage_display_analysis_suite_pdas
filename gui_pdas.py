import sys
import os
import json
from PyQt6.QtWidgets import QApplication, QLabel, QMainWindow, QVBoxLayout, QWidget, QHBoxLayout, QLineEdit, QFormLayout, QTextEdit, QFileDialog, QPushButton, QComboBox
from PyQt6.QtGui import QFont
from seq_analysis import PdasSeqAnalysis


class PdasWindow(QMainWindow):
    """Main window."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("PDAS")
        self.json_data = None

        # Define layouts
        mainLayout = QVBoxLayout()
        topLayout = QVBoxLayout()
        bottom_layout = QFormLayout()

        # Create lines and buttons
        self.folder_path_edit = QLineEdit()
        self.library_file_edit = QLineEdit()
        self.library_choice_combo = QComboBox()
        browse_folders_button = QPushButton("Browse...")
        browse_files_button = QPushButton("Browse...")
        run_button = QPushButton("Run")
        exit_button = QPushButton("Exit")
        button_layout = QHBoxLayout()
        button_layout.addWidget(run_button)
        button_layout.addWidget(exit_button)

        # Connect buttons
        browse_folders_button.clicked.connect(
            lambda: self.browse_folders(self.folder_path_edit))
        browse_files_button.clicked.connect(
            lambda: self.browse_files(self.library_file_edit))
        run_button.clicked.connect(self.run_button_click)
        exit_button.clicked.connect(self.close)

        # Connect combo box selection change
        self.library_choice_combo.currentIndexChanged.connect(self.update_display)

        # Create horizontal layouts
        folder_layout = QHBoxLayout()
        folder_layout.addWidget(self.folder_path_edit)
        folder_layout.addWidget(browse_folders_button)

        library_layout = QHBoxLayout()
        library_layout.addWidget(self.library_file_edit)
        library_layout.addWidget(browse_files_button)

        # Add widgets
        font = QFont()
        font.setPointSize(16)
        introText = QLabel("Phage Display Analysis Suite")
        introText.setFixedSize(400, 50)
        introText.setFont(font)
        topLayout.addWidget(introText)
        bottom_layout.addRow("Input folder:", folder_layout)
        bottom_layout.addRow("Library file:", library_layout)
        bottom_layout.addRow("Library choice", self.library_choice_combo)
        self._createDisplay(bottom_layout)
        bottom_layout.addRow(button_layout)

        # Organize layouts
        mainLayout.addLayout(topLayout)
        mainLayout.addLayout(bottom_layout)

        # Create a central widget and set the main layout
        central_widget = QWidget()
        central_widget.setLayout(mainLayout)
        self.setCentralWidget(central_widget)

    def _createDisplay(self, layout):
        """
        Create and add the display text edit.
        """
        self.display = QTextEdit()
        self.display.setReadOnly(True)
        self.display.setPlainText(
            "Library data will appear here when selected.")
        layout.addWidget(self.display)

    def browse_folders(self, line_edit):
        """
        Function for browsing and selecting folders.
        """
        home_dir = os.path.expanduser("~")
        dir_path = QFileDialog.getExistingDirectory(
            self,
            "Select Directory",
            home_dir,
            QFileDialog.Option.ShowDirsOnly
        )
        if dir_path:
            line_edit.setText(dir_path)

    def browse_files(self, line_edit):
        """
        Function for browsing and selecting a single file.
        """
        home_dir = os.path.expanduser("~")
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select File",
            home_dir,
            "JSON files (*.json);;All files (*.*)"
        )
        if file_path:
            line_edit.setText(file_path)
            self.populate_combo_box(file_path)

    def populate_combo_box(self, file_path):
        """
        Populate the combo box with the names from the JSON file.
        """
        self.library_choice_combo.clear()
        try:
            with open(file_path, 'r') as file:
                self.json_data = json.load(file)
                libraries = self.json_data.get("libraries", [])
                names = [lib["name"] for lib in libraries if "name" in lib]
                self.library_choice_combo.addItems(names)
        except Exception as e:
            print(f"Error reading or parsing file: {e}")

    def update_display(self):
        """
        Update the display with all the information of the selected library.
        """
        selected_index = self.library_choice_combo.currentIndex()
        if selected_index >= 0 and self.json_data:
            library = self.json_data["libraries"][selected_index]
            display_text = self.format_library_info(library)
            self.display.setPlainText(display_text)
        else:
            self.display.setPlainText(
                "Library data will appear here when selected.")

    def format_library_info(self, library):
        """
        Format all library information for display.
        """
        info = (
            f"Name: {library.get('name', 'N/A')}\n"
            f"Reference: {library.get('reference', 'N/A')}\n"
            f"Trim Motif: {library.get('trim motif', 'N/A')}\n"
            f"Trim Length: {library.get('trim length', 'N/A')}\n"
            f"DNA Reference: {library.get('dna reference', 'N/A')}\n"
            f"Protein Reference: {library.get('protein reference', 'N/A')}\n"
            f"Diversified Residues:\n"
        )
        diversified_residues = library.get("diversified_residues", {})
        for region, positions in diversified_residues.items():
            positions_str = ", ".join(map(str, positions))
            info += f"  {region} - {positions_str}\n"
        return info

    def run_button_click(self):
        """
        Handle the Run button click event.
        """
        selected_index = self.library_choice_combo.currentIndex()
        folder_path = self.folder_path_edit.text()
        library_file_path = self.library_file_edit.text()
        if selected_index >= 0 and self.json_data and folder_path and library_file_path:
            self.display.setPlainText("Running analysis...")
            self.run_process()
            self.display.setPlainText("Analysis complete.")
        else:
            self.display.setPlainText(
                "Please select a folder and a library before running.")

    def run_process(self):
        """
        Perform analysis on the selected data.
        """
        input_folder = self.folder_path_edit.text()
        lib_file_path = self.library_file_edit.text()
        selected_library_name = self.library_choice_combo.currentText()
        if not input_folder or not lib_file_path or not selected_library_name:
            self.display.setPlainText(
                "Please fill in all required fields and select a library.")
            return
        analysis = PdasSeqAnalysis(
            input_path=input_folder,
            lib_file=lib_file_path,
            lib_name=selected_library_name
        )
        analysis.run_analysis(analysis.master_path)


def main():
    """Main function."""
    pdasApp = QApplication([])
    pdasApp.setStyle('fusion')
    pdasWindow = PdasWindow()
    pdasWindow.show()
    sys.exit(pdasApp.exec())


if __name__ == "__main__":
    main()
