#!/usr/bin/env python3

"""
gui_pdas.py

Simple GUI for seq_analysis.py.
"""


import sys
import os
import json
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QWidget, QHBoxLayout,
    QLineEdit, QFormLayout, QTextEdit, QFileDialog, QPushButton,
    QComboBox, QTabWidget
    )
from PyQt6.QtGui import QIcon
from seq_analysis import PdasSeqAnalysis


class PdasWindow(QMainWindow):
    """Main window."""

    def __init__(self):
        super().__init__()
        self.setWindowIcon(QIcon('C:\\Users\\spyro\\Documents\\work\\programming\\phage_display_analysis_suite_pdas\\files\\images\\icon.ico'))
        self.setWindowTitle("Phage Display Analysis Suite")
        self.resize(200, 450)
        self.json_data = None

        # Create the tab widget tabs
        self.tabs = QTabWidget()
        self.tab1 = QWidget()
        self.tab2 = QWidget()
        self._setup_tab1()
        self._setup_tab2()
        self.tabs.addTab(self.tab1, "Main")
        self.tabs.addTab(self.tab2, "Troubleshooting")
        self.setCentralWidget(self.tabs)

    def _setup_tab1(self):
        """Sets up content for tab 1 (main)."""
        # Define layouts
        mainLayout = QVBoxLayout()
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

        # Add widgets to layouts
        bottom_layout.addRow("Input folder:", folder_layout)
        bottom_layout.addRow("Library file:", library_layout)
        bottom_layout.addRow("Library choice", self.library_choice_combo)
        self._createDisplay(bottom_layout)
        bottom_layout.addRow(button_layout)

        # Set layouts to tabs
        mainLayout.addLayout(bottom_layout)
        self.tab1.setLayout(mainLayout)

    def _setup_tab2(self):
        """Sets up content for tab 2 (troubleshooting)."""
        layout = QVBoxLayout()

        # Create a QTextEdit widget and set HTML content
        html_text = QTextEdit()
        html_text.setHtml("""
        <h2>Troubleshooting</h2>
        <p>This will be updated in the future with common issues.</p>
        <p>To report any bugs, please do so on GitHub or <a href='mailto:spyrolivia@gmail.com'>email me</a>.</p>
        <p>- Liv</p>
        """)

        # Make sure the QTextEdit widget is not editable
        html_text.setReadOnly(True)

        # Set the minimum height and width to make sure it is visible
        html_text.setMinimumHeight(200)
        html_text.setMinimumWidth(400)

        # Add the QTextEdit to the layout
        layout.addWidget(html_text)
        self.tab2.setLayout(layout)

    def _createDisplay(self, layout):
        """
        Create and add the display text edit.
        """
        self.display = QTextEdit()
        self.display.setReadOnly(True)
        self.display.setHtml(
            "<p><i>Library data will appear here when selected.</i></p>")
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
            self.display.setHtml(display_text)
        else:
            self.display.setHtml(
                "<p><i>Library data will appear here when selected.</i></p>")

    def format_library_info(self, library):
        """
        Format all library information for display.
        """
        info = (
            f"<p><b>Name:</b> {library.get('name', 'N/A')}</p>"
            f"<p><b>Reference:</b> {library.get('reference', 'N/A')}</p>"
            f"<p><b>Trim Motif:</b> {library.get('trim motif', 'N/A')}</p>"
            f"<p><b>Trim Length:</b> {library.get('trim length', 'N/A')}</p>"
            f"<p><b>DNA Reference:</b> {library.get('dna reference', 'N/A')}</p>"
            f"<p><b>Protein Reference:</b> {library.get('protein reference', 'N/A')}</p>"
            f"<p><b>Diversified Residues:</b><br>"
        )
        diversified_residues = library.get("diversified_residues", {})
        for region, positions in diversified_residues.items():
            positions_str = ", ".join(map(str, positions))
            info += f"<i>{region}</i> - {positions_str}<br>"
        info += "</p>"
        return info

    def run_button_click(self):
        """
        Handle the Run button click event.
        """
        selected_index = self.library_choice_combo.currentIndex()
        folder_path = self.folder_path_edit.text()
        library_file_path = self.library_file_edit.text()
        if selected_index >= 0 and self.json_data and folder_path and library_file_path:
            self.display.setHtml("<p><i>Running analysis...</i></p>")
            self.run_process()
            self.display.setHtml("<p><i>Analysis complete.</i></p>")
        else:
            self.display.setHtml(
                "<p><i>Please select a folder and a library before running.</i></p>")

    def run_process(self):
        """
        Perform analysis on the selected data.
        """
        input_folder = self.folder_path_edit.text()
        lib_file_path = self.library_file_edit.text()
        selected_library_name = self.library_choice_combo.currentText()
        if not input_folder or not lib_file_path or not selected_library_name:
            self.display.setHtml(
                "<p><i>Please fill in all required fields and select a library.</i></p>")
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
