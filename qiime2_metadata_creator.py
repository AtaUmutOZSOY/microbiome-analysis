#!/usr/bin/env python3

import sys
import os
import pandas as pd
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                           QHBoxLayout, QPushButton, QFileDialog, QLabel, 
                           QTableWidget, QTableWidgetItem, QMessageBox,
                           QInputDialog, QComboBox, QHeaderView)
from PyQt5.QtCore import Qt

class QIIME2MetadataCreator(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('QIIME2 Metadata Creator')
        self.setGeometry(100, 100, 1200, 800)
        
        # Create main widget and layout
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        layout = QVBoxLayout(main_widget)
        
        # Create toolbar
        toolbar = QHBoxLayout()
        layout.addLayout(toolbar)
        
        # Add buttons
        new_btn = QPushButton('New Metadata')
        new_btn.clicked.connect(self.create_new_metadata)
        toolbar.addWidget(new_btn)
        
        load_btn = QPushButton('Load Existing Metadata')
        load_btn.clicked.connect(self.load_metadata)
        toolbar.addWidget(load_btn)
        
        add_col_btn = QPushButton('Add Column')
        add_col_btn.clicked.connect(self.add_column)
        toolbar.addWidget(add_col_btn)
        
        add_row_btn = QPushButton('Add Row')
        add_row_btn.clicked.connect(self.add_row)
        toolbar.addWidget(add_row_btn)
        
        save_btn = QPushButton('Save Metadata')
        save_btn.clicked.connect(self.save_metadata)
        toolbar.addWidget(save_btn)
        
        validate_btn = QPushButton('Validate Metadata')
        validate_btn.clicked.connect(self.validate_metadata)
        toolbar.addWidget(validate_btn)
        
        # Add status label
        self.status_label = QLabel('No metadata loaded')
        toolbar.addWidget(self.status_label)
        
        # Create table widget
        self.table = QTableWidget()
        layout.addWidget(self.table)
        
        # Initialize with basic columns
        self.initialize_table()
        
    def initialize_table(self):
        # Basic required columns for QIIME2
        columns = ['#SampleID', 'sample-type', 'group']
        self.table.setColumnCount(len(columns))
        self.table.setHorizontalHeaderLabels(columns)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        
        # Add one empty row
        self.table.setRowCount(1)
        for col in range(len(columns)):
            self.table.setItem(0, col, QTableWidgetItem(''))
            
    def create_new_metadata(self):
        reply = QMessageBox.question(self, 'Confirm', 
                                   'This will clear current data. Continue?',
                                   QMessageBox.Yes | QMessageBox.No)
        
        if reply == QMessageBox.Yes:
            self.table.clear()
            self.initialize_table()
            self.status_label.setText('New metadata created')
            
    def load_metadata(self):
        file_path, _ = QFileDialog.getOpenFileName(self, 'Open Metadata File', '', 
                                                 'TSV files (*.tsv);;All Files (*)')
        if file_path:
            try:
                df = pd.read_csv(file_path, sep='\t')
                self.table.setRowCount(len(df))
                self.table.setColumnCount(len(df.columns))
                self.table.setHorizontalHeaderLabels(df.columns)
                
                for i in range(len(df)):
                    for j in range(len(df.columns)):
                        self.table.setItem(i, j, QTableWidgetItem(str(df.iloc[i, j])))
                        
                self.status_label.setText(f'Loaded: {file_path}')
            except Exception as e:
                QMessageBox.critical(self, 'Error', f'Error loading file: {str(e)}')
                
    def add_column(self):
        column_types = ['categorical', 'numeric', 'datetime']
        column_name, ok = QInputDialog.getText(self, 'Add Column', 
                                             'Enter column name:')
        if ok and column_name:
            type_dialog = QInputDialog()
            type_dialog.setComboBoxItems(column_types)
            type_dialog.setWindowTitle('Column Type')
            type_dialog.setLabelText('Select column type:')
            if type_dialog.exec_() == QInputDialog.Accepted:
                col_type = type_dialog.textValue()
                
                # Add new column
                current_cols = self.table.columnCount()
                self.table.setColumnCount(current_cols + 1)
                self.table.setHorizontalHeaderItem(current_cols, 
                                                 QTableWidgetItem(f'{column_name} [{col_type}]'))
                
                # Fill with empty items
                for row in range(self.table.rowCount()):
                    self.table.setItem(row, current_cols, QTableWidgetItem(''))
                    
    def add_row(self):
        current_rows = self.table.rowCount()
        self.table.setRowCount(current_rows + 1)
        
        # Fill with empty items
        for col in range(self.table.columnCount()):
            self.table.setItem(current_rows, col, QTableWidgetItem(''))
            
    def save_metadata(self):
        file_path, _ = QFileDialog.getSaveFileName(self, 'Save Metadata File', '', 
                                                 'TSV files (*.tsv);;All Files (*)')
        if file_path:
            try:
                # Convert table to dataframe
                rows = self.table.rowCount()
                cols = self.table.columnCount()
                headers = [self.table.horizontalHeaderItem(i).text() for i in range(cols)]
                data = []
                
                for row in range(rows):
                    row_data = []
                    for col in range(cols):
                        item = self.table.item(row, col)
                        row_data.append('' if item is None else item.text())
                    data.append(row_data)
                    
                df = pd.DataFrame(data, columns=headers)
                df.to_csv(file_path, sep='\t', index=False)
                self.status_label.setText(f'Saved: {file_path}')
                
            except Exception as e:
                QMessageBox.critical(self, 'Error', f'Error saving file: {str(e)}')
                
    def validate_metadata(self):
        # Basic validation rules
        errors = []
        
        # Check if #SampleID column exists
        if '#SampleID' not in [self.table.horizontalHeaderItem(i).text() 
                              for i in range(self.table.columnCount())]:
            errors.append('Missing #SampleID column')
            
        # Check for empty cells
        for row in range(self.table.rowCount()):
            for col in range(self.table.columnCount()):
                item = self.table.item(row, col)
                if item is None or not item.text().strip():
                    errors.append(f'Empty cell at row {row+1}, column {col+1}')
                    
        # Check for duplicate sample IDs
        sample_ids = []
        sample_id_col = [self.table.horizontalHeaderItem(i).text() 
                        for i in range(self.table.columnCount())].index('#SampleID')
        
        for row in range(self.table.rowCount()):
            item = self.table.item(row, sample_id_col)
            if item is not None:
                sample_id = item.text()
                if sample_id in sample_ids:
                    errors.append(f'Duplicate sample ID: {sample_id}')
                sample_ids.append(sample_id)
                
        # Show validation results
        if errors:
            QMessageBox.warning(self, 'Validation Errors', 
                              'Found the following errors:\n' + '\n'.join(errors))
        else:
            QMessageBox.information(self, 'Validation Success', 
                                  'Metadata is valid for QIIME2!')

def main():
    app = QApplication(sys.argv)
    creator = QIIME2MetadataCreator()
    creator.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main() 