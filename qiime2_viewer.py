#!/usr/bin/env python3

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                           QHBoxLayout, QPushButton, QFileDialog, QLabel, 
                           QTabWidget, QTableWidget, QTableWidgetItem)
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

class QIIME2Viewer(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('QIIME2 Results Viewer')
        self.setGeometry(100, 100, 1200, 800)
        
        # Create main widget and layout
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        layout = QVBoxLayout(main_widget)
        
        # Create toolbar
        toolbar = QHBoxLayout()
        layout.addLayout(toolbar)
        
        # Add load button
        load_btn = QPushButton('Load Results Directory')
        load_btn.clicked.connect(self.load_results)
        toolbar.addWidget(load_btn)
        
        # Add status label
        self.status_label = QLabel('No data loaded')
        toolbar.addWidget(self.status_label)
        
        # Create tab widget
        self.tabs = QTabWidget()
        layout.addWidget(self.tabs)
        
        # Initialize data storage
        self.data_dir = None
        
    def load_results(self):
        self.data_dir = QFileDialog.getExistingDirectory(self, 'Select Results Directory')
        if self.data_dir:
            self.status_label.setText(f'Loading data from: {self.data_dir}')
            self.load_all_data()
            
    def load_all_data(self):
        # Clear existing tabs
        while self.tabs.count():
            self.tabs.removeTab(0)
            
        # Load feature table
        feature_table_path = os.path.join(self.data_dir, 'feature_table', 'feature-table.tsv')
        if os.path.exists(feature_table_path):
            self.add_table_tab('Feature Table', feature_table_path)
            
        # Load taxonomy
        taxonomy_path = os.path.join(self.data_dir, 'taxonomy', 'taxonomy.tsv')
        if os.path.exists(taxonomy_path):
            self.add_table_tab('Taxonomy', taxonomy_path)
            
        # Load alpha diversity metrics
        diversity_dir = os.path.join(self.data_dir, 'diversity')
        if os.path.exists(diversity_dir):
            self.add_diversity_tab('Alpha Diversity', diversity_dir)
            
        # Load taxonomic level tables
        taxonomy_dir = os.path.join(self.data_dir, 'taxonomy')
        if os.path.exists(taxonomy_dir):
            for level in range(2, 8):
                level_path = os.path.join(taxonomy_dir, f'level{level}', 'feature-table.tsv')
                if os.path.exists(level_path):
                    self.add_table_tab(f'Level {level} Taxonomy', level_path)
                    
        self.status_label.setText('Data loaded successfully')
        
    def add_table_tab(self, name, file_path):
        try:
            # Read data
            df = pd.read_csv(file_path, sep='\t')
            
            # Create tab widget
            tab = QWidget()
            layout = QVBoxLayout(tab)
            
            # Create table widget
            table = QTableWidget()
            layout.addWidget(table)
            
            # Set table dimensions
            table.setRowCount(len(df))
            table.setColumnCount(len(df.columns))
            
            # Set headers
            table.setHorizontalHeaderLabels(df.columns)
            
            # Fill data
            for i in range(len(df)):
                for j in range(len(df.columns)):
                    item = QTableWidgetItem(str(df.iloc[i, j]))
                    table.setItem(i, j, item)
            
            # Add tab
            self.tabs.addTab(tab, name)
            
        except Exception as e:
            print(f"Error loading {name}: {str(e)}")
            
    def add_diversity_tab(self, name, directory):
        try:
            # Create tab widget
            tab = QWidget()
            layout = QVBoxLayout(tab)
            
            # Create figure
            fig, ax = plt.subplots(figsize=(10, 6))
            canvas = FigureCanvas(fig)
            layout.addWidget(canvas)
            
            # Plot alpha diversity metrics
            metrics = ['observed_features', 'shannon', 'evenness', 'faith_pd']
            data = []
            labels = []
            
            for metric in metrics:
                metric_path = os.path.join(directory, metric, 'alpha-diversity.tsv')
                if os.path.exists(metric_path):
                    df = pd.read_csv(metric_path, sep='\t')
                    data.append(df.iloc[:, 1].values)
                    labels.append(metric)
            
            if data:
                ax.boxplot(data, labels=labels)
                ax.set_title('Alpha Diversity Metrics')
                ax.set_ylabel('Value')
                
                # Add tab
                self.tabs.addTab(tab, name)
                
        except Exception as e:
            print(f"Error loading diversity tab: {str(e)}")

def main():
    app = QApplication(sys.argv)
    viewer = QIIME2Viewer()
    viewer.show()
    sys.exit(app.exec_()) 