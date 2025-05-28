import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
import matplotlib.backends.backend_tkagg as backend_tkagg
from typing import Optional, Tuple

class VisualizationPanel:
    """Panel for displaying quantum simulation visualizations"""
    
    def __init__(self, parent: ttk.Frame):
        self.frame = ttk.Frame(parent)
        self.frame.pack(fill=tk.BOTH, expand=True)
        
        self.fig = plt.Figure(figsize=(8, 6), dpi=100)
        self.canvas = backend_tkagg.FigureCanvasTkAgg(self.fig, master=self.frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        self.toolbar = backend_tkagg.NavigationToolbar2Tk(self.canvas, self.frame)
        self.toolbar.update()
        
        self.ax = self.fig.add_subplot(111)
        self.setup_plot()
    
    def setup_plot(self, title: str = "Quantum Wave Packet Simulation",
                  xlabel: str = "Position", ylabel: str = "Probability Density") -> None:
        """Set up the basic plot configuration"""
        self.ax.set_title(title)
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        self.ax.grid(True)
        self.canvas.draw()
    
    def update_plot(self, x: list, y: list, clear: bool = True) -> None:
        """Update the plot with new data"""
        if clear:
            self.ax.clear()
            self.setup_plot()
        self.ax.plot(x, y)
        self.canvas.draw()
    
    def clear_plot(self) -> None:
        """Clear the current plot"""
        self.ax.clear()
        self.setup_plot()
