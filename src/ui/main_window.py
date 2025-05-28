import tkinter as tk
from tkinter import ttk
from typing import Optional, Callable

class MainWindow:
    """Main application window class"""
    
    def __init__(self, title: str = "Quantum Wave Packet Simulation"):
        self.root = tk.Tk()
        self.root.title(title)
        self.root.geometry("1200x800")
        self.root.configure(bg="#f0f0f0")
        
        self.main_frame = ttk.Frame(self.root, padding="10")
        self.main_frame.pack(fill=tk.BOTH, expand=True)
        
        self.notebook = ttk.Notebook(self.main_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True)
        
        self.tab_1d = ttk.Frame(self.notebook)
        self.tab_2d = ttk.Frame(self.notebook)
        
        self.notebook.add(self.tab_1d, text="1D Simulation")
        self.notebook.add(self.tab_2d, text="2D Simulation")
        
        self.status_bar = ttk.Label(self.root, text="Ready", relief=tk.SUNKEN, anchor=tk.W)
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)
    
    def set_status(self, message: str) -> None:
        """Update the status bar message"""
        self.status_bar.config(text=message)
    
    def run(self) -> None:
        """Start the main event loop"""
        self.root.mainloop()
