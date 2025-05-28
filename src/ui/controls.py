import tkinter as tk
from tkinter import ttk
from typing import Callable, List, Tuple

class ControlPanel:
    """Control panel for simulation parameters"""
    
    def __init__(self, parent: ttk.Frame, title: str):
        self.frame = ttk.LabelFrame(parent, text=title, padding="10")
        self.frame.pack(fill=tk.X, pady=10)
        
        self.controls: List[Tuple[str, tk.Widget]] = []
    
    def add_button(self, text: str, command: Callable, tooltip: str = "") -> ttk.Button:
        """Add a button to the control panel"""
        btn = ttk.Button(self.frame, text=text, command=command)
        btn.pack(fill=tk.X, pady=5)
        if tooltip:
            self._create_tooltip(btn, tooltip)
        self.controls.append((text, btn))
        return btn
    
    def add_slider(self, text: str, from_: float, to: float, 
                  command: Callable, tooltip: str = "") -> ttk.Scale:
        """Add a slider to the control panel"""
        frame = ttk.Frame(self.frame)
        frame.pack(fill=tk.X, pady=5)
        
        label = ttk.Label(frame, text=text)
        label.pack(side=tk.LEFT)
        
        slider = ttk.Scale(frame, from_=from_, to=to, orient=tk.HORIZONTAL, command=command)
        slider.pack(side=tk.RIGHT, fill=tk.X, expand=True)
        
        if tooltip:
            self._create_tooltip(slider, tooltip)
        
        self.controls.append((text, slider))
        return slider
    
    def _create_tooltip(self, widget: tk.Widget, text: str) -> None:
        """Create a tooltip for a widget"""
        def show_tooltip(event):
            tooltip = tk.Toplevel()
            tooltip.wm_overrideredirect(True)
            tooltip.wm_geometry(f"+{event.x_root+10}+{event.y_root+10}")
            
            label = ttk.Label(tooltip, text=text, background="#ffffe0", 
                            relief=tk.SOLID, borderwidth=1)
            label.pack()
            
            def hide_tooltip():
                tooltip.destroy()
            
            widget.tooltip = tooltip
            widget.bind("<Leave>", lambda e: hide_tooltip())
        
        widget.bind("<Enter>", show_tooltip)
