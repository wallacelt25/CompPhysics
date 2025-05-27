import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_tkagg as backend_tkagg
import threading
import time
from matplotlib import style

from .simulation.quantum_core import QuantumWavePacket, QuantumWavePacket2D
from .simulation.functions.free_particle import FreqParticleFunction
from .simulation.functions.barrier_tunnel import BarrierTunnelingFunction
from .simulation.functions.double_slit import DoubleSlitFunction
from .simulation.functions.interactive_observables import InteractiveObservablesFunction
from .simulation.functions.tunneling import TunnelingFunction

# Set dark theme for matplotlib
style.use('dark_background')

class QuantumApp:
    """
    Main application class for the Quantum Wave Packet Simulation GUI
    Implements a modern dark theme with purple accents and optimized performance
    """
    def __init__(self, master):
        # Configure the main window
        self.master = master
        self.master.title("Quantum Wave Packet Simulation")
        self.master.geometry("1200x800")
        
        # Define color scheme
        self.colors = {
            'bg': '#1a1a1a',  # Dark background
            'fg': '#ffffff',  # White text
            'accent': '#9b4dca',  # Purple accent
            'button': '#2d2d2d',  # Dark button background
            'button_active': '#3d3d3d',  # Active button state
            'plot_bg': '#000000',  # Black plot background
            'grid': '#333333'  # Dark grid lines
        }
        
        # Configure ttk styles
        self.configure_styles()
        
        # Set up the main container with dark background
        self.main_frame = ttk.Frame(self.master, padding="10", style='Dark.TFrame')
        self.main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Create a notebook (tabbed interface) for 1D and 2D
        self.notebook = ttk.Notebook(self.main_frame, style='Dark.TNotebook')
        self.notebook.pack(fill=tk.BOTH, expand=True)
        
        # Create tabs for 1D and 2D
        self.tab_1d = ttk.Frame(self.notebook, style='Dark.TFrame')
        self.tab_2d = ttk.Frame(self.notebook, style='Dark.TFrame')
        
        self.notebook.add(self.tab_1d, text="1D Simulation")
        self.notebook.add(self.tab_2d, text="2D Simulation")
        
        # Create interfaces for both tabs
        self.create_1d_interface()
        self.create_2d_interface()
        
        # Initialize simulation objects with optimized parameters
        self.sim_1d = QuantumWavePacket()
        self.sim_1d.initialize_gaussian_wavepacket()
        
        self.sim_2d = QuantumWavePacket2D()
        self.sim_2d.initialize_gaussian_wavepacket()
        
        # Initialize function modules
        self.initialize_function_modules()
        
        # Animation control variables (separate for 1D and 2D)
        self.animation_running_1d = False
        self.animation_thread_1d = None
        self.animation_running_2d = False
        self.animation_thread_2d = None
        self.frame_time = 1/30  # Target 30 FPS for smooth animation
        
        # Current active function
        self.current_function_1d = self.free_particle_1d
        self.current_function_2d = self.free_particle_2d
        print('[DEBUG] Setting default current_function_1d and 2d and calling setup()')
        self.current_function_1d.setup()
        self.current_function_2d.setup()

    def configure_styles(self):
        """Configure ttk styles for the dark theme"""
        style = ttk.Style()
        
        # Configure the dark theme
        style.configure('Dark.TFrame', background=self.colors['bg'])
        style.configure('Dark.TLabel', background=self.colors['bg'], foreground=self.colors['fg'])
        style.configure('Dark.TButton', 
                       background=self.colors['button'],
                       foreground=self.colors['fg'],
                       font=("Arial", 10, "bold"),
                       padding=5)
        style.map('Dark.TButton',
                 background=[('active', self.colors['accent'])],
                 foreground=[('active', self.colors['fg'])])
        style.configure('Dark.TNotebook', background=self.colors['bg'])
        style.configure('Dark.TLabelframe', 
                       background=self.colors['bg'],
                       foreground=self.colors['fg'])
        style.configure('Dark.TLabelframe.Label', 
                       background=self.colors['bg'],
                       foreground=self.colors['fg'])

    def initialize_function_modules(self):
        """Initialize all function modules for both 1D and 2D simulations"""
        # 1D function modules
        self.free_particle_1d = FreqParticleFunction(self, is_2d=False)
        self.barrier_tunnel_1d = BarrierTunnelingFunction(self, is_2d=False)
        self.double_slit_1d = DoubleSlitFunction(self, is_2d=False)
        self.interactive_obs_1d = InteractiveObservablesFunction(self, is_2d=False)
        self.tunneling_1d = TunnelingFunction(self, is_2d=False)

        # 2D function modules
        self.free_particle_2d = FreqParticleFunction(self, is_2d=True)
        self.barrier_tunnel_2d = BarrierTunnelingFunction(self, is_2d=True)
        self.double_slit_2d = DoubleSlitFunction(self, is_2d=True)
        self.interactive_obs_2d = InteractiveObservablesFunction(self, is_2d=True)
        self.tunneling_2d = TunnelingFunction(self, is_2d=True)

    def create_1d_interface(self):
        """Create the interface for the 1D tab with dark theme"""
        # Create sidebar frame with dark theme
        self.sidebar_1d = ttk.Frame(self.tab_1d, padding="10", style='Dark.TFrame')
        self.sidebar_1d.pack(side=tk.LEFT, fill=tk.Y)
        
        # Create visualization frame
        self.viz_frame_1d = ttk.Frame(self.tab_1d, style='Dark.TFrame')
        self.viz_frame_1d.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        
        # Create function buttons in sidebar with dark theme
        ttk.Label(self.sidebar_1d, 
                 text="Simulation Functions", 
                 font=("Arial", 12, "bold"),
                 style='Dark.TLabel').pack(pady=10)
        
        # Create buttons for each function with dark theme
        functions = [
            ("Free Particle", lambda: self.switch_function_1d(self.free_particle_1d)),
            ("Barrier Tunneling", lambda: self.switch_function_1d(self.barrier_tunnel_1d)),
            ("Double Slit", lambda: self.switch_function_1d(self.double_slit_1d)),
            ("Interactive Observables", lambda: self.switch_function_1d(self.interactive_obs_1d)),
            ("Tunneling Analysis", lambda: self.switch_function_1d(self.tunneling_1d))
        ]
        
        self.function_buttons_1d = []
        for text, command in functions:
            btn = tk.Button(self.sidebar_1d, text=text, command=command, width=20,
                            bg=self.colors['button'], fg=self.colors['fg'],
                            activebackground=self.colors['accent'], activeforeground=self.colors['fg'],
                            font=("Arial", 10, "bold"),
                            relief=tk.FLAT, highlightthickness=0, bd=0)
            btn.pack(pady=5, padx=10, fill=tk.X)
            self.function_buttons_1d.append(btn)
            
        # Separator
        ttk.Separator(self.sidebar_1d, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        
        # Controls frame with dark theme
        self.controls_frame_1d = ttk.LabelFrame(self.sidebar_1d, 
                                              text="Simulation Controls", 
                                              padding="10",
                                              style='Dark.TLabelframe')
        self.controls_frame_1d.pack(fill=tk.X, pady=10)
        
        # Add control buttons with dark theme
        self.btn_start_1d = tk.Button(self.controls_frame_1d, text="Start Simulation", 
                                      command=self.start_simulation_1d,
                                      bg=self.colors['button'], fg=self.colors['fg'],
                                      activebackground=self.colors['accent'], activeforeground=self.colors['fg'],
                                      font=("Arial", 10, "bold"),
                                      relief=tk.FLAT, highlightthickness=0, bd=0)
        self.btn_start_1d.pack(fill=tk.X, pady=5)
        self.btn_stop_1d = tk.Button(self.controls_frame_1d, text="Stop Simulation", 
                                     command=self.stop_simulation_1d,
                                     bg=self.colors['button'], fg=self.colors['fg'],
                                     activebackground=self.colors['accent'], activeforeground=self.colors['fg'],
                                     font=("Arial", 10, "bold"),
                                     relief=tk.FLAT, highlightthickness=0, bd=0)
        self.btn_stop_1d.pack(fill=tk.X, pady=5)
        self.btn_reset_1d = tk.Button(self.controls_frame_1d, text="Reset Simulation", 
                                      command=self.reset_simulation_1d,
                                      bg=self.colors['button'], fg=self.colors['fg'],
                                      activebackground=self.colors['accent'], activeforeground=self.colors['fg'],
                                      font=("Arial", 10, "bold"),
                                      relief=tk.FLAT, highlightthickness=0, bd=0)
        self.btn_reset_1d.pack(fill=tk.X, pady=5)
        
        # Create a frame for function-specific controls
        self.function_controls_frame_1d = ttk.Frame(self.sidebar_1d, style='Dark.TFrame')
        self.function_controls_frame_1d.pack(fill=tk.X, pady=10, expand=True)
        
        # Create matplotlib figure for 1D visualization with dark theme
        self.fig_1d = plt.Figure(figsize=(8, 6), dpi=100, facecolor=self.colors['plot_bg'])
        self.canvas_1d = backend_tkagg.FigureCanvasTkAgg(self.fig_1d, master=self.viz_frame_1d)
        self.canvas_1d.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Add toolbar with dark theme
        self.toolbar_1d = backend_tkagg.NavigationToolbar2Tk(self.canvas_1d, self.viz_frame_1d)
        self.toolbar_1d.update()
        
        # Initialize with an empty plot with dark theme
        self.ax_1d = self.fig_1d.add_subplot(111)
        self.ax_1d.set_title("Quantum Wave Packet Simulation (1D)", color=self.colors['fg'])
        self.ax_1d.set_xlabel("Position", color=self.colors['fg'])
        self.ax_1d.set_ylabel("Probability Density", color=self.colors['fg'])
        self.ax_1d.grid(True, color=self.colors['grid'])
        self.ax_1d.set_facecolor(self.colors['plot_bg'])
        self.canvas_1d.draw()
    
    def create_2d_interface(self):
        """Create the interface for the 2D tab with dark theme"""
        # Create sidebar frame with dark theme
        self.sidebar_2d = ttk.Frame(self.tab_2d, padding="10", style='Dark.TFrame')
        self.sidebar_2d.pack(side=tk.LEFT, fill=tk.Y)
        
        # Create visualization frame
        self.viz_frame_2d = ttk.Frame(self.tab_2d, style='Dark.TFrame')
        self.viz_frame_2d.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        
        # Create function buttons in sidebar with dark theme
        ttk.Label(self.sidebar_2d, 
                 text="Simulation Functions", 
                 font=("Arial", 12, "bold"),
                 style='Dark.TLabel').pack(pady=10)
        
        # Create buttons for each function with dark theme
        functions = [
            ("Free Particle", lambda: self.switch_function_2d(self.free_particle_2d)),
            ("Barrier Tunneling", lambda: self.switch_function_2d(self.barrier_tunnel_2d)),
            ("Double Slit", lambda: self.switch_function_2d(self.double_slit_2d)),
            ("Interactive Observables", lambda: self.switch_function_2d(self.interactive_obs_2d)),
            ("Tunneling Analysis", lambda: self.switch_function_2d(self.tunneling_2d))
        ]
        
        self.function_buttons_2d = []
        for text, command in functions:
            btn = tk.Button(self.sidebar_2d, text=text, command=command, width=20,
                            bg=self.colors['button'], fg=self.colors['fg'],
                            activebackground=self.colors['accent'], activeforeground=self.colors['fg'],
                            font=("Arial", 10, "bold"),
                            relief=tk.FLAT, highlightthickness=0, bd=0)
            btn.pack(pady=5, padx=10, fill=tk.X)
            self.function_buttons_2d.append(btn)
            
        # Separator
        ttk.Separator(self.sidebar_2d, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        
        # Controls frame with dark theme
        self.controls_frame_2d = ttk.LabelFrame(self.sidebar_2d, 
                                              text="Simulation Controls", 
                                              padding="10",
                                              style='Dark.TLabelframe')
        self.controls_frame_2d.pack(fill=tk.X, pady=10)
        
        # Add control buttons with dark theme
        self.btn_start_2d = tk.Button(self.controls_frame_2d, text="Start Simulation", 
                                      command=self.start_simulation_2d,
                                      bg=self.colors['button'], fg=self.colors['fg'],
                                      activebackground=self.colors['accent'], activeforeground=self.colors['fg'],
                                      font=("Arial", 10, "bold"),
                                      relief=tk.FLAT, highlightthickness=0, bd=0)
        self.btn_start_2d.pack(fill=tk.X, pady=5)
        self.btn_stop_2d = tk.Button(self.controls_frame_2d, text="Stop Simulation", 
                                     command=self.stop_simulation_2d,
                                     bg=self.colors['button'], fg=self.colors['fg'],
                                     activebackground=self.colors['accent'], activeforeground=self.colors['fg'],
                                     font=("Arial", 10, "bold"),
                                     relief=tk.FLAT, highlightthickness=0, bd=0)
        self.btn_stop_2d.pack(fill=tk.X, pady=5)
        self.btn_reset_2d = tk.Button(self.controls_frame_2d, text="Reset Simulation", 
                                      command=self.reset_simulation_2d,
                                      bg=self.colors['button'], fg=self.colors['fg'],
                                      activebackground=self.colors['accent'], activeforeground=self.colors['fg'],
                                      font=("Arial", 10, "bold"),
                                      relief=tk.FLAT, highlightthickness=0, bd=0)
        self.btn_reset_2d.pack(fill=tk.X, pady=5)
        
        # Create a frame for function-specific controls
        self.function_controls_frame_2d = ttk.Frame(self.sidebar_2d, style='Dark.TFrame')
        self.function_controls_frame_2d.pack(fill=tk.X, pady=10, expand=True)
        
        # Create matplotlib figure for 2D visualization with dark theme
        self.fig_2d = plt.Figure(figsize=(8, 6), dpi=100, facecolor=self.colors['plot_bg'])
        self.canvas_2d = backend_tkagg.FigureCanvasTkAgg(self.fig_2d, master=self.viz_frame_2d)
        self.canvas_2d.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Add toolbar with dark theme
        self.toolbar_2d = backend_tkagg.NavigationToolbar2Tk(self.canvas_2d, self.viz_frame_2d)
        self.toolbar_2d.update()
        
        # Initialize with an empty 2D plot with dark theme
        self.ax_2d = self.fig_2d.add_subplot(111)
        self.ax_2d.set_title("Quantum Wave Packet Simulation (2D)", color=self.colors['fg'])
        self.ax_2d.set_xlabel("X Position", color=self.colors['fg'])
        self.ax_2d.set_ylabel("Y Position", color=self.colors['fg'])
        self.ax_2d.grid(True, color=self.colors['grid'])
        self.ax_2d.set_facecolor(self.colors['plot_bg'])
        self.canvas_2d.draw()
    
    def switch_function_1d(self, function):
        """Switch to a different 1D function"""
        self.stop_simulation_1d()
        
        # Clear function-specific controls
        for widget in self.function_controls_frame_1d.winfo_children():
            widget.destroy()
        
        self.current_function_1d = function
        function.setup()
    
    def switch_function_2d(self, function):
        """Switch to a different 2D function"""
        self.stop_simulation_2d()
        
        # Clear function-specific controls
        for widget in self.function_controls_frame_2d.winfo_children():
            widget.destroy()
        
        self.current_function_2d = function
        function.setup()
    
    def start_simulation_1d(self):
        """Start the 1D simulation"""
        if self.current_function_1d and not self.animation_running_1d:
            self.animation_running_1d = True
            self.animation_thread_1d = threading.Thread(target=self.run_animation_1d)
            self.animation_thread_1d.daemon = True
            self.animation_thread_1d.start()
    
    def stop_simulation_1d(self):
        """Stop the 1D simulation"""
        self.animation_running_1d = False
        if self.animation_thread_1d is not None:
            self.animation_thread_1d.join(timeout=0.5)
    
    def reset_simulation_1d(self):
        """Reset the 1D simulation"""
        if self.current_function_1d:
            self.stop_simulation_1d()
            self.current_function_1d.reset()
    
    def start_simulation_2d(self):
        """Start the 2D simulation"""
        if self.current_function_2d and not self.animation_running_2d:
            self.animation_running_2d = True
            self.animation_thread_2d = threading.Thread(target=self.run_animation_2d)
            self.animation_thread_2d.daemon = True
            self.animation_thread_2d.start()
    
    def stop_simulation_2d(self):
        """Stop the 2D simulation"""
        self.animation_running_2d = False
        if self.animation_thread_2d is not None:
            self.animation_thread_2d.join(timeout=0.5)
    
    def reset_simulation_2d(self):
        """Reset the 2D simulation"""
        if self.current_function_2d:
            self.stop_simulation_2d()
            self.current_function_2d.reset()
    
    def run_animation_1d(self):
        """Run the 1D animation with optimized performance"""
        print("[DEBUG] 1D animation thread started")
        while self.animation_running_1d and self.current_function_1d:
            start_time = time.time()
            self.current_function_1d.update()
            print(f"[DEBUG] 1D time: {getattr(self.sim_1d, 'time', 'N/A')}")
            self.canvas_1d.draw()
            self.master.update_idletasks()
            elapsed = time.time() - start_time
            if elapsed < self.frame_time:
                time.sleep(self.frame_time - elapsed)
        print("[DEBUG] 1D animation thread ended")
    
    def run_animation_2d(self):
        """Run the 2D animation with robust update method selection"""
        print("[DEBUG] 2D animation thread started")
        try:
            while self.animation_running_2d and self.current_function_2d:
                start_time = time.time()
                try:
                    if hasattr(self.current_function_2d, 'update'):
                        self.current_function_2d.update()
                    elif hasattr(self.current_function_2d, 'step'):
                        self.current_function_2d.step()
                    elif hasattr(self.current_function_2d, 'update_plot'):
                        self.current_function_2d.update_plot()
                    print(f"[DEBUG] 2D time: {getattr(self.sim_2d, 'time', 'N/A')}")
                    self.canvas_2d.draw()
                    self.master.update_idletasks()
                except Exception as e:
                    import traceback
                    print("[2D Animation Exception]", e)
                    traceback.print_exc()
                    break
                elapsed = time.time() - start_time
                if elapsed < self.frame_time:
                    time.sleep(self.frame_time - elapsed)
        finally:
            print("[DEBUG] 2D animation thread ended")