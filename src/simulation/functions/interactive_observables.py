import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from ..quantum_core import QuantumWavePacket, QuantumWavePacket2D

class InteractiveObservablesFunction:
    """
    Combined interactive mode and observables dashboard for quantum wave packet
    """
    def __init__(self, app, is_2d=False):
        self.app = app
        self.is_2d = is_2d
        self.controls_frame = self.app.function_controls_frame_2d if is_2d else self.app.function_controls_frame_1d
        self.fig = self.app.fig_2d if is_2d else self.app.fig_1d
        self.canvas = self.app.canvas_2d if is_2d else self.app.canvas_1d
        
        # Parameters
        self.position_x_var = tk.DoubleVar(value=-5.0)
        self.position_y_var = tk.DoubleVar(value=0.0) if is_2d else None
        self.momentum_x_var = tk.DoubleVar(value=2.0)
        self.momentum_y_var = tk.DoubleVar(value=0.0) if is_2d else None
        self.width_var = tk.DoubleVar(value=1.0 if is_2d else 0.5)
        self.potential_var = tk.StringVar(value="none")
        self.pot_height_var = tk.DoubleVar(value=1.0)
        self.colormap_var = tk.StringVar(value="viridis") if is_2d else None
        self.visualization_var = tk.StringVar(value="probability")
        self.barrier_width_var = tk.DoubleVar(value=1.0)
        self.barrier_position_var = tk.DoubleVar(value=0.0)
        
        # Interactive mode needs live preview
        self.preview_simulation = QuantumWavePacket2D() if is_2d else QuantumWavePacket()
        
        # For time evolution tracking
        self.position_history = []
        self.momentum_history = []
        self.energy_history = []
        self.time_history = []
        
        self.last_potential_func_1d = None
        self.last_potential_func_2d = None
        
        self.transmission_history = []
        self.reflection_history = []
    
    def setup(self):
        """Set up the combined interactive mode and observables dashboard"""
        # Clear the plot
        self.fig.clear()
        
        # Create the plot layout
        if self.is_2d:
            # 2D layout
            gs = gridspec.GridSpec(2, 2, figure=self.fig)
            
            # Wave packet plot
            self.ax_wave = self.fig.add_subplot(gs[0, 0])
            # Position trajectory plot
            self.ax_position = self.fig.add_subplot(gs[0, 1])
            # Energy plot
            self.ax_energy = self.fig.add_subplot(gs[1, :])
            
            # Initialize wave packet plot
            self.img_wave = None  # Will be created in update_preview
            
            # Initialize position trajectory plot
            self.line_position, = self.ax_position.plot([], [], 'r-o', markersize=3)
            self.ax_position.set_title("Position Expectation")
            self.ax_position.set_xlabel("X Position")
            self.ax_position.set_ylabel("Y Position")
            self.ax_position.grid(True)
            # Set initial axis limits
            self.ax_position.set_xlim(self.app.sim_2d.x_min, self.app.sim_2d.x_max)
            self.ax_position.set_ylim(self.app.sim_2d.y_min, self.app.sim_2d.y_max)
            
            # Initialize energy plot
            self.line_energy_total, = self.ax_energy.plot([], [], 'm-', label='Total Energy')
            self.line_energy_kin, = self.ax_energy.plot([], [], 'b-', label='Kinetic Energy')
            self.line_energy_pot, = self.ax_energy.plot([], [], 'g-', label='Potential Energy')
            self.ax_energy.set_title("Energy Expectation")
            self.ax_energy.set_xlabel("Time")
            self.ax_energy.set_ylabel("Energy")
            self.ax_energy.grid(True)
            self.ax_energy.legend()
            
        else:
            # 1D layout
            gs = gridspec.GridSpec(3, 3, figure=self.fig)
            
            # Wave packet plot
            self.ax_wave = self.fig.add_subplot(gs[0, 0])
            # Position expectation plot
            self.ax_pos = self.fig.add_subplot(gs[0, 1])
            # Momentum expectation plot
            self.ax_mom = self.fig.add_subplot(gs[1, 0])
            # Energy expectation plot
            self.ax_energy = self.fig.add_subplot(gs[1, 1])
            # Probability current plot
            self.ax_current = self.fig.add_subplot(gs[0, 2])
            
            # Transmission/reflection plot
            self.ax_tr = self.fig.add_subplot(gs[2, 2])
            
            # Initialize wave packet plot
            self.line_wave, = self.ax_wave.plot([], [], 'b-', label='Wave Packet')
            self.line_potential, = self.ax_wave.plot([], [], 'r--', label='Potential')
            self.ax_wave.legend()
            self.ax_wave.set_title("Wave Packet")
            self.ax_wave.set_xlabel("Position")
            self.ax_wave.set_ylabel("Probability Density")
            
            # Initialize observables plots
            self.line_pos, = self.ax_pos.plot([], [], 'g-')
            self.ax_pos.set_title("Position Expectation")
            self.ax_pos.set_xlabel("Time")
            self.ax_pos.set_ylabel("<x>")
            self.ax_pos.grid(True)
            
            self.line_mom, = self.ax_mom.plot([], [], 'r-')
            self.ax_mom.set_title("Momentum Expectation")
            self.ax_mom.set_xlabel("Time")
            self.ax_mom.set_ylabel("<p>")
            self.ax_mom.grid(True)
            
            # Initialize energy plot
            self.line_energy_total, = self.ax_energy.plot([], [], 'm-', label='Total Energy')
            self.line_energy_kin, = self.ax_energy.plot([], [], 'b-', label='Kinetic Energy')
            self.line_energy_pot, = self.ax_energy.plot([], [], 'g-', label='Potential Energy')
            self.ax_energy.set_title("Energy Expectation")
            self.ax_energy.set_xlabel("Time")
            self.ax_energy.set_ylabel("Energy")
            self.ax_energy.grid(True)
            self.ax_energy.legend()
            
            # Initialize probability current plot
            self.line_current, = self.ax_current.plot([], [], 'c-', label='Probability Current')
            self.ax_current.set_title("Probability Current j(x)")
            self.ax_current.set_xlabel("Position")
            self.ax_current.set_ylabel("j(x)")
            self.ax_current.grid(True)
            self.ax_current.legend()
            
            # Transmission/reflection plot
            self.line_trans, = self.ax_tr.plot([], [], 'g-', label='Transmission')
            self.line_refl, = self.ax_tr.plot([], [], 'r-', label='Reflection')
            self.ax_tr.set_title("Transmission/Reflection")
            self.ax_tr.set_xlabel("Time")
            self.ax_tr.set_ylabel("Probability")
            self.ax_tr.set_ylim(0, 1)
            self.ax_tr.grid(True)
            self.ax_tr.legend()
        
        # Create interactive controls
        self.create_controls()
        
        # Initialize the preview
        self.update_preview()
        # Apply to main simulation as well
        if self.is_2d:
            self.app.sim_2d.set_potential(self.last_potential_func_2d if self.last_potential_func_2d else (lambda x, y: np.zeros_like(x)))
            self.app.sim_2d.initialize_gaussian_wavepacket(x0=self.position_x_var.get(), y0=self.position_y_var.get(), kx0=self.momentum_x_var.get(), ky0=self.momentum_y_var.get(), sigma=self.width_var.get())
        else:
            self.app.sim_1d.set_potential(self.last_potential_func_1d if self.last_potential_func_1d else (lambda x: np.zeros_like(x)))
            self.app.sim_1d.initialize_gaussian_wavepacket(x0=self.position_x_var.get(), k0=self.momentum_x_var.get(), sigma=self.width_var.get())
        
        # Update the canvas
        self.fig.tight_layout()
        self.canvas.draw()
    
    def create_controls(self):
        """Create the interactive controls with scrollable sections"""
        # Clear existing controls
        for widget in self.controls_frame.winfo_children():
            widget.destroy()
        
        # Create a canvas and scrollbar for scrollable content
        control_canvas = tk.Canvas(self.controls_frame, highlightthickness=0)
        scrollbar = ttk.Scrollbar(self.controls_frame, orient="vertical", command=control_canvas.yview)
        scrollable_frame = ttk.Frame(control_canvas)
        
        # Configure canvas and scrollbar
        scrollable_frame.bind(
            "<Configure>",
            lambda e: control_canvas.configure(scrollregion=control_canvas.bbox("all"))
        )
        control_canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        control_canvas.configure(yscrollcommand=scrollbar.set)
        
        # Pack canvas and scrollbar
        control_canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Main controls title
        ttk.Label(scrollable_frame, text="Interactive Observables", 
                font=("Arial", 12, "bold")).pack(pady=10)
        
        # Wave packet parameters frame
        params_frame = ttk.LabelFrame(scrollable_frame, text="Wave Packet Parameters", padding="10")
        params_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # Initial position X
        ttk.Label(params_frame, text="Initial Position X:").grid(row=0, column=0, sticky=tk.W, pady=5)
        position_x_slider = ttk.Scale(params_frame, from_=-8.0, to=8.0, orient=tk.HORIZONTAL, 
                                     variable=self.position_x_var, length=150,
                                     command=lambda _: self.update_preview())
        position_x_slider.grid(row=0, column=1, padx=5, pady=5)
        ttk.Label(params_frame, textvariable=self.position_x_var, width=5).grid(row=0, column=2)
        
        if self.is_2d:
            # Initial position Y (2D only)
            ttk.Label(params_frame, text="Initial Position Y:").grid(row=1, column=0, sticky=tk.W, pady=5)
            position_y_slider = ttk.Scale(params_frame, from_=-8.0, to=8.0, orient=tk.HORIZONTAL, 
                                        variable=self.position_y_var, length=150,
                                        command=lambda _: self.update_preview())
            position_y_slider.grid(row=1, column=1, padx=5, pady=5)
            ttk.Label(params_frame, textvariable=self.position_y_var, width=5).grid(row=1, column=2)
        
        # Initial momentum X
        row = 1 if not self.is_2d else 2
        ttk.Label(params_frame, text="Initial Momentum X:").grid(row=row, column=0, sticky=tk.W, pady=5)
        momentum_x_slider = ttk.Scale(params_frame, from_=0.0, to=5.0, orient=tk.HORIZONTAL, 
                                     variable=self.momentum_x_var, length=150,
                                     command=lambda _: self.update_preview())
        momentum_x_slider.grid(row=row, column=1, padx=5, pady=5)
        ttk.Label(params_frame, textvariable=self.momentum_x_var, width=5).grid(row=row, column=2)
        
        if self.is_2d:
            # Initial momentum Y (2D only)
            ttk.Label(params_frame, text="Initial Momentum Y:").grid(row=3, column=0, sticky=tk.W, pady=5)
            momentum_y_slider = ttk.Scale(params_frame, from_=-5.0, to=5.0, orient=tk.HORIZONTAL, 
                                        variable=self.momentum_y_var, length=150,
                                        command=lambda _: self.update_preview())
            momentum_y_slider.grid(row=3, column=1, padx=5, pady=5)
            ttk.Label(params_frame, textvariable=self.momentum_y_var, width=5).grid(row=3, column=2)
        
        # Width (sigma)
        row = 2 if not self.is_2d else 4
        ttk.Label(params_frame, text="Width (sigma):").grid(row=row, column=0, sticky=tk.W, pady=5)
        width_slider = ttk.Scale(params_frame, from_=0.1, to=3.0, orient=tk.HORIZONTAL, 
                               variable=self.width_var, length=150,
                               command=lambda _: self.update_preview())
        width_slider.grid(row=row, column=1, padx=5, pady=5)
        ttk.Label(params_frame, textvariable=self.width_var, width=5).grid(row=row, column=2)
        
        # Potential type frame - now using vertical radiobuttons for better spacing
        pot_frame = ttk.LabelFrame(scrollable_frame, text="Potential Type", padding="10")
        pot_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # Define potential options with better vertical arrangement
        options = [
            ("None", "none"),
            ("Barrier", "barrier"),
            ("Double Slit", "double_slit"),
            ("Harmonic Oscillator", "harmonic")
        ]
        
        # Add 2D-specific options
        if self.is_2d:
            options.append(("Circular Barrier", "circular"))
        
        # Create radio buttons in vertical layout
        for i, (text, value) in enumerate(options):
            ttk.Radiobutton(pot_frame, text=text, variable=self.potential_var, 
                           value=value, command=self.update_preview).pack(anchor=tk.W, pady=2)
        
        # Potential parameters frame
        pot_params_frame = ttk.LabelFrame(scrollable_frame, text="Potential Parameters", padding="10")
        pot_params_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # Potential height
        ttk.Label(pot_params_frame, text="Potential Height:").grid(row=0, column=0, sticky=tk.W, pady=5)
        pot_height_slider = ttk.Scale(pot_params_frame, from_=0.1, to=5.0, orient=tk.HORIZONTAL, 
                                    variable=self.pot_height_var, length=150,
                                    command=lambda _: self.update_preview())
        pot_height_slider.grid(row=0, column=1, padx=5, pady=5)
        ttk.Label(pot_params_frame, textvariable=self.pot_height_var, width=5).grid(row=0, column=2)
        # Barrier width
        ttk.Label(pot_params_frame, text="Barrier Width:").grid(row=1, column=0, sticky=tk.W, pady=5)
        barrier_width_slider = ttk.Scale(pot_params_frame, from_=0.1, to=5.0, orient=tk.HORIZONTAL, 
                                    variable=self.barrier_width_var, length=150,
                                    command=lambda _: self.update_preview())
        barrier_width_slider.grid(row=1, column=1, padx=5, pady=5)
        ttk.Label(pot_params_frame, textvariable=self.barrier_width_var, width=5).grid(row=1, column=2)
        # Barrier position
        ttk.Label(pot_params_frame, text="Barrier Position:").grid(row=2, column=0, sticky=tk.W, pady=5)
        barrier_position_slider = ttk.Scale(pot_params_frame, from_=-8.0, to=8.0, orient=tk.HORIZONTAL, 
                                    variable=self.barrier_position_var, length=150,
                                    command=lambda _: self.update_preview())
        barrier_position_slider.grid(row=2, column=1, padx=5, pady=5)
        ttk.Label(pot_params_frame, textvariable=self.barrier_position_var, width=5).grid(row=2, column=2)
        
        # Visualization selector
        viz_options = [
            ("Probability Density |ψ(x)|²", "probability"),
            ("Momentum Space |ψ̃(p)|²", "momentum")
        ]
        if self.is_2d:
            viz_options = [
                ("Probability Density |ψ(x,y)|²", "probability"),
                ("Momentum Space |ψ̃(px,py)|²", "momentum")
            ]
        viz_frame2 = ttk.LabelFrame(scrollable_frame, text="Visualization Type", padding="10")
        viz_frame2.pack(fill=tk.X, padx=5, pady=5)
        for i, (label, value) in enumerate(viz_options):
            ttk.Radiobutton(viz_frame2, text=label, variable=self.visualization_var, value=value, command=self.update_preview).pack(anchor=tk.W, pady=2)
        
        # Apply button
        ttk.Button(scrollable_frame, text="Apply to Simulation", 
                 command=self.apply_to_simulation).pack(pady=10, padx=5, fill=tk.X)
        
        # Reset history button
        ttk.Button(scrollable_frame, text="Reset Observables History", 
                 command=self.reset_history).pack(pady=5, padx=5, fill=tk.X)

        # --- Observables Dashboard ---
        self.dashboard_frame = ttk.LabelFrame(scrollable_frame, text="Current Observables", padding="10")
        self.dashboard_frame.pack(fill=tk.X, padx=5, pady=10)
        self.dashboard_labels = {}
        obs_labels = [
            ("Time", "time"),
            ("Position <x>", "pos_x"),
            ("Position <y>", "pos_y"),
            ("Momentum <p>", "mom_x"),
            ("Momentum <py>", "mom_y"),
            ("Δx", "unc_x"),
            ("Δy", "unc_y"),
            ("Δpx", "unc_px"),
            ("Δpy", "unc_py"),
            ("Δx·Δpx", "heis_x"),
            ("Δy·Δpy", "heis_y"),
            ("Heisenberg Limit (ħ/2)", "heis_limit"),
            ("Transmission T", "trans"),
            ("Reflection R", "refl"),
        ]
        for i, (label, key) in enumerate(obs_labels):
            self.dashboard_labels[key] = ttk.Label(self.dashboard_frame, text=f"{label}: ")
            self.dashboard_labels[key].grid(row=i, column=0, sticky=tk.W)
        # Hide y/py/Δy/Δpy/Δy·Δpy for 1D
        if not self.is_2d:
            for key in ["pos_y", "mom_y", "unc_y", "unc_py", "heis_y"]:
                self.dashboard_labels[key].grid_remove()
        # Heisenberg limit
        self.dashboard_labels["heis_limit"].config(text="Heisenberg Limit (ħ/2): {:.2f}".format(0.5))
        if self.is_2d:
            for key in ["trans", "refl"]:
                self.dashboard_labels[key].grid_remove()

        # --- Formula Info Button ---
        def show_formula_popup():
            popup = tk.Toplevel()
            popup.title("Formulas & Methods")
            text = tk.Text(popup, wrap=tk.WORD, width=70, height=18)
            text.pack(padx=10, pady=10)
            text.insert(tk.END, """
Formulas and Code Locations:

Position expectation:
  <x> = ∫ x |ψ(x)|² dx   (1D)   [calculate_observables in quantum_core.py]
  <x> = ∬ x |ψ(x,y)|² dxdy   (2D)
Momentum expectation:
  <p> = ∫ p |ψ(p)|² dp   (1D)
  <px> = ∬ px |ψ(px,py)|² dpxdpy   (2D)
Uncertainty:
  Δx = sqrt(<x²> - <x>²)
  Δp = sqrt(<p²> - <p>²)
Heisenberg:
  Δx·Δp ≥ ħ/2
  Δx·Δp is calculated and compared to ħ/2 (ħ=1 in code)

Transmission coefficient:
  T = ∫_{x > barrier} |ψ(x)|² dx
Reflection coefficient:
  R = ∫_{x < barrier} |ψ(x)|² dx
  [transmission_coefficient, reflection_coefficient in quantum_core.py]

All calculated in: calculate_observables (quantum_core.py)

Colormap (2D):
  Maps probability density |ψ(x,y)|² to color.
            """)
            text.config(state=tk.DISABLED)
        ttk.Button(scrollable_frame, text="Show Formulas & Methods", command=show_formula_popup).pack(pady=5, padx=5, fill=tk.X)
        # Colormap explanation for 2D
        if self.is_2d:
            ttk.Label(scrollable_frame, text="Colormap: maps probability density |ψ(x,y)|² to color.", foreground="#555").pack(pady=2)
    
    def update_preview(self, *args):
        """Update the preview based on current parameters"""
        # Get parameters
        if self.is_2d:
            pos_x = self.position_x_var.get()
            pos_y = self.position_y_var.get()
            mom_x = self.momentum_x_var.get()
            mom_y = self.momentum_y_var.get()
            width = self.width_var.get()
            potential_type = self.potential_var.get()
            pot_height = self.pot_height_var.get()
            barrier_width = self.barrier_width_var.get()
            barrier_location = self.barrier_position_var.get()
            colormap = self.colormap_var.get()
            viz_type = self.visualization_var.get()
            
            # Set up the potential based on selection
            if potential_type == "none":
                potential_func = lambda x, y: np.zeros_like(x)
            elif potential_type == "barrier":
                def barrier(x, y):
                    barrier_height = pot_height
                    return barrier_height * ((np.abs(x - barrier_location) < barrier_width/2) & 
                                          (np.abs(y) < 10.0))
                potential_func = barrier
            elif potential_type == "double_slit":
                def double_slit(x, y):
                    barrier_height = pot_height
                    slit_width = 1.0
                    slit_separation = 3.0
                    slit1_center_y = slit_separation/2
                    slit2_center_y = -slit_separation/2
                    barrier = barrier_height * ((np.abs(x - barrier_location) < barrier_width/2) & 
                                             (np.abs(y) < 10.0))
                    slit1 = ((np.abs(x - barrier_location) < barrier_width/2) & 
                           (np.abs(y - slit1_center_y) < slit_width/2))
                    slit2 = ((np.abs(x - barrier_location) < barrier_width/2) & 
                           (np.abs(y - slit2_center_y) < slit_width/2))
                    barrier[slit1 | slit2] = 0
                    return barrier
                potential_func = double_slit
            elif potential_type == "circular":
                def circular_barrier(x, y):
                    barrier_height = pot_height
                    barrier_radius = 3.0
                    barrier_thickness = 0.5
                    center_x, center_y = 0.0, 0.0
                    r = np.sqrt((x - center_x)**2 + (y - center_y)**2)
                    return barrier_height * ((r > barrier_radius - barrier_thickness/2) & 
                                           (r < barrier_radius + barrier_thickness/2))
                potential_func = circular_barrier
            elif potential_type == "harmonic":
                def harmonic_potential(x, y):
                    center_x, center_y = 0.0, 0.0
                    freq = pot_height / 5.0  # Scale the frequency
                    return 0.5 * freq * ((x - center_x)**2 + (y - center_y)**2)
                potential_func = harmonic_potential
            
            self.last_potential_func_2d = potential_func
            
            # Initialize the preview simulation
            self.preview_simulation = QuantumWavePacket2D()
            self.preview_simulation.set_potential(potential_func)
            self.preview_simulation.initialize_gaussian_wavepacket(x0=pos_x, y0=pos_y, 
                                                               kx0=mom_x, ky0=mom_y, 
                                                               sigma=width)
            
            # Clear the previous plot data
            self.ax_wave.clear()
            
            # Plot the selected visualization
            if viz_type == "momentum":
                prob = self.preview_simulation.momentum_density()
                title = "Momentum Space |ψ̃(px,py)|² (momentum_density)"
            else:
                prob = self.preview_simulation.probability_density()
                title = "Probability Density |ψ(x,y)|² (probability_density)"
            self.img_wave = self.ax_wave.imshow(prob, 
                                             extent=[self.preview_simulation.x_min, 
                                                    self.preview_simulation.x_max, 
                                                    self.preview_simulation.y_min, 
                                                    self.preview_simulation.y_max],
                                            origin='lower', cmap=colormap, aspect='auto')
            
            # Plot the potential contours if not "none"
            if potential_type != "none":
                self.ax_wave.contour(self.preview_simulation.X, self.preview_simulation.Y, 
                                   self.preview_simulation.potential,
                                   levels=[0.5 * pot_height], colors='red', alpha=0.8)
            
            # Update plot labels
            self.ax_wave.set_title(title)
            self.ax_wave.set_xlabel("X Position")
            self.ax_wave.set_ylabel("Y Position")
            
            # Autoscale axes and color limits
            self.ax_wave.relim()
            self.ax_wave.autoscale_view()
            self.img_wave.set_clim(0, np.max(prob))
            
        else:  # 1D case
            pos_x = self.position_x_var.get()
            mom_x = self.momentum_x_var.get()
            width = self.width_var.get()
            potential_type = self.potential_var.get()
            pot_height = self.pot_height_var.get()
            barrier_width = self.barrier_width_var.get()
            barrier_location = self.barrier_position_var.get()
            viz_type = self.visualization_var.get()
            
            # Set up the potential based on selection
            if potential_type == "none":
                potential_func = lambda x: np.zeros_like(x)
            elif potential_type == "barrier":
                potential_func = lambda x: pot_height * (np.abs(x - barrier_location) < barrier_width/2)
            elif potential_type == "double_slit":
                def double_slit(x):
                    barrier = pot_height * np.ones_like(x)
                    barrier[np.abs(x - (0.7 + barrier_location)) < 0.2] = 0
                    barrier[np.abs(x + (0.7 - barrier_location)) < 0.2] = 0
                    return barrier
                potential_func = double_slit
            elif potential_type == "harmonic":
                potential_func = lambda x: 0.5 * pot_height * x**2
            
            self.last_potential_func_1d = potential_func
            
            # Initialize the preview simulation
            self.preview_simulation = QuantumWavePacket()
            self.preview_simulation.set_potential(potential_func)
            self.preview_simulation.initialize_gaussian_wavepacket(x0=pos_x, k0=mom_x, sigma=width)
            
            # Update the wave packet plot
            self.ax_wave.clear()
            
            # Plot the selected visualization
            if viz_type == "momentum":
                prob = self.preview_simulation.momentum_density()
                title = "Momentum Space |ψ̃(p)|² (momentum_density)"
                xvals = self.preview_simulation.k * self.preview_simulation.hbar
                xlabel = "Momentum p"
            else:
                prob = self.preview_simulation.probability_density()
                title = "Probability Density |ψ(x)|² (probability_density)"
                xvals = self.preview_simulation.x
                xlabel = "Position"
            self.line_wave, = self.ax_wave.plot(xvals, prob, 'b-', label='Wave Packet')
            
            # Scale potential for better visibility
            pot_scale = np.max(prob) / max(1.0, np.max(self.preview_simulation.potential))
            self.line_potential, = self.ax_wave.plot(self.preview_simulation.x, 
                                                  self.preview_simulation.potential * pot_scale, 
                                                  'r--', label='Potential')
            
            # Update plot labels
            self.ax_wave.set_title(title)
            self.ax_wave.set_xlabel(xlabel)
            self.ax_wave.set_ylabel("Probability Density")
            self.ax_wave.legend()
            
            # Autoscale axes
            self.ax_wave.relim()
            self.ax_wave.autoscale_view()
        
        # Update the figure
        self.fig.tight_layout()
        self.canvas.draw()
    
    def apply_to_simulation(self):
        """Apply the current parameters to the main simulation"""
        # Reset history
        self.reset_history()
        
        # Apply parameters to the actual simulation
        if self.is_2d:
            sim = self.app.sim_2d
            sim.set_potential(self.last_potential_func_2d)
            sim.initialize_gaussian_wavepacket(
                x0=self.position_x_var.get(),
                y0=self.position_y_var.get(),
                kx0=self.momentum_x_var.get(),
                ky0=self.momentum_y_var.get(),
                sigma=self.width_var.get()
            )
        else:
            sim = self.app.sim_1d
            sim.set_potential(self.last_potential_func_1d)
            sim.initialize_gaussian_wavepacket(
                x0=self.position_x_var.get(),
                k0=self.momentum_x_var.get(),
                sigma=self.width_var.get()
            )
    
    def step(self):
        """Advance the simulation by one time step"""
        if self.is_2d:
            self.app.sim_2d.step()
            # Record observables
            px = self.app.sim_2d.position_expectation_x[-1]
            py = self.app.sim_2d.position_expectation_y[-1]
            print(f"[DEBUG] Position expectation: ({px}, {py})")
            self.position_history.append((px, py))
            self.energy_history.append(self.app.sim_2d.energy_expectation[-1])
            self.time_history.append(self.app.sim_2d.time)
            # Do NOT update transmission/reflection histories for 2D
        else:
            self.app.sim_1d.step()
            # Record observables
            self.position_history.append(self.app.sim_1d.position_expectation[-1])
            self.momentum_history.append(self.app.sim_1d.momentum_expectation[-1])
            self.energy_history.append(self.app.sim_1d.energy_expectation[-1])
            self.time_history.append(self.app.sim_1d.time)
            # Update transmission/reflection histories (1D only)
            sim = self.app.sim_1d
            barrier_pos = self.barrier_position_var.get()
            barrier_width = self.barrier_width_var.get()
            self.transmission_history.append(sim.transmission_coefficient(barrier_pos, barrier_width))
            self.reflection_history.append(sim.reflection_coefficient(barrier_pos, barrier_width))
    
    def update_plot(self):
        """Update the plot with current simulation state"""
        if self.is_2d:
            # --- 2D plotting block ---
            prob = self.app.sim_2d.probability_density()
            self.img_wave.set_data(prob)
            self.img_wave.set_clim(0, np.max(prob))
            
            # Update the trajectory plot if we have data
            if len(self.position_history) > 0:
                x_pos = [p[0] for p in self.position_history]
                y_pos = [p[1] for p in self.position_history]
                self.line_position.set_data(x_pos, y_pos)
                # Update axis limits to ensure all points are visible
                x_min, x_max = min(x_pos), max(x_pos)
                y_min, y_max = min(y_pos), max(y_pos)
                x_range = x_max - x_min
                y_range = y_max - y_min
                # Add some padding
                padding = 0.1
                self.ax_position.set_xlim(x_min - x_range * padding, x_max + x_range * padding)
                self.ax_position.set_ylim(y_min - y_range * padding, y_max + y_range * padding)
            
            # Update the energy plot
            sim = self.app.sim_2d
            self.line_energy_total.set_data(sim.time_points, sim.energy_expectation)
            self.line_energy_kin.set_data(sim.time_points, sim.kinetic_energy)
            self.line_energy_pot.set_data(sim.time_points, sim.potential_energy)
            self.ax_energy.relim()
            self.ax_energy.autoscale_view()
            
            # Update titles and labels
            viz_type = self.visualization_var.get()
            if viz_type == "momentum":
                title = "Momentum Space |ψ̃(px,py)|² (momentum_density)"
            else:
                title = "Probability Density |ψ(x,y)|² (probability_density)"
            self.ax_wave.set_title(f"{title} - Time: {self.app.sim_2d.time:.2f}")
            
            # Update dashboard values
            t = sim.time
            px = sim.position_expectation_x[-1] if sim.position_expectation_x else 0
            py = sim.position_expectation_y[-1] if sim.position_expectation_y else 0
            dx = sim.uncertainty_x[-1] if sim.uncertainty_x else 0
            dy = sim.uncertainty_y[-1] if sim.uncertainty_y else 0
            dpx = sim.uncertainty_px[-1] if sim.uncertainty_px else 0
            dpy = sim.uncertainty_py[-1] if sim.uncertainty_py else 0
            heis_x = sim.heisenberg_x[-1] if sim.heisenberg_x else 0
            heis_y = sim.heisenberg_y[-1] if sim.heisenberg_y else 0
            
            self.dashboard_labels["time"].config(text=f"Time: {t:.2f}")
            self.dashboard_labels["pos_x"].config(text=f"<x>: {px:.2f}")
            self.dashboard_labels["pos_y"].config(text=f"<y>: {py:.2f}")
            self.dashboard_labels["unc_x"].config(text=f"Δx: {dx:.2f}")
            self.dashboard_labels["unc_y"].config(text=f"Δy: {dy:.2f}")
            self.dashboard_labels["unc_px"].config(text=f"Δpx: {dpx:.2f}")
            self.dashboard_labels["unc_py"].config(text=f"Δpy: {dpy:.2f}")
            self.dashboard_labels["heis_x"].config(text=f"Δx·Δpx: {heis_x:.2f}")
            self.dashboard_labels["heis_y"].config(text=f"Δy·Δpy: {heis_y:.2f}")
            self.dashboard_labels["mom_x"].config(text=f"<px>: -")
            self.dashboard_labels["mom_y"].config(text=f"<py>: -")
            
            # Single canvas draw at the end
            self.canvas.draw()
        else:
            # --- 1D plotting block ---
            if not hasattr(self, 'line_trans') or not hasattr(self, 'line_refl'):
                return
            prob = self.app.sim_1d.probability_density()
            self.line_wave.set_ydata(prob)
            self.canvas.draw()  # Draw after updating wave packet
            
            # Update the observables plots
            self.line_pos.set_data(self.time_history, self.position_history)
            self.ax_pos.relim()
            self.ax_pos.autoscale_view()
            self.canvas.draw()  # Draw after updating position plot
            
            self.line_mom.set_data(self.time_history, self.momentum_history)
            self.ax_mom.relim()
            self.ax_mom.autoscale_view()
            self.canvas.draw()  # Draw after updating momentum plot
            
            # Update the energy plots
            sim = self.app.sim_1d
            self.line_energy_total.set_data(sim.time_points, sim.energy_expectation)
            self.line_energy_kin.set_data(sim.time_points, sim.kinetic_energy)
            self.line_energy_pot.set_data(sim.time_points, sim.potential_energy)
            self.ax_energy.relim()
            self.ax_energy.autoscale_view()
            self.canvas.draw()  # Draw after updating energy plot
            
            viz_type = self.visualization_var.get()
            if viz_type == "momentum":
                prob = self.app.sim_1d.momentum_density()
                title = "Momentum Space |ψ̃(p)|² (momentum_density)"
                xvals = self.app.sim_1d.k * self.app.sim_1d.hbar
                xlabel = "Momentum p"
            else:
                prob = self.app.sim_1d.probability_density()
                title = "Probability Density |ψ(x)|² (probability_density)"
                xvals = self.app.sim_1d.x
                xlabel = "Position"
            self.ax_wave.set_title(f"{title} - Time: {self.app.sim_1d.time:.2f}")
            self.ax_wave.set_xlabel(xlabel)
            
            # --- Update dashboard values ---
            sim = self.app.sim_1d
            t = sim.time
            px = sim.position_expectation[-1] if sim.position_expectation else 0
            mom = sim.momentum_expectation[-1] if sim.momentum_expectation else 0
            dx = sim.uncertainty_x[-1] if sim.uncertainty_x else 0
            dp = sim.uncertainty_p[-1] if sim.uncertainty_p else 0
            heis = sim.heisenberg_product[-1] if sim.heisenberg_product else 0
            trans = self.transmission_history[-1] if self.transmission_history else 0
            refl = self.reflection_history[-1] if self.reflection_history else 0
            self.dashboard_labels["time"].config(text=f"Time: {t:.2f}")
            self.dashboard_labels["pos_x"].config(text=f"<x>: {px:.2f}")
            self.dashboard_labels["mom_x"].config(text=f"<p>: {mom:.2f}")
            self.dashboard_labels["unc_x"].config(text=f"Δx: {dx:.2f}")
            self.dashboard_labels["unc_px"].config(text=f"Δp: {dp:.2f}")
            self.dashboard_labels["heis_x"].config(text=f"Δx·Δp: {heis:.2f}")
            self.dashboard_labels["trans"].config(text=f"T: {trans:.2f}")
            self.dashboard_labels["refl"].config(text=f"R: {refl:.2f}")
            
            # Update the probability current plot (1D)
            if sim.probability_current_history:
                jx = sim.probability_current_history[-1]
                self.line_current.set_data(sim.x, jx)
                self.ax_current.relim()
                self.ax_current.autoscale_view()
                self.canvas.draw()  # Draw after updating current plot
            
            # Update transmission/reflection plot (1D only)
            self.line_trans.set_data(self.time_history, self.transmission_history)
            self.line_refl.set_data(self.time_history, self.reflection_history)
            self.ax_tr.relim()
            self.ax_tr.autoscale_view()
            self.canvas.draw()  # Draw after updating transmission/reflection plot
    
    def reset_history(self):
        """Reset the observables history"""
        self.position_history = []
        self.momentum_history = []
        self.energy_history = []
        self.time_history = []
        
        # Reset the plots
        if self.is_2d:
            self.line_position.set_data([], [])
            self.line_energy_total.set_data([], [])
            self.line_energy_kin.set_data([], [])
            self.line_energy_pot.set_data([], [])
            self.line_trans.set_data([], [])
            self.line_refl.set_data([], [])
        else:
            self.line_pos.set_data([], [])
            self.line_mom.set_data([], [])
            self.line_energy_total.set_data([], [])
            self.line_energy_kin.set_data([], [])
            self.line_energy_pot.set_data([], [])
            self.line_current.set_data([], [])
            self.line_trans.set_data([], [])
            self.line_refl.set_data([], [])
        
        self.canvas.draw()
        
        # Reset transmission/reflection histories
        self.transmission_history = []
        self.reflection_history = []
    
    def reset(self):
        """Reset the simulation"""
        self.update_preview()
        self.apply_to_simulation()

    def update(self):
        print('[DEBUG] InteractiveObservablesFunction.update() called')
        # Call step() on the function itself to calculate observables
        self.step()
        # Update the plot
        self.update_plot()