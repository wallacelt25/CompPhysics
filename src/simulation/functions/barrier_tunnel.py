import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from ..quantum_core import QuantumWavePacket, QuantumWavePacket2D

class BarrierTunnelingFunction:
    """
    Barrier tunneling simulation
    """
    def __init__(self, app, is_2d=False):
        self.app = app
        self.is_2d = is_2d
        self.fig = self.app.fig_2d if is_2d else self.app.fig_1d
        self.canvas = self.app.canvas_2d if is_2d else self.app.canvas_1d
        self.controls_frame = self.app.function_controls_frame_2d if is_2d else self.app.function_controls_frame_1d
        
        # Barrier parameters
        self.barrier_location = 0.0
        self.barrier_width = 0.5
        self.barrier_height = 1.0
        
        # For tracking transmission probability
        self.transmission_probs = []
        self.transmission_times = []
    
    def setup(self):
        """Set up barrier tunneling simulation"""
        # Clear the figure
        self.fig.clear()
        if self.is_2d:
            # 2D setup with two subplots (wave function and transmission probability)
            gs = gridspec.GridSpec(1, 2, figure=self.fig, width_ratios=[4, 1])
            self.ax_wave = self.fig.add_subplot(gs[0])
            def barrier_potential_2d(x, y):
                return self.barrier_height * ((np.abs(x - self.barrier_location) < self.barrier_width/2) & 
                                           (np.abs(y) < 10.0))
            self.app.sim_2d.set_potential(barrier_potential_2d)
            self.app.sim_2d.initialize_gaussian_wavepacket(x0=-5.0, y0=0.0, kx0=2.0, ky0=0.0, sigma=1.0)
            prob = self.app.sim_2d.probability_density()
            self.img_wave = self.ax_wave.imshow(prob, 
                                             extent=[self.app.sim_2d.x_min, self.app.sim_2d.x_max, 
                                                   self.app.sim_2d.y_min, self.app.sim_2d.y_max],
                                            origin='lower', cmap='viridis', aspect='auto')
            self.barrier_rectangle = plt.Rectangle(
                (self.barrier_location - self.barrier_width/2, self.app.sim_2d.y_min),
                self.barrier_width, self.app.sim_2d.y_max - self.app.sim_2d.y_min,
                color='red', alpha=0.3)
            self.ax_wave.add_patch(self.barrier_rectangle)
            self.cbar = self.fig.colorbar(self.img_wave, ax=self.ax_wave)
            self.cbar.set_label('Probability Density')
            self.ax_wave.set_title("Barrier Tunneling (2D)")
            self.ax_wave.set_xlabel("X Position")
            self.ax_wave.set_ylabel("Y Position")
            self.ax_trans = self.fig.add_subplot(gs[1])
            self.line_trans, = self.ax_trans.plot([], [], 'g-o')
            self.ax_trans.set_title("Transmission")
            self.ax_trans.set_xlabel("Time")
            self.ax_trans.set_ylabel("Probability")
            self.ax_trans.set_ylim(0, 1)
            self.transmission_probs = []
            self.transmission_times = []
        else:
            self.ax = self.fig.add_subplot(111)
            def barrier_potential(x):
                return self.barrier_height * (np.abs(x - self.barrier_location) < self.barrier_width/2)
            self.app.sim_1d.set_potential(barrier_potential)
            self.app.sim_1d.initialize_gaussian_wavepacket(x0=-2.0, k0=1.0, sigma=0.5)
            prob = self.app.sim_1d.probability_density()
            self.line_wave, = self.ax.plot(self.app.sim_1d.x, prob, 'b-', label='Probability Density')
            self.line_potential, = self.ax.plot(self.app.sim_1d.x, self.app.sim_1d.potential, 'r--', 
                                             label='Potential Barrier')
            self.ax.set_title("Quantum Tunneling Through a Barrier (1D)")
            self.ax.set_xlabel("Position")
            self.ax.set_ylabel("Probability Density")
            self.ax.set_ylim(0, np.max(prob) * 1.2)
            self.ax.legend()
            self.ax.grid(True)
        self.fig.tight_layout()
        self.canvas.draw()
    
    def step(self):
        """Advance the simulation by one time step"""
        if self.is_2d:
            self.app.sim_2d.step()
            
            # Calculate transmission probability (2D)
            prob = self.app.sim_2d.probability_density()
            right_region_mask = self.app.sim_2d.X > (self.barrier_location + self.barrier_width/2)
            trans_prob = np.sum(prob[right_region_mask]) * self.app.sim_2d.dx * self.app.sim_2d.dy
            
            self.transmission_probs.append(trans_prob)
            self.transmission_times.append(self.app.sim_2d.time)
        else:
            self.app.sim_1d.step()
    
    def update_plot(self):
        """Update the plot with current simulation state"""
        if self.is_2d:
            # Update the 2D plot
            prob = self.app.sim_2d.probability_density()
            self.img_wave.set_data(prob)
            self.img_wave.set_clim(0, np.max(prob))
            
            # Update transmission probability plot
            if len(self.transmission_times) > 0:
                self.line_trans.set_data(self.transmission_times, self.transmission_probs)
                self.ax_trans.set_xlim(min(self.transmission_times), max(self.transmission_times))
                self.ax_trans.set_ylim(0, max(1.0, max(self.transmission_probs) * 1.1))
                self.ax_trans.relim()
                self.ax_trans.autoscale_view()
            
            self.ax_wave.set_title(f"Barrier Tunneling (2D) - Time: {self.app.sim_2d.time:.2f}")
            self.canvas.draw()
        else:
            # Update the 1D plot
            prob = self.app.sim_1d.probability_density()
            self.line_wave.set_ydata(prob)
            self.ax.set_title(f"Quantum Tunneling - Time: {self.app.sim_1d.time:.2f}")
            self.canvas.draw()
    
    def reset(self):
        """Reset the simulation"""
        if self.is_2d:
            self.app.sim_2d.initialize_gaussian_wavepacket(x0=-5.0, y0=0.0, kx0=2.0, ky0=0.0, sigma=1.0)
            self.transmission_probs = []
            self.transmission_times = []
            self.line_trans.set_data([], [])
        else:
            self.app.sim_1d.initialize_gaussian_wavepacket(x0=-2.0, k0=1.0, sigma=0.5)
        self.update_plot()

    def update(self):
        print('[DEBUG] BarrierTunnelingFunction.update() called')
        # Call step() on the function itself to calculate transmission
        self.step()
        # Update the plot
        self.update_plot()