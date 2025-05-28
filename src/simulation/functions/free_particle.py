import numpy as np
import matplotlib.pyplot as plt
from ..quantum_core import QuantumWavePacket, QuantumWavePacket2D

class FreqParticleFunction:
    """
    Free particle simulation
    """
    def __init__(self, app, is_2d=False):
        self.app = app
        self.is_2d = is_2d
        self.fig = self.app.fig_2d if is_2d else self.app.fig_1d
        self.canvas = self.app.canvas_2d if is_2d else self.app.canvas_1d
    
    def setup(self):
        """Set up free particle simulation"""
        # Clear the figure
        self.fig.clear()
        if self.is_2d:
            # 2D setup
            self.app.sim_2d.set_potential(lambda x, y: np.zeros_like(x))
            self.app.sim_2d.initialize_gaussian_wavepacket(x0=-5.0, y0=0.0, kx0=2.0, ky0=0.0, sigma=1.0)
            # Create plot
            self.ax = self.fig.add_subplot(111)
            prob = self.app.sim_2d.probability_density()
            self.img = self.ax.imshow(prob, extent=[self.app.sim_2d.x_min, self.app.sim_2d.x_max, 
                                              self.app.sim_2d.y_min, self.app.sim_2d.y_max],
                                     origin='lower', cmap='viridis', aspect='auto')
            self.cbar = self.fig.colorbar(self.img, ax=self.ax)
            self.cbar.set_label('Probability Density')
            self.ax.set_title("Free Particle Wave Packet (2D)")
            self.ax.set_xlabel("X Position")
            self.ax.set_ylabel("Y Position")
        else:
            # 1D setup
            self.app.sim_1d.set_potential(lambda x: np.zeros_like(x))
            self.app.sim_1d.initialize_gaussian_wavepacket(x0=-5, k0=2.0, sigma=0.5)
            # Create plot
            self.ax = self.fig.add_subplot(111)
            prob = self.app.sim_1d.probability_density()
            self.line, = self.ax.plot(self.app.sim_1d.x, prob, 'b-', label='Probability Density')
            self.ax.set_title("Free Particle Wave Packet (1D)")
            self.ax.set_xlabel("Position")
            self.ax.set_ylabel("Probability Density")
            self.ax.set_ylim(0, np.max(prob) * 1.2)
            self.ax.legend()
            self.ax.grid(True)
        # Update the canvas
        self.fig.tight_layout()
        self.canvas.draw()
    
    def step(self):
        """Advance the simulation by one time step"""
        if self.is_2d:
            self.app.sim_2d.step()
        else:
            self.app.sim_1d.step()
    
    def update_plot(self):
        """Update the plot with current simulation state"""
        if self.is_2d:
            # Update the 2D plot
            prob = self.app.sim_2d.probability_density()
            self.img.set_data(prob)
            self.img.set_clim(0, np.max(prob))
            self.ax.set_title(f"Free Particle Wave Packet (2D) - Time: {self.app.sim_2d.time:.2f}")
            self.canvas.draw()
        else:
            # Update the 1D plot
            prob = self.app.sim_1d.probability_density()
            self.line.set_ydata(prob)
            self.ax.set_title(f"Free Particle Wave Packet (1D) - Time: {self.app.sim_1d.time:.2f}")
            self.canvas.draw()
    
    def reset(self):
        """Reset the simulation"""
        if self.is_2d:
            self.app.sim_2d.initialize_gaussian_wavepacket(x0=-5.0, y0=0.0, kx0=2.0, ky0=0.0, sigma=1.0)
        else:
            self.app.sim_1d.initialize_gaussian_wavepacket(x0=-5, k0=2.0, sigma=0.5)
        self.update_plot()

    def update(self):
        print('[DEBUG] FreqParticleFunction.update() called')
        if self.is_2d:
            self.app.sim_2d.step()
        else:
            self.app.sim_1d.step()
        self.update_plot()