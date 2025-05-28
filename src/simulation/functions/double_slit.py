import numpy as np
import matplotlib.pyplot as plt
from ..quantum_core import QuantumWavePacket, QuantumWavePacket2D

class DoubleSlitFunction:
    """
    Double slit simulation
    """
    def __init__(self, app, is_2d=False):
        self.app = app
        self.is_2d = is_2d
        self.fig = self.app.fig_2d if is_2d else self.app.fig_1d
        self.canvas = self.app.canvas_2d if is_2d else self.app.canvas_1d
    
    def setup(self):
        """Set up double slit simulation"""
        # Clear the figure
        self.fig.clear()
        
        if self.is_2d:
            # 2D setup
            self.ax = self.fig.add_subplot(111)
            
            # Set up a double slit potential
            def double_slit_potential_2d(x, y):
                barrier_height = 5.0
                barrier_width = 0.5
                barrier_length = 20.0  # Full y-extent
                barrier_location = 0.0
                
                slit_width = 1.0
                slit_separation = 3.0
                slit1_center_y = slit_separation/2
                slit2_center_y = -slit_separation/2
                
                # Create a barrier at x = barrier_location
                barrier = barrier_height * ((np.abs(x - barrier_location) < barrier_width/2) & 
                                          (np.abs(y) < barrier_length/2))
                
                # Create slits
                slit1 = ((np.abs(x - barrier_location) < barrier_width/2) & 
                       (np.abs(y - slit1_center_y) < slit_width/2))
                slit2 = ((np.abs(x - barrier_location) < barrier_width/2) & 
                       (np.abs(y - slit2_center_y) < slit_width/2))
                
                # Remove slits from barrier
                barrier[slit1 | slit2] = 0
                
                return barrier
            
            self.app.sim_2d.set_potential(double_slit_potential_2d)
            self.app.sim_2d.initialize_gaussian_wavepacket(x0=-10.0, y0=0.0, kx0=3.0, ky0=0.0, sigma=1.5)
            
            # Plot the initial state
            prob = self.app.sim_2d.probability_density()
            self.img = self.ax.imshow(prob, extent=[self.app.sim_2d.x_min, self.app.sim_2d.x_max, 
                                              self.app.sim_2d.y_min, self.app.sim_2d.y_max],
                                     origin='lower', cmap='viridis', aspect='auto')
            
            # Plot the potential
            self.contour = self.ax.contour(self.app.sim_2d.X, self.app.sim_2d.Y, self.app.sim_2d.potential,
                                          levels=[0.5], colors='red', alpha=0.8)
            
            self.cbar = self.fig.colorbar(self.img, ax=self.ax)
            self.cbar.set_label('Probability Density')
            
            self.ax.set_title("Double Slit Interference (2D)")
            self.ax.set_xlabel("X Position")
            self.ax.set_ylabel("Y Position")
            
        else:
            # 1D setup
            self.ax = self.fig.add_subplot(111)
            
            # Set up a double slit potential
            def double_slit_potential(x):
                barrier_height = 5.0
                barrier_width = 4.0
                slit_width = 0.4
                slit_separation = 1.5
                slit_center = 0.0
                
                # Create a barrier with two slits
                barrier = barrier_height * np.ones_like(x)
                
                # First slit
                slit1_center = slit_center - slit_separation/2
                barrier[np.abs(x - slit1_center) < slit_width/2] = 0
                
                # Second slit
                slit2_center = slit_center + slit_separation/2
                barrier[np.abs(x - slit2_center) < slit_width/2] = 0
                
                return barrier
            
            self.app.sim_1d.set_potential(double_slit_potential)
            self.app.sim_1d.initialize_gaussian_wavepacket(x0=-5.0, k0=3.0, sigma=0.5)
            
            # Plot the initial state
            prob = self.app.sim_1d.probability_density()
            self.line_wave, = self.ax.plot(self.app.sim_1d.x, prob, 'b-', label='Probability Density')
            
            # Scale potential for better visibility
            scaled_potential = self.app.sim_1d.potential/5.0
            self.line_potential, = self.ax.plot(self.app.sim_1d.x, scaled_potential, 'r--', 
                                             label='Potential (scaled)')
            
            self.ax.set_title("Double Slit Interference (1D)")
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
            self.ax.set_title(f"Double Slit Interference (2D) - Time: {self.app.sim_2d.time:.2f}")
            self.canvas.draw()
        else:
            # Update the 1D plot
            prob = self.app.sim_1d.probability_density()
            self.line_wave.set_ydata(prob)
            self.ax.set_title(f"Double Slit Interference (1D) - Time: {self.app.sim_1d.time:.2f}")
            self.canvas.draw()
    
    def reset(self):
        """Reset the simulation"""
        if self.is_2d:
            self.app.sim_2d.initialize_gaussian_wavepacket(x0=-10.0, y0=0.0, kx0=3.0, ky0=0.0, sigma=1.5)
        else:
            self.app.sim_1d.initialize_gaussian_wavepacket(x0=-5.0, k0=3.0, sigma=0.5)
        self.update_plot()

    def update(self):
        print('[DEBUG] DoubleSlitFunction.update() called')
        if self.is_2d:
            self.app.sim_2d.step()
        else:
            self.app.sim_1d.step()
        self.update_plot()