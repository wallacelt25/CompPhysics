import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import tkinter as tk
from tkinter import ttk
import threading
from ..quantum_core import QuantumWavePacket, QuantumWavePacket2D

class TunnelingFunction:
    """
    Tunneling analysis for quantum wave packet
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
        
        # For tracking transmission probability
        self.transmission_probs = []
        self.transmission_times = []
        
        # Analysis in progress flag
        self.analysis_running = False
    
    def setup(self):
        """Set up tunneling analysis"""
        # Clear the figure
        self.fig.clear()
        
        if self.is_2d:
            # 2D setup with 3 subplots
            gs = gridspec.GridSpec(2, 2, figure=self.fig, height_ratios=[3, 1])
            
            # Wave function plot
            self.ax_wave = self.fig.add_subplot(gs[0, 0])
            # Barrier cross-section plot
            self.ax_barrier = self.fig.add_subplot(gs[0, 1])
            # Transmission probability plot
            self.ax_trans = self.fig.add_subplot(gs[1, :])
            
            # Set up a barrier potential
            def barrier_potential_2d(x, y):
                barrier_height = 1.0
                barrier_length = 10.0  # Full y-extent
                
                return barrier_height * ((np.abs(x - self.barrier_location) < self.barrier_width/2) & 
                                       (np.abs(y) < barrier_length/2))
            
            self.app.sim_2d.set_potential(barrier_potential_2d)
            self.app.sim_2d.initialize_gaussian_wavepacket(x0=-5.0, y0=0.0, kx0=2.0, ky0=0.0, sigma=1.0)
            
            # Plot the initial wave function
            prob = self.app.sim_2d.probability_density()
            self.img_wave = self.ax_wave.imshow(prob, 
                                             extent=[self.app.sim_2d.x_min, self.app.sim_2d.x_max, 
                                                   self.app.sim_2d.y_min, self.app.sim_2d.y_max],
                                            origin='lower', cmap='viridis', aspect='auto')
            
            # Draw the barrier
            self.barrier_rectangle = plt.Rectangle(
                (self.barrier_location - self.barrier_width/2, self.app.sim_2d.y_min),
                self.barrier_width, self.app.sim_2d.y_max - self.app.sim_2d.y_min,
                color='red', alpha=0.3)
            self.ax_wave.add_patch(self.barrier_rectangle)
            
            self.ax_wave.set_title("Wave Function")
            self.ax_wave.set_xlabel("X Position")
            self.ax_wave.set_ylabel("Y Position")
            
            # Set up the barrier cross-section plot (integrated along y)
            y_center_idx = self.app.sim_2d.ny // 2
            barrier_cross_section = prob[:, y_center_idx]
            self.line_barrier, = self.ax_barrier.plot(self.app.sim_2d.x, barrier_cross_section, 'b-')
            self.ax_barrier.axvspan(self.barrier_location - self.barrier_width/2, 
                                  self.barrier_location + self.barrier_width/2,
                                  alpha=0.3, color='red')
            
            self.ax_barrier.set_title("Barrier Cross-section (y=0)")
            self.ax_barrier.set_xlabel("X Position")
            self.ax_barrier.set_ylabel("Probability Density")
            
            # Initialize the transmission probability plot
            self.line_trans, = self.ax_trans.plot([], [], 'g-o')
            self.ax_trans.set_title("Transmission Probability")
            self.ax_trans.set_xlabel("Time")
            self.ax_trans.set_ylabel("Probability")
            self.ax_trans.set_ylim(0, 1)
            self.ax_trans.grid(True)
            
        else:
            # 1D setup with 2 subplots
            gs = gridspec.GridSpec(2, 1, figure=self.fig, height_ratios=[3, 1])
            
            # Wave function and potential plot
            self.ax_tunnel = self.fig.add_subplot(gs[0])
            # Transmission probability plot
            self.ax_trans = self.fig.add_subplot(gs[1])
            
            # Set up a barrier potential
            def barrier_potential(x):
                barrier_height = 1.0
                return barrier_height * (np.abs(x - self.barrier_location) < self.barrier_width/2)
            
            self.app.sim_1d.set_potential(barrier_potential)
            self.app.sim_1d.initialize_gaussian_wavepacket(x0=-2.0, k0=1.0, sigma=0.5)
            
            # Plot the initial state and potential
            prob = self.app.sim_1d.probability_density()
            self.line_tunnel, = self.ax_tunnel.plot(self.app.sim_1d.x, prob, 'b-', 
                                                 label='Probability Density')
            self.ax_tunnel.plot(self.app.sim_1d.x, self.app.sim_1d.potential, 'r--', 
                              label='Potential Barrier')
            
            # Add vertical lines to mark the regions
            self.ax_tunnel.axvline(x=self.barrier_location - self.barrier_width/2, color='k', linestyle=':')
            self.ax_tunnel.axvline(x=self.barrier_location + self.barrier_width/2, color='k', linestyle=':')
            
            # Shade the regions
            self.ax_tunnel.axvspan(self.app.sim_1d.x_min, self.barrier_location - self.barrier_width/2, 
                                 alpha=0.2, color='blue', label='Incident')
            self.ax_tunnel.axvspan(self.barrier_location - self.barrier_width/2, 
                                 self.barrier_location + self.barrier_width/2, 
                                 alpha=0.2, color='red', label='Barrier')
            self.ax_tunnel.axvspan(self.barrier_location + self.barrier_width/2, self.app.sim_1d.x_max, 
                                 alpha=0.2, color='green', label='Transmission')
            
            self.ax_tunnel.set_title("Quantum Tunneling Through a Barrier")
            self.ax_tunnel.set_xlabel("Position")
            self.ax_tunnel.set_ylabel("Probability Density")
            self.ax_tunnel.legend(loc='upper right')
            self.ax_tunnel.grid(True)
            
            # Prepare the transmission probability plot
            self.line_trans, = self.ax_trans.plot([], [], 'g-')
            self.ax_trans.set_title("Transmission Probability")
            self.ax_trans.set_xlabel("Time")
            self.ax_trans.set_ylabel("Probability")
            self.ax_trans.set_ylim(0, 1)
            self.ax_trans.grid(True)
        
        # Reset transmission data
        self.transmission_probs = []
        self.transmission_times = []
        
        # Add the "Analyze Tunneling" button to the controls
        # Note: We're clearing all widgets first to ensure we don't add multiple buttons
        for widget in self.controls_frame.winfo_children():
            widget.destroy()
        
        ttk.Label(self.controls_frame, text="Tunneling Analysis", 
                font=("Arial", 10, "bold")).pack(pady=5)
        
        self.analyze_button = ttk.Button(self.controls_frame, text="Analyze Tunneling", 
                                       command=self.analyze_tunneling)
        self.analyze_button.pack(fill=tk.X, pady=5)
        
        # Update the canvas
        self.fig.tight_layout()
        self.canvas.draw()
    
    def analyze_tunneling(self):
        """Run tunneling analysis"""
        if self.analysis_running:
            return  # Prevent multiple analyses running simultaneously
        
        self.analysis_running = True
        
        # Disable the analyze button during analysis
        self.analyze_button.config(state="disabled")
        
        # Reset time and observables
        if self.is_2d:
            self.app.sim_2d.time = 0.0
        else:
            self.app.sim_1d.time = 0.0
        
        self.transmission_probs = []
        self.transmission_times = []
        
        # Create a progress window
        progress_window = tk.Toplevel(self.app.master)
        progress_window.title("Tunneling Analysis Progress")
        progress_window.geometry("300x100")
        
        ttk.Label(progress_window, text="Running tunneling analysis...").pack(pady=10)
        progress_var = tk.DoubleVar()
        progress_bar = ttk.Progressbar(progress_window, variable=progress_var, length=250)
        progress_bar.pack(pady=10)
        
        # Number of steps to simulate
        num_steps = 100
        
        # Function to run the analysis in a separate thread
        def run_analysis():
            try:
                for i in range(num_steps):
                    # Update simulation
                    if self.is_2d:
                        self.app.sim_2d.step()
                        
                        # Calculate probability in the right region (transmission)
                        prob = self.app.sim_2d.probability_density()
                        right_region_mask = self.app.sim_2d.X > (self.barrier_location + self.barrier_width/2)
                        trans_prob = np.sum(prob[right_region_mask]) * self.app.sim_2d.dx * self.app.sim_2d.dy
                        
                        self.transmission_probs.append(trans_prob)
                        self.transmission_times.append(self.app.sim_2d.time)
                    else:
                        self.app.sim_1d.step()
                        
                        # Calculate probability in the right region (transmission)
                        prob = self.app.sim_1d.probability_density()
                        right_region = self.app.sim_1d.x > (self.barrier_location + self.barrier_width/2)
                        trans_prob = np.sum(prob[right_region]) * self.app.sim_1d.dx
                        
                        self.transmission_probs.append(trans_prob)
                        self.transmission_times.append(self.app.sim_1d.time)
                    
                    # Update progress bar
                    progress_var.set((i+1) / num_steps * 100)
                    
                    # Update the plots periodically (every 5 steps)
                    if i % 5 == 0:
                        self.app.master.after(0, lambda: self.update_analysis_plots())
                
                # Final update
                self.app.master.after(0, lambda: self.update_analysis_plots())
                
                # Close progress window
                self.app.master.after(100, progress_window.destroy)
                
                # Re-enable the analyze button
                self.app.master.after(100, lambda: self.analyze_button.config(state="normal"))
                
            finally:
                self.analysis_running = False
        
        # Start analysis in a separate thread
        threading.Thread(target=run_analysis).start()
    
    def update_analysis_plots(self):
        """Update the tunneling analysis plots"""
        if self.is_2d:
            # Update the wave function plot
            prob = self.app.sim_2d.probability_density()
            self.img_wave.set_data(prob)
            self.img_wave.set_clim(0, np.max(prob))
            
            # Update the barrier cross-section plot
            y_center_idx = self.app.sim_2d.ny // 2
            self.line_barrier.set_ydata(prob[:, y_center_idx])
            self.ax_barrier.relim()
            self.ax_barrier.autoscale_view()
            
            # Update the transmission probability plot
            self.line_trans.set_data(self.transmission_times, self.transmission_probs)
            self.ax_trans.relim()
            self.ax_trans.autoscale_view()
            
            self.ax_wave.set_title(f"Wave Function - Time: {self.app.sim_2d.time:.2f}")
        else:
            # Update the wave function plot
            prob = self.app.sim_1d.probability_density()
            self.line_tunnel.set_ydata(prob)
            
            # Update the transmission probability plot
            self.line_trans.set_data(self.transmission_times, self.transmission_probs)
            self.ax_trans.relim()
            self.ax_trans.autoscale_view()
            
            self.ax_tunnel.set_title(f"Quantum Tunneling - Time: {self.app.sim_1d.time:.2f}")
        
        # Update the canvas
        self.canvas.draw()
    
    def step(self):
        """Advance the simulation by one time step"""
        if not self.analysis_running:  # Don't step if an analysis is running
            if self.is_2d:
                self.app.sim_2d.step()
                
                # Calculate transmission probability
                prob = self.app.sim_2d.probability_density()
                right_region_mask = self.app.sim_2d.X > (self.barrier_location + self.barrier_width/2)
                trans_prob = np.sum(prob[right_region_mask]) * self.app.sim_2d.dx * self.app.sim_2d.dy
                
                self.transmission_probs.append(trans_prob)
                self.transmission_times.append(self.app.sim_2d.time)
            else:
                self.app.sim_1d.step()
                
                # Calculate transmission probability
                prob = self.app.sim_1d.probability_density()
                right_region = self.app.sim_1d.x > (self.barrier_location + self.barrier_width/2)
                trans_prob = np.sum(prob[right_region]) * self.app.sim_1d.dx
                
                self.transmission_probs.append(trans_prob)
                self.transmission_times.append(self.app.sim_1d.time)
    
    def update_plot(self):
        """Update the plot with current simulation state"""
        if not self.analysis_running:  # Don't update if an analysis is running
            if self.is_2d:
                # Update the wave function plot
                prob = self.app.sim_2d.probability_density()
                self.img_wave.set_data(prob)
                self.img_wave.set_clim(0, np.max(prob))
                
                # Update the barrier cross-section plot
                y_center_idx = self.app.sim_2d.ny // 2
                self.line_barrier.set_ydata(prob[:, y_center_idx])
                self.ax_barrier.relim()
                self.ax_barrier.autoscale_view()
                
                # Update the transmission probability plot
                self.line_trans.set_data(self.transmission_times, self.transmission_probs)
                self.ax_trans.relim()
                self.ax_trans.autoscale_view()
                
                self.ax_wave.set_title(f"Wave Function - Time: {self.app.sim_2d.time:.2f}")
                self.canvas.draw()
            else:
                # Update the wave function plot
                prob = self.app.sim_1d.probability_density()
                self.line_tunnel.set_ydata(prob)
                
                # Update the transmission probability plot
                self.line_trans.set_data(self.transmission_times, self.transmission_probs)
                self.ax_trans.relim()
                self.ax_trans.autoscale_view()
                
                self.ax_tunnel.set_title(f"Quantum Tunneling - Time: {self.app.sim_1d.time:.2f}")
                self.canvas.draw()
    
    def reset(self):
        """Reset the simulation"""
        if self.is_2d:
            self.app.sim_2d.initialize_gaussian_wavepacket(x0=-5.0, y0=0.0, kx0=2.0, ky0=0.0, sigma=1.0)
        else:
            self.app.sim_1d.initialize_gaussian_wavepacket(x0=-2.0, k0=1.0, sigma=0.5)
        
        # Reset transmission data
        self.transmission_probs = []
        self.transmission_times = []
        
        # Update the plots
        self.update_plot()

    def update(self):
        print('[DEBUG] TunnelingFunction.update() called')
        if self.is_2d:
            self.app.sim_2d.step()
        else:
            self.app.sim_1d.step()
        self.update_plot()