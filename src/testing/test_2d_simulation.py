import numpy as np
import pandas as pd
import time
from datetime import datetime
import os
import tkinter as tk
from src.simulation.quantum_core import QuantumWavePacket2D
from src.simulation.functions.free_particle import FreqParticleFunction
from src.simulation.functions.barrier_tunnel import BarrierTunnelingFunction
from src.simulation.functions.double_slit import DoubleSlitFunction
from src.app import QuantumApp

class Quantum2DTestRunner:
    def __init__(self, output_dir="test_results"):
        self.output_dir = output_dir
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Test parameters
        self.positions = [-5, 0, 5]
        self.momenta = [0.5, 1.0, 2.0]
        self.widths = [0.5, 1.0, 2.0]
        self.grid_sizes = [50, 100, 200]
        self.time_steps = [0.005, 0.01, 0.02]
        
        # Initialize results storage
        self.results = []
        
        # Initialize app for potential functions
        self.root = tk.Tk()
        self.app = QuantumApp(self.root)
        
    def run_basic_wavepacket_tests(self):
        """Run tests for basic wave packet evolution"""
        print("Running basic wave packet tests...")
        
        for x0 in self.positions:
            for y0 in self.positions:
                for kx0 in self.momenta:
                    for ky0 in self.momenta:
                        for sigma in self.widths:
                            # Initialize simulation
                            sim = QuantumWavePacket2D(nx=100, ny=100)
                            sim.initialize_gaussian_wavepacket(x0=x0, y0=y0, 
                                                            kx0=kx0, ky0=ky0, 
                                                            sigma=sigma)
                            
                            # Run simulation for 100 steps
                            start_time = time.time()
                            for _ in range(100):
                                sim.step()
                            computation_time = time.time() - start_time
                            
                            # Record results
                            self.results.append({
                                'test_type': 'basic_wavepacket',
                                'x0': x0,
                                'y0': y0,
                                'kx0': kx0,
                                'ky0': ky0,
                                'sigma': sigma,
                                'final_x_exp': sim.position_expectation_x[-1],
                                'final_y_exp': sim.position_expectation_y[-1],
                                'final_energy': sim.energy_expectation[-1],
                                'computation_time': computation_time,
                                'uncertainty_x': sim.uncertainty_x[-1],
                                'uncertainty_y': sim.uncertainty_y[-1],
                                'heisenberg_x': sim.heisenberg_x[-1],
                                'heisenberg_y': sim.heisenberg_y[-1]
                            })
    
    def run_grid_resolution_tests(self):
        """Run tests for different grid resolutions"""
        print("Running grid resolution tests...")
        
        for nx in self.grid_sizes:
            for ny in self.grid_sizes:
                for dt in self.time_steps:
                    # Initialize simulation
                    sim = QuantumWavePacket2D(nx=nx, ny=ny, dt=dt)
                    sim.initialize_gaussian_wavepacket()
                    
                    # Run simulation for 50 steps
                    start_time = time.time()
                    for _ in range(50):
                        sim.step()
                    computation_time = time.time() - start_time
                    
                    # Record results
                    self.results.append({
                        'test_type': 'grid_resolution',
                        'nx': nx,
                        'ny': ny,
                        'dt': dt,
                        'computation_time': computation_time,
                        'final_energy': sim.energy_expectation[-1],
                        'energy_conservation': np.std(sim.energy_expectation) / np.mean(sim.energy_expectation)
                    })
    
    def run_potential_tests(self):
        """Run tests for different potential types"""
        print("Running potential tests...")
        
        # Test different potential types
        potential_types = {
            'free_particle': FreqParticleFunction,
            'barrier': BarrierTunnelingFunction,
            'double_slit': DoubleSlitFunction
        }
        
        for pot_name, pot_class in potential_types.items():
            for x0 in [-5, 0, 5]:
                # Initialize simulation
                sim = QuantumWavePacket2D()
                sim.initialize_gaussian_wavepacket(x0=x0)
                
                # Set up potential
                pot = pot_class(self.app, is_2d=True)
                pot.setup()  # Initialize the potential
                
                # Get the potential function based on the type
                if pot_name == 'free_particle':
                    potential_func = lambda x, y: np.zeros_like(x)
                elif pot_name == 'barrier':
                    potential_func = lambda x, y: np.where(np.abs(x) < 1, 1.0, 0.0)
                else:  # double_slit
                    potential_func = lambda x, y: np.where((np.abs(x) > 0.5) & (np.abs(x) < 1.5), 1.0, 0.0)
                
                sim.set_potential(potential_func)
                
                # Run simulation for 100 steps
                start_time = time.time()
                for _ in range(100):
                    sim.step()
                computation_time = time.time() - start_time
                
                # Record results
                self.results.append({
                    'test_type': 'potential',
                    'potential': pot_name,
                    'x0': x0,
                    'computation_time': computation_time,
                    'final_energy': sim.energy_expectation[-1],
                    'energy_conservation': np.std(sim.energy_expectation) / np.mean(sim.energy_expectation)
                })
    
    def save_results(self):
        """Save test results to CSV files"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Convert results to DataFrame
        df = pd.DataFrame(self.results)
        
        # Save all results
        df.to_csv(f"{self.output_dir}/all_results_{timestamp}.csv", index=False)
        
        # Save results by test type
        for test_type in df['test_type'].unique():
            test_df = df[df['test_type'] == test_type]
            test_df.to_csv(f"{self.output_dir}/{test_type}_results_{timestamp}.csv", index=False)
        
        print(f"Results saved to {self.output_dir}")
    
    def run_all_tests(self):
        """Run all test categories"""
        try:
            self.run_basic_wavepacket_tests()
            self.run_grid_resolution_tests()
            self.run_potential_tests()
            self.save_results()
        finally:
            # Clean up
            self.root.destroy()

if __name__ == "__main__":
    # Create test runner and execute all tests
    runner = Quantum2DTestRunner()
    runner.run_all_tests() 