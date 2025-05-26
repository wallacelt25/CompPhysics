import numpy as np
from scipy.fftpack import fft, ifft, fftfreq
from numba import jit, complex128, float64
import numba as nb

# Remove jit from fft2/ifft2, use numpy directly

def fft2(x):
    """2D FFT using numpy"""
    return np.fft.fft2(x)

def ifft2(x):
    """2D inverse FFT using numpy"""
    return np.fft.ifft2(x)

def calculate_probability_density(psi):
    """Optimized probability density calculation"""
    return np.abs(psi)**2

def calculate_momentum_density(psi_k):
    """Optimized momentum density calculation"""
    return np.abs(psi_k)**2

def calculate_probability_current(psi, dx, hbar, mass):
    """Optimized probability current calculation"""
    dpsi_dx = np.gradient(psi, dx)
    return (hbar / mass) * np.imag(np.conj(psi) * dpsi_dx)

class QuantumWavePacket:
    """
    Core class for quantum wave packet simulation and evolution
    using the split-step Fourier method to solve the time-dependent Schrödinger equation.
    Optimized for performance using numba and vectorized operations.
    """
    def __init__(self, grid_points=1000, x_min=-10, x_max=10, dt=0.01, mass=1.0, hbar=1.0):
        # Spatial grid parameters
        self.grid_points = grid_points
        self.x_min = x_min
        self.x_max = x_max
        self.x = np.linspace(x_min, x_max, grid_points)
        self.dx = (x_max - x_min) / (grid_points - 1)
        
        # Time parameters
        self.dt = dt
        self.time = 0.0
        
        # Physical constants
        self.mass = mass
        self.hbar = hbar
        
        # Pre-compute momentum space and kinetic energy operator
        self.k = 2 * np.pi * fftfreq(grid_points, self.dx)
        self.kinetic_operator = np.exp(-1j * self.hbar * self.k**2 * self.dt / (2 * self.mass))
        
        # Initialize wave function to None
        self.psi = None
        self.potential = np.zeros(grid_points)
        
        # Pre-allocate arrays for observables
        self.initialize_observable_arrays()
        
    def initialize_observable_arrays(self):
        """Pre-allocate arrays for storing observables"""
        self.position_expectation = []
        self.momentum_expectation = []
        self.energy_expectation = []
        self.time_points = []
        self.uncertainty_x = []
        self.uncertainty_p = []
        self.heisenberg_product = []
        self.kinetic_energy = []
        self.potential_energy = []
        self.probability_current_history = []
        
    def initialize_gaussian_wavepacket(self, x0=0, k0=1, sigma=1.0, normalize=True):
        """Initialize a Gaussian wave packet with specified parameters"""
        # Vectorized initialization
        self.psi = np.exp(-(self.x - x0)**2 / (2 * sigma**2)) * np.exp(1j * k0 * self.x)
        
        if normalize:
            norm = np.sqrt(np.sum(np.abs(self.psi)**2) * self.dx)
            self.psi = self.psi / norm
            
        self.time = 0.0
        self.initialize_observable_arrays()
        
    def set_potential(self, potential_func):
        """Set the potential energy function"""
        self.potential = potential_func(self.x)
        
    def step(self):
        print(f'[DEBUG] QuantumWavePacket.step() called, time before: {self.time}')
        # Pre-compute potential phase factor
        potential_phase = np.exp(-1j * self.potential * self.dt / (2 * self.hbar))
        
        # Half-step in position space (applying potential)
        self.psi *= potential_phase
        
        # Full step in momentum space using pre-computed kinetic operator
        psi_k = fft(self.psi)
        psi_k *= self.kinetic_operator
        self.psi = ifft(psi_k)
        
        # Another half-step in position space
        self.psi *= potential_phase
        
        # Update time
        self.time += self.dt
        print(f'[DEBUG] QuantumWavePacket.step() called, time after: {self.time}')
        
        # Calculate and store observables
        self.calculate_observables()
        
    def probability_density(self):
        """Calculate the probability density |ψ|² using optimized function"""
        return calculate_probability_density(self.psi)
    
    def momentum_density(self):
        """Calculate the momentum space probability density |ψ̃(p)|² using optimized function"""
        psi_k = fft(self.psi)
        return calculate_momentum_density(psi_k)
    
    def probability_current(self):
        """Calculate the probability current density using optimized function"""
        return calculate_probability_current(self.psi, self.dx, self.hbar, self.mass)
    
    def calculate_observables(self):
        """Calculate expectation values for position, momentum, and energy using vectorized operations"""
        # Probability density
        prob = self.probability_density()
        
        # Position expectation value using vectorized operations
        x_exp = np.sum(self.x * prob) * self.dx
        self.position_expectation.append(x_exp)
        
        # Momentum space wave function
        psi_k = fft(self.psi)
        prob_k = calculate_momentum_density(psi_k)
        
        # Momentum expectation value using vectorized operations
        k_exp = np.sum(self.k * prob_k) * (self.k[1] - self.k[0])
        self.momentum_expectation.append(k_exp * self.hbar)
        
        # Energy calculations using vectorized operations
        kinetic = np.sum((self.hbar**2 * self.k**2 / (2 * self.mass)) * prob_k) * (self.k[1] - self.k[0])
        potential = np.sum(self.potential * prob) * self.dx
        
        self.kinetic_energy.append(kinetic)
        self.potential_energy.append(potential)
        self.energy_expectation.append(kinetic + potential)
        self.time_points.append(self.time)
        
        # Probability current profile
        self.probability_current_history.append(self.probability_current())
        
        # Uncertainty calculations using vectorized operations
        x2_exp = np.sum((self.x ** 2) * prob) * self.dx
        dx_val = x2_exp - x_exp ** 2
        dx = np.sqrt(dx_val) if dx_val >= 0 else 0.0
        self.uncertainty_x.append(dx)
        
        p = self.k * self.hbar
        p2_exp = np.sum((p ** 2) * prob_k) * (self.k[1] - self.k[0])
        p_exp = np.sum(p * prob_k) * (self.k[1] - self.k[0])
        dp_val = p2_exp - p_exp ** 2
        dp = np.sqrt(dp_val) if dp_val >= 0 else 0.0
        self.uncertainty_p.append(dp)
        
        # Heisenberg uncertainty product
        self.heisenberg_product.append(dx * dp)

    def transmission_coefficient(self, barrier_position=0.0, barrier_width=1.0):
        """Compute the transmission coefficient (probability to the right of the barrier)"""
        right_region = self.x > (barrier_position + barrier_width/2)
        prob = self.probability_density()
        return np.sum(prob[right_region]) * self.dx

    def reflection_coefficient(self, barrier_position=0.0, barrier_width=1.0):
        """Compute the reflection coefficient (probability to the left of the barrier)"""
        left_region = self.x < (barrier_position - barrier_width/2)
        prob = self.probability_density()
        return np.sum(prob[left_region]) * self.dx


class QuantumWavePacket2D:
    """
    Class for 2D quantum wave packet simulation using the split-step Fourier method.
    Optimized for performance using numba and vectorized operations.
    """
    def __init__(self, nx=100, ny=100, x_min=-10, x_max=10, y_min=-10, y_max=10, dt=0.01, mass=1.0, hbar=1.0):
        # Spatial grid parameters
        self.nx = nx
        self.ny = ny
        self.x_min, self.x_max = x_min, x_max
        self.y_min, self.y_max = y_min, y_max
        
        self.x = np.linspace(x_min, x_max, nx)
        self.y = np.linspace(y_min, y_max, ny)
        self.dx = (x_max - x_min) / (nx - 1)
        self.dy = (y_max - y_min) / (ny - 1)
        
        # Create 2D mesh
        self.X, self.Y = np.meshgrid(self.x, self.y)
        
        # Time parameters
        self.dt = dt
        self.time = 0.0
        
        # Physical constants
        self.mass = mass
        self.hbar = hbar
        
        # Pre-compute momentum space and kinetic energy operator
        self.kx = 2 * np.pi * fftfreq(nx, self.dx)
        self.ky = 2 * np.pi * fftfreq(ny, self.dy)
        self.KX, self.KY = np.meshgrid(self.kx, self.ky)
        self.kinetic_operator = np.exp(-1j * self.hbar * (self.KX**2 + self.KY**2) * self.dt / (2 * self.mass))
        
        # Initialize wave function to None
        self.psi = None
        self.potential = np.zeros((ny, nx))
        
        # Pre-allocate arrays for observables
        self.initialize_observable_arrays()
        
    def initialize_observable_arrays(self):
        """Pre-allocate arrays for storing observables"""
        self.position_expectation_x = []
        self.position_expectation_y = []
        self.energy_expectation = []
        self.time_points = []
        self.uncertainty_x = []
        self.uncertainty_y = []
        self.uncertainty_px = []
        self.uncertainty_py = []
        self.heisenberg_x = []
        self.heisenberg_y = []
        self.kinetic_energy = []
        self.potential_energy = []
        self.probability_current_history = []
        
    def initialize_gaussian_wavepacket(self, x0=0, y0=0, kx0=1, ky0=1, sigma=1.0, normalize=True):
        """Initialize a 2D Gaussian wave packet with specified parameters using vectorized operations"""
        r_squared = (self.X - x0)**2 + (self.Y - y0)**2
        self.psi = np.exp(-r_squared / (2 * sigma**2)) * np.exp(1j * (kx0 * self.X + ky0 * self.Y))
        
        if normalize:
            norm = np.sqrt(np.sum(np.abs(self.psi)**2) * self.dx * self.dy)
            self.psi = self.psi / norm
            
        self.time = 0.0
        self.initialize_observable_arrays()
        
    def set_potential(self, potential_func):
        """Set the potential energy function"""
        self.potential = potential_func(self.X, self.Y)
        
    def step(self):
        print(f'[DEBUG] QuantumWavePacket2D.step() called, time before: {self.time}')
        # Pre-compute potential phase factor
        potential_phase = np.exp(-1j * self.potential * self.dt / (2 * self.hbar))
        
        # Half-step in position space (applying potential)
        self.psi *= potential_phase
        
        # Full step in momentum space using pre-computed kinetic operator
        psi_k = fft2(self.psi)
        psi_k *= self.kinetic_operator
        self.psi = ifft2(psi_k)
        
        # Another half-step in position space
        self.psi *= potential_phase
        
        # Update time
        self.time += self.dt
        print(f'[DEBUG] QuantumWavePacket2D.step() called, time after: {self.time}')
        
        # Calculate and store observables
        self.calculate_observables()
        
    def probability_density(self):
        """Calculate the probability density |ψ|² using optimized function"""
        return calculate_probability_density(self.psi)
    
    def momentum_density(self):
        """Calculate the momentum space probability density |ψ̃(p)|² using optimized function"""
        psi_k = fft2(self.psi)
        return calculate_momentum_density(psi_k)
    
    def calculate_observables(self):
        """Calculate expectation values for position, momentum, and energy using vectorized operations"""
        # Probability density
        prob = self.probability_density()
        
        # Position expectation values using vectorized operations
        x_exp = np.sum(self.X * prob) * self.dx * self.dy
        y_exp = np.sum(self.Y * prob) * self.dx * self.dy
        self.position_expectation_x.append(x_exp)
        self.position_expectation_y.append(y_exp)
        
        # Momentum space wave function
        psi_k = fft2(self.psi)
        prob_k = calculate_momentum_density(psi_k)
        
        # Energy calculations using vectorized operations
        kinetic = np.sum((self.hbar**2 * (self.KX**2 + self.KY**2) / (2 * self.mass)) * prob_k) * self.dx * self.dy
        potential = np.sum(self.potential * prob) * self.dx * self.dy
        
        self.kinetic_energy.append(kinetic)
        self.potential_energy.append(potential)
        self.energy_expectation.append(kinetic + potential)
        self.time_points.append(self.time)
        
        # Uncertainty calculations using vectorized operations
        x2_exp = np.sum((self.X ** 2) * prob) * self.dx * self.dy
        y2_exp = np.sum((self.Y ** 2) * prob) * self.dx * self.dy
        
        dx_val = x2_exp - x_exp ** 2
        dy_val = y2_exp - y_exp ** 2
        
        dx = np.sqrt(dx_val) if dx_val >= 0 else 0.0
        dy = np.sqrt(dy_val) if dy_val >= 0 else 0.0
        
        self.uncertainty_x.append(dx)
        self.uncertainty_y.append(dy)
        
        # Momentum uncertainties
        px = self.KX * self.hbar
        py = self.KY * self.hbar
        
        px2_exp = np.sum((px ** 2) * prob_k) * self.dx * self.dy
        py2_exp = np.sum((py ** 2) * prob_k) * self.dx * self.dy
        
        px_exp = np.sum(px * prob_k) * self.dx * self.dy
        py_exp = np.sum(py * prob_k) * self.dx * self.dy
        
        dpx_val = px2_exp - px_exp ** 2
        dpy_val = py2_exp - py_exp ** 2
        
        dpx = np.sqrt(dpx_val) if dpx_val >= 0 else 0.0
        dpy = np.sqrt(dpy_val) if dpy_val >= 0 else 0.0
        
        self.uncertainty_px.append(dpx)
        self.uncertainty_py.append(dpy)
        
        # Heisenberg uncertainty products
        self.heisenberg_x.append(dx * dpx)
        self.heisenberg_y.append(dy * dpy)