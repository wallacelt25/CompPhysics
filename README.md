# Quantum Wave Packet Simulation

A comprehensive Python-based application for simulating and visualizing quantum wave packet phenomena in both 1D and 2D. This interactive educational tool provides real-time quantum mechanics simulations with a modern dark-themed interface. For our project, you only need to focus on the 2D simulation because that is the goal of our project.

## 🚀 Quick Start

**The easiest way to run this application is through the Jupyter notebook launcher:**

1. **Open the Jupyter Notebook**: `QuantumSimulation_Launcher.ipynb`
2. **Run All Cells**: Click "Run All" in your Jupyter environment (Jupyter Lab, Jupyter Notebook, or VS Code)
3. **Follow the automated setup**: The notebook will install dependencies and launch the application automatically

That's it! The notebook handles everything for you.

## 📋 Prerequisites

- Python 3.8 or higher
- Jupyter Lab, Jupyter Notebook, or VS Code with Jupyter extension
- Internet connection (for automatic dependency installation)

## 🎯 Features

### Simulation Types
- **Free Particle**: Quantum particle evolution in free space
- **Barrier Tunneling**: Quantum tunneling through potential barriers
- **Double Slit**: Wave-particle duality and interference patterns
- **Interactive Observables**: Real-time quantum observables monitoring
- **Tunneling Analysis**: Detailed transmission and reflection analysis

### Advanced Features
- **1D and 2D Simulations**: Complete quantum wave packet evolution (please focus on the 2D simulation)
- **Multiple Potential Types**: Barrier, double slit, harmonic oscillator, circular barriers
- **Real-time Visualization**: Live plotting of wave functions and probability densities
- **Interactive Controls**: Adjust parameters with immediate visual feedback
- **Observables Dashboard**: Position, momentum, energy, and uncertainty tracking
- **Modern UI**: Dark theme interface optimized for scientific visualization

## 📖 How to Use the Jupyter Notebook

### Step 1: Launch Jupyter
Open your preferred Jupyter environment:
- **Jupyter Lab**: `jupyter lab`
- **Jupyter Notebook**: `jupyter notebook`
- **VS Code**: Open the `.ipynb` file directly

### Step 2: Open the Launcher
Navigate to and open `QuantumSimulation_Launcher.ipynb`

### Step 3: Run All Cells
- **Jupyter Lab/Notebook**: Click "Run" → "Run All Cells"
- **VS Code**: Click "Run All" at the top of the notebook

### Step 4: Follow the Automated Process
The notebook will:
1. ✅ Install all required dependencies automatically
2. ✅ Check your system compatibility
3. ✅ Launch the quantum simulation application
4. ✅ Provide troubleshooting if needed

### Step 5: Use the Application
Once launched:
1. Choose between **1D** and **2D** simulation tabs (for our project we only use the **2D simulation**)
2. Select **"Interactive Observables"** for the full feature set
3. Adjust potential types and parameters
4. Click **"Reset Simulation"** when changing potential types
5. Use **"Start Simulation"** to begin the animation


## 🎮 Application Controls

### Main Interface
- **Start Simulation**: Begin the quantum wave packet evolution
- **Stop Simulation**: Pause the current simulation
- **Reset Simulation**: Apply new parameters and reset to initial state

### Interactive Observables Features
- **Real-time Parameter Adjustment**: Sliders for position, momentum, and width
- **Potential Type Selection**: Choose from multiple quantum potentials
- **Live Preview**: See changes immediately before applying
- **Observables Dashboard**: Monitor quantum mechanical properties
- **Warning System**: Alerts when potential changes need reset


### Key Quantum Concepts Demonstrated
- Wave-particle duality
- Quantum tunneling
- Heisenberg uncertainty principle
- Quantum interference
- Energy conservation in quantum systems
- Probability density evolution

## 🛠️ Technical Details

### Built With
- **Python**: Core programming language
- **NumPy**: Numerical computations and FFT operations
- **Matplotlib**: Scientific visualization and plotting
- **Tkinter**: Cross-platform GUI framework
- **Numba**: High-performance numerical optimization
- **SciPy**: Advanced scientific computing

### Simulation Methods
- Finite difference time evolution
- Split-operator method for 2D simulations
- Real-time Fourier transforms for momentum space
- Optimized numerical algorithms for performance

## 📁 Project Structure

```
CompPhysics/
├── QuantumSimulation_Launcher.ipynb  # 🎯 Main launcher (START HERE)
├── main.py                           # Direct application entry
├── requirements.txt                  # Python dependencies
├── src/                             # Source code
│   ├── app.py                       # Main application
│   ├── simulation/                  # Quantum simulation engine
│   │   ├── quantum_core.py         # Core quantum mechanics
│   │   └── functions/              # Simulation scenarios
│   └── utils/                      # Utility functions
└── README.md                       # This file
```

## 👥 Authors

- **Wallace Louis Tjang** - 2602169705
- **Ari Jaya Teguh** - 2702403996

**🎯 Remember: Just open `QuantumSimulation_Launcher.ipynb` and click "Run All" to get started!**
