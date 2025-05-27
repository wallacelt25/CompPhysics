# Quantum Wave Packet Simulation

A Python-based application for simulating and visualizing quantum wave packet phenomena in both 1D and 2D.

## Features

- 1D and 2D quantum wave packet simulations
- Multiple simulation scenarios:
  - Free Particle
  - Barrier Tunneling
  - Double Slit
  - Interactive Observables
  - Tunneling Analysis
- Real-time visualization
- Interactive controls
- Customizable parameters

## Installation

1. Clone the repository:
```bash
git clone [repository-url]
cd quantumWaveApp
```

2. Create a virtual environment (recommended):
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

Run the application:
```bash
python main.py
```

### Controls

- **Start Simulation**: Begins the quantum wave packet simulation
- **Stop Simulation**: Pauses the current simulation
- **Reset Simulation**: Resets all parameters to default values

### Simulation Types

1. **Free Particle**: Simulates a quantum particle in free space
2. **Barrier Tunneling**: Demonstrates quantum tunneling through potential barriers
3. **Double Slit**: Shows interference patterns in double-slit experiments
4. **Interactive Observables**: Allows real-time manipulation of quantum observables
5. **Tunneling Analysis**: Provides detailed analysis of tunneling phenomena

## Development

### Project Structure
quantumWaveApp/
├── src/
│ ├── ui/ # User interface components
│ ├── simulation/ # Quantum simulation logic
│ └── utils/ # Utility functions
├── tests/ # Test files
└── main.py # Application entry point
```

### Running Tests

```bash
pytest tests/
```

### Code Style

This project follows PEP 8 guidelines. To check code style:

```bash
pylint src/
```

## Screenshots

![1D Simulation](screenshots/1d_simulation.png)
![2D Simulation](screenshots/2d_simulation.png)

## Troubleshooting

- If you encounter issues with the UI not displaying, ensure you have `tkinter` installed.
- For best performance, use Python 3.8+ and ensure all dependencies are up to date.

## License

This project is licensed under the MIT License.

## Authors

1. Ari Jaya Teguh - 2702403996
2. Wallace Louis Tjang - 2602169705