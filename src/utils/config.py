from dataclasses import dataclass
from typing import Dict, Any
import json
import os

@dataclass
class SimulationConfig:
    """Configuration for quantum simulations"""
    # General parameters
    dt: float = 0.1
    num_steps: int = 1000
    grid_size: int = 1000
    
    # Wave packet parameters
    initial_position: float = 0.0
    initial_momentum: float = 1.0
    width: float = 1.0
    
    # Potential parameters
    barrier_height: float = 1.0
    barrier_width: float = 1.0
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'SimulationConfig':
        """Create a config from a dictionary"""
        return cls(**data)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert config to dictionary"""
        return {
            'dt': self.dt,
            'num_steps': self.num_steps,
            'grid_size': self.grid_size,
            'initial_position': self.initial_position,
            'initial_momentum': self.initial_momentum,
            'width': self.width,
            'barrier_height': self.barrier_height,
            'barrier_width': self.barrier_width
        }
    
    def save(self, filename: str) -> None:
        """Save configuration to file"""
        with open(filename, 'w') as f:
            json.dump(self.to_dict(), f, indent=4)
    
    @classmethod
    def load(cls, filename: str) -> 'SimulationConfig':
        """Load configuration from file"""
        if not os.path.exists(filename):
            return cls()
        
        with open(filename, 'r') as f:
            data = json.load(f)
        return cls.from_dict(data)
