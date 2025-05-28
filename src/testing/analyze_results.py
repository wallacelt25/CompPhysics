import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from datetime import datetime

class TestResultAnalyzer:
    def __init__(self, results_dir="test_results"):
        self.results_dir = results_dir
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.plots_dir = os.path.join(results_dir, "plots")
        if not os.path.exists(self.plots_dir):
            os.makedirs(self.plots_dir)
    
    def load_latest_results(self):
        """Load the most recent test results"""
        files = [f for f in os.listdir(self.results_dir) if f.endswith('.csv')]
        if not files:
            raise FileNotFoundError("No test results found")
        
        # Get the most recent timestamp
        latest_file = max(files, key=lambda x: x.split('_')[-1])
        return pd.read_csv(os.path.join(self.results_dir, latest_file))
    
    def analyze_basic_wavepacket(self, df):
        """Analyze basic wave packet test results"""
        wavepacket_df = df[df['test_type'] == 'basic_wavepacket']
        
        # Plot position expectation values
        plt.figure(figsize=(10, 6))
        sns.scatterplot(data=wavepacket_df, x='x0', y='final_x_exp', hue='kx0')
        plt.title('Final X Position vs Initial Position')
        plt.savefig(os.path.join(self.plots_dir, f'position_x_{self.timestamp}.png'))
        plt.close()
        
        # Plot energy vs initial momentum
        plt.figure(figsize=(10, 6))
        sns.scatterplot(data=wavepacket_df, x='kx0', y='final_energy', hue='sigma')
        plt.title('Final Energy vs Initial Momentum')
        plt.savefig(os.path.join(self.plots_dir, f'energy_momentum_{self.timestamp}.png'))
        plt.close()
        
        # Plot uncertainty products
        plt.figure(figsize=(10, 6))
        sns.scatterplot(data=wavepacket_df, x='sigma', y='heisenberg_x', hue='kx0')
        plt.title('Heisenberg Uncertainty Product vs Wave Packet Width')
        plt.savefig(os.path.join(self.plots_dir, f'uncertainty_{self.timestamp}.png'))
        plt.close()
    
    def analyze_grid_resolution(self, df):
        """Analyze grid resolution test results"""
        grid_df = df[df['test_type'] == 'grid_resolution']
        
        # Plot computation time vs grid size
        plt.figure(figsize=(10, 6))
        sns.scatterplot(data=grid_df, x='nx', y='computation_time', hue='dt')
        plt.title('Computation Time vs Grid Size')
        plt.savefig(os.path.join(self.plots_dir, f'computation_time_{self.timestamp}.png'))
        plt.close()
        
        # Plot energy conservation
        plt.figure(figsize=(10, 6))
        sns.scatterplot(data=grid_df, x='dt', y='energy_conservation', hue='nx')
        plt.title('Energy Conservation vs Time Step')
        plt.savefig(os.path.join(self.plots_dir, f'energy_conservation_{self.timestamp}.png'))
        plt.close()
    
    def analyze_potential_tests(self, df):
        """Analyze potential test results"""
        potential_df = df[df['test_type'] == 'potential']
        
        # Plot energy vs potential type
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=potential_df, x='potential', y='final_energy')
        plt.title('Energy Distribution by Potential Type')
        plt.savefig(os.path.join(self.plots_dir, f'potential_energy_{self.timestamp}.png'))
        plt.close()
        
        # Plot energy conservation
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=potential_df, x='potential', y='energy_conservation')
        plt.title('Energy Conservation by Potential Type')
        plt.savefig(os.path.join(self.plots_dir, f'potential_conservation_{self.timestamp}.png'))
        plt.close()
    
    def generate_summary_report(self, df):
        """Generate a summary report of all test results"""
        report = []
        
        # Basic wave packet statistics
        wavepacket_df = df[df['test_type'] == 'basic_wavepacket']
        report.append("Basic Wave Packet Statistics:")
        report.append(f"Average final energy: {wavepacket_df['final_energy'].mean():.4f}")
        report.append(f"Average uncertainty product: {wavepacket_df['heisenberg_x'].mean():.4f}")
        report.append("")
        
        # Grid resolution statistics
        grid_df = df[df['test_type'] == 'grid_resolution']
        report.append("Grid Resolution Statistics:")
        report.append(f"Average computation time: {grid_df['computation_time'].mean():.4f} seconds")
        report.append(f"Average energy conservation: {grid_df['energy_conservation'].mean():.4f}")
        report.append("")
        
        # Potential test statistics
        potential_df = df[df['test_type'] == 'potential']
        report.append("Potential Test Statistics:")
        for pot in potential_df['potential'].unique():
            pot_data = potential_df[potential_df['potential'] == pot]
            report.append(f"{pot}:")
            report.append(f"  Average energy: {pot_data['final_energy'].mean():.4f}")
            report.append(f"  Energy conservation: {pot_data['energy_conservation'].mean():.4f}")
        
        # Save report
        with open(os.path.join(self.results_dir, f'summary_report_{self.timestamp}.txt'), 'w') as f:
            f.write('\n'.join(report))
    
    def analyze_all(self):
        """Run all analysis"""
        df = self.load_latest_results()
        self.analyze_basic_wavepacket(df)
        self.analyze_grid_resolution(df)
        self.analyze_potential_tests(df)
        self.generate_summary_report(df)
        print(f"Analysis complete. Results saved to {self.plots_dir}")

if __name__ == "__main__":
    analyzer = TestResultAnalyzer()
    analyzer.analyze_all() 