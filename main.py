import tkinter as tk
from src.app import QuantumApp

def main():
    """Main function to run the application"""
    root = tk.Tk()
    app = QuantumApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()