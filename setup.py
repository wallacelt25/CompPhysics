from setuptools import setup, find_packages

setup(
    name="quantumWaveApp",
    version="1.0.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.21.0",
        "matplotlib>=3.4.0",
        "pytest>=7.0.0",
        "pylint>=2.12.0",
        "mypy>=0.910",
    ],
    author="[Your Name]",
    author_email="[Your Email]",
    description="A quantum wave packet simulation application",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="[Your Repository URL]",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)
