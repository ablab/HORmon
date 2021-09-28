from setuptools import setup, find_packages

setup(
    name="HORmon",
    version="1.0",
    description="Centromere annotation tool",
    packages=find_packages(include=['HORmon', 'HORmon.*']),
    entry_points={
        'console_scripts': ['HORmon=HORmon.HORmon:main', 'monomer_inference=HORmon.monomer_inference:main']
    }
)