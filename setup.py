from setuptools import setup, find_packages

setup(name='biopython-utils',
      author='Guillaume Poirier-Morency',
      author_email='guillaumepoiriermorency@gmail.com',
      packages=find_packages(),
      install_requires=['biopython', 'pandas'])
