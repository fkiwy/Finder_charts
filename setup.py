from setuptools import setup

setup(name='Finder_charts',
      version='1.1.0',
      description='Multi-bands finder charts from image data of various sky surveys',
      url='https://github.com/fkiwy/Finder_charts',
      author='Frank Kiwy',
      author_email='frank.kiwy@outlook.com',
      license='MIT',
      py_modules=['Finder_charts', 'Gaia_finder_chart'],
      install_requires=['numpy', 'matplotlib', 'pillow', 'requests', 'astropy', 'astroquery', 'reproject', 'tdqm'])
