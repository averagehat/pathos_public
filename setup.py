from glob import glob
from setuptools import setup, find_packages
import setuptools
import pathos


# Run setuptools setup
setup(

    name =    pathos.__projectname__,
    version = pathos.__version__,
    packages = find_packages(),
    scripts = glob('bin/*'),
    entry_points = {
        'console_scripts': [
            'pathos_single = pathos.pipeline:main',
            'pathos_sheet =  pathos.pathos_sheet:main',
            'pathos_summary =  pathos.summary:main' 
          ]
        } )
    
