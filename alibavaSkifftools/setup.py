from distutils.core import setup, Extension
from glob import glob


setup(name='alibavaSkifftools',
        version='v0.1',
        description='Configure the Alibava Marlin processors',
        author='Jordi Duarte-Campderros',
        author_email='Jordi.Duarte.Campderros@cern.ch',
        url='https://github.com/duartej/postproc-alibava/alibavaSkifftools',
        # See https://docs.python.org/2/distutils/setupscript.html#listing-whole-packages
        # for changes in the package distribution
        package_dir={'alibavaSkifftools':'python'},
        # Additional steering files used as templates
        package_data = { 'alibavaSkifftools': ['steering_files/01-ab_converter.xml'] },
        packages = ['alibavaSkifftools'],
        scripts=['bin/open_sesame'],
        )
