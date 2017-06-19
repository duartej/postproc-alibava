from distutils.core import setup, Extension


setup(name='alibavaSkifftools',
        version='v0.1',
        description='Configure the Alibava Marlin processors',
        author='Jordi Duarte-Campderros',
        author_email='Jordi.Duarte.Campderros@cern.ch',
        url='https://github.com/duartej/postproc-alibava/alibavaSkifftools',
        # See https://docs.python.org/2/distutils/setupscript.html#listing-whole-packages
        # for changes in the package distribution
        package_dir={'alibavaSkifftools':'python'},
        packages = ['alibavaSkifftools'],
        scripts=['bin/open_sesame'],
        )
