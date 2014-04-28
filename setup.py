from distutils.core import setup, Extension

setup(
    name='Amplishot',
    version='0.3.3',
    author='Connor Skennerton',
    author_email='c.tkennerton@gmail.com',
    packages=['amplishot', 'amplishot.test', 'amplishot.search',
        'amplishot.parse', 'amplishot.app'],
    package_data={'amplishot': ['data/logWordPrior.txt']},
    scripts=['bin/amplishot', 'bin/link16SToContigs.py'],
    url='http://pypi.python.org/pypi/Amplishot/',
    license='LICENSE.txt',
    description='Amplishot',
    long_description=open('README.md').read(),
    requires=['biopython'],
    ext_modules = [Extension('_wumanber',['amplishot/search/WuManber.cpp'])]
)
