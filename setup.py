from distutils.core import setup, Extension

setup(
    name='Amplishot',
    version='0.3.2',
    author='Connor Skennerton',
    author_email='c.tkennerton@gmail.com',
    packages=['amplishot', 'amplishot.test', 'amplishot.search',
        'amplishot.parse', 'amplishot.app'],
    scripts=['bin/amplishot'],
    url='http://pypi.python.org/pypi/Amplishot/',
    license='LICENSE.txt',
    description='Amplishot',
    long_description=open('README.md').read(),
    requires=['biopython'],
    ext_modules = [Extension('_wumanber',['amplishot/search/WuManber.cpp'])]
)
