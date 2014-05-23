from distutils.core import setup, Extension

setup(
    name='Amplishot',
    version='0.9.0',
    author='Connor Skennerton',
    author_email='c.skennerton@gmail.com',
    packages=['amplishot', 'amplishot.test',
        'amplishot.parse', 'amplishot.app'],
    package_data={'amplishot': ['data/logWordPrior.txt']},
    scripts=['bin/amplishot'],
    url='http://pypi.python.org/pypi/Amplishot/',
    license='LICENSE.txt',
    description='Amplishot',
    long_description=open('README.md').read(),
    requires=['qiime', 'pyyaml'],
)
