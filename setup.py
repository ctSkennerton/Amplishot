from distutils.core import setup, Extension

setup(
    name='Amplishot',
    version='0.0.1',
    author='Connor Skennerton',
    author_email='c.tkennerton@gmail.com',
    packages=['amplishot', 'amplishot.test'],
    scripts=['bin/amplishot.py'],
    url='http://pypi.python.org/pypi/Amplishot/',
    license='LICENSE.txt',
    description='Amplishot',
    long_description=open('README.txt').read(),
    install_requires=[],
    ext_modules = [Extension('pywumanber',['amplishot/pywumanber/wumanber_impl.c'])]
)
