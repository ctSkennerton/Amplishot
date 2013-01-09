from distutils.core import setup, Extension

setup(
    name='Amplishot',
    version='0.0.1',
    author='Connor Skennerton',
    author_email='c.tkennerton@gmail.com',
    packages=['amplishot', 'amplishot.test', 'amplishot.pywumanber'],
    scripts=['bin/amplishot'],
    url='http://pypi.python.org/pypi/Amplishot/',
    license='LICENSE.txt',
    description='Amplishot',
    long_description=open('README.md').read(),
    install_requires=['cogent >= 1.5.1'],
    ext_modules = [Extension('pywumanber',['amplishot/pywumanber/wumanber_impl.c'])]
)
