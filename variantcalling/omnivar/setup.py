import setuptools
import os

here = os.path.abspath(os.path.dirname(__file__))
package_dir = os.path.join(here, 'omnivar')

def _get_version():
    versionpath = os.path.join(package_dir, 'VERSION')
    with open(versionpath) as versionfile:
        return versionfile.read().strip()

version = _get_version()

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="omnivar",
    version="0.0.1",
    author="Nathan Hammond",
    author_email="nathan.hammond@gmail.com",
    description="A local assembly-based variant caller",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://stanfordhealthcare.org/medical-clinics/clinical-genomics.html",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'pydot',
        'pyfaidx',
        'pysam',
    ],
    entry_points={
         'console_scripts': [
             'omnivar=omnivar.main:main',
         ],
     },
)
