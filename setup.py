#!/usr/bin/env python


from setuptools import setup, find_packages

# define distribution
setup(
    name = "snschops",
    version = "0.1",
    packages = find_packages(".", exclude=['setup.py', 'notebooks']),
    package_dir = {'': "."},
    # test_suite = 'tests',
    install_requires = [
        'numpy',
    ],
    dependency_links = [
    ],
    author = "ORNL/NDAV",
    description = "SNS chopper spectrometer utils",
    license = 'BSD',
    keywords = "neutron",
    url = "https://github.com/sns-chops/snschops",
    # download_url = '',
)

# End of file
