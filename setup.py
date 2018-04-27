#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-03-18

try:
    from setuptools import setup, find_packages
except ImportError:
    raise ImportError("Please install setuptools.")

setup(
    name="sphere",
    version='1.0.0',
    packages=find_packages(exclude=[
        "tests"
    ]),

    install_requires=[
        "numpy>=1.13.1",
        "pandas>=0.20.3",
        "matplotlib>=2.0.2",
        "pystan>=2.16.0.0",
        "scipy>=1.0.0"
    ],
    extras_require={
        "test": ["nose"]
    },
    test_suite="nosetests",

    author="Shinya SUZUKI",
    author_email="sshinya@bio.titech.ac.jp",
    description='Synthetic PHasE Rate Estimator by single metagenome sequence',
    long_description="""
        Synthetic PHasE Rate Estimator by single metagenome sequence
    """,
    license="BSD 3-Clause License",
    keywords=["Bioinformatics", "Metagenome", "Microbiome", "Data analysis"],
    url="https://github.com/TaskeHAMANO/SPHERE",
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Infomation Analysis",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: BSD License"
    ],
    entry_points={
        "console_scripts": [
            "sphere_estimate = sphere.sphere_estimate:main_wrapper",
            "sphere_compile = sphere.sphere_compile:main_wrapper",
            "sphere_dplot = sphere.sphere_dplot:main_wrapper",
            "sphere_mcplot = sphere.sphere_mcplot:main_wrapper",
            "sphere_cstats= sphere.sphere_cstats:main_wrapper",
            "sphere_filter= sphere.sphere_filter:main_wrapper",
            "sphere_symtest= sphere.sphere_symtest:main_wrapper"
        ]
    }
)
