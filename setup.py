#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-03-18

from pathlib import Path
import pickle
import sys
from pkg_resources import (
    normalize_path,
    working_set,
    add_activation_listener,
    require,
)
try:
    from setuptools import setup, find_packages
    from setuptools.command.build_py import build_py
    from setuptools.command.develop import develop
    from setuptools.command.test import test as test_command
except ImportError:
    raise ImportError("Please install setuptools.")

SETUP_DIR = Path(__file__).parent.resolve()
MODELS_DIR = SETUP_DIR / "stan"
MODELS_TARGET_DIR = Path("sphere/stan_models")


def take_package_name(name):
    if name.startswith("-e"):
        return name[name.find("=")+1:name.rfind("-")]
    else:
        return name.strip()


def load_requires_from_file(filepath):
    with open(filepath) as fp:
        return [take_package_name(pkg_name) for pkg_name in fp.readlines()]


def compile_stan_models(target_dir, models_dir=MODELS_DIR):
    from pystan import StanModel
    model_path_list = MODELS_DIR.glob("*.stan")
    for model_path in model_path_list:
        model_type = model_path.stem
        target_name = model_type + ".pkl"
        target_path = target_dir / target_name
        with open(model_path) as f:
            model_code = f.read()
        model = StanModel(model_code=model_code, model_name=model_type)
        with open(target_path, "wb") as f:
            pickle.dump(model, f, protocol=pickle.HIGHEST_PROTOCOL)


class BuildPyCommand(build_py):
    """Custom build command to pre-compile Stan models."""
    def run(self):
        if not self.dry_run:
            build_lib_path = Path(self.build_lib).resolve()
            target_dir = build_lib_path / MODELS_TARGET_DIR
            self.mkpath(str(target_dir))
            compile_stan_models(target_dir)
        build_py.run(self)


class DevelopCommand(develop):
    """Custom develop command to pre-compile Stan models in-place."""

    def run(self):
        if not self.dry_run:
            build_lib_path = Path(self.setup_path).resolve()
            target_dir = build_lib_path / MODELS_TARGET_DIR
            self.mkpath(str(target_dir))
            compile_stan_models(target_dir)
        develop.run(self)


class TestCommand(test_command):
    """We must run tests on the build directory, not source."""

    def with_project_on_sys_path(self, func):
        # Ensure metadata is up-to-date
        self.reinitialize_command('build_py', inplace=0)
        self.run_command('build_py')
        bpy_cmd = self.get_finalized_command("build_py")
        build_path = normalize_path(bpy_cmd.build_lib)

        # Build extensions
        self.reinitialize_command('egg_info', egg_base=build_path)
        self.run_command('egg_info')

        self.reinitialize_command('build_ext', inplace=0)
        self.run_command('build_ext')

        ei_cmd = self.get_finalized_command("egg_info")

        old_path = sys.path[:]
        old_modules = sys.modules.copy()

        try:
            sys.path.insert(0, normalize_path(ei_cmd.egg_base))
            working_set.__init__()
            add_activation_listener(lambda dist: dist.activate())
            require('%s==%s' % (ei_cmd.egg_name, ei_cmd.egg_version))
            func()
        finally:
            sys.path[:] = old_path
            sys.modules.clear()
            sys.modules.update(old_modules)
            working_set.__init__()


setup(
    name="sphere",
    version='1.0.0',
    packages=find_packages(),

    install_requires=load_requires_from_file("requirements.txt"),

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
            "sphere_dplot = sphere.sphere_dplot:main_wrapper",
            "sphere_mcplot = sphere.sphere_mcplot:main_wrapper",
            "sphere_cstats= sphere.sphere_cstats:main_wrapper",
            "sphere_filter= sphere.sphere_filter:main_wrapper",
            "sphere_symtest= sphere.sphere_symtest:main_wrapper",
            "sphere_rdump= sphere.sphere_rdump:main_wrapper"
        ]
    },
    cmdclass={
        "build_py": BuildPyCommand,
        "develop": DevelopCommand,
        "test": TestCommand
    },
    test_suite="sphere.tests",
    include_package_data=True,
    zip_safe=False
)
