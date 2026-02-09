# This code is part of kartograf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

"""
Tools for integration with miscellaneous non-required packages.
shamelessly borrowed from openff.toolkit
"""

import functools
from collections.abc import (
    Callable,
)


def requires_package(package_name: str) -> Callable:
    """Optional package detection function
    Helper function to denote that a funciton requires some optional
    dependency. A function decorated with this decorator will raise
    `MissingDependencyError` if the package is not found by
    `importlib.import_module()`.
    Parameters
    ----------
    package_name : str
        The directory path to enter within the context
    Raises
    ------
    MissingDependencyError
    """

    def test_import_for_require_package(function: Callable) -> Callable:
        @functools.wraps(function)
        def wrapper(*args, **kwargs):
            import importlib

            try:
                importlib.import_module(package_name)
            except (ImportError, ModuleNotFoundError):
                raise ImportError(f"{function.__name__} requires package: {package_name}")
            except Exception as e:
                raise e

            return function(*args, **kwargs)

        return wrapper

    return test_import_for_require_package
