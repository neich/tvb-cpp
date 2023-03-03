import sys

try:
    from skbuild import setup
except ImportError:
    print(
        "Please update pip, you need pip 10 or greater,\n"
        " or you need to install the PEP 518 requirements in pyproject.toml yourself",
        file=sys.stderr,
    )
    raise

from setuptools import find_packages

setup(
    name="tvbcpp",
    version="0.0.3",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    cmake_install_dir="src/tvbcpp",
    include_package_data=True,
    python_requires=">=3.6",
)
