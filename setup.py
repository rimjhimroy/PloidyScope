from setuptools import find_packages
from setuptools import setup

setup(
    name="ploidyscope",
    version="0.0.1",
    description="Windowed population-genetic summaries for mixed-ploidy data",
    packages=find_packages(include=["ploidyscope", "ploidyscope.*"]),
    include_package_data=True,
)
