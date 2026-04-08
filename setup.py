from pathlib import Path

from setuptools import find_packages
from setuptools import setup


README = (Path(__file__).resolve().parent / "README.md").read_text(encoding="utf-8")

setup(
    name="ploidyscope",
    version="0.0.1",
    description="Windowed population-genetic summaries for mixed-ploidy data",
    long_description=README,
    long_description_content_type="text/markdown",
    packages=find_packages(include=["ploidyscope", "ploidyscope.*"]),
    include_package_data=True,
)
