from setuptools import setup, find_packages

with open("requirement.txt") as f:
    required_packages = f.read().splitlines()

setup(
    name="breathXplorer",
    version="0.1.4",
    packages=find_packages(),
    install_requires=required_packages,
    author="wykswr",
    author_email="bifocal.above.0y@icloud.com",
    description="A toolkit for breath metabolomics analysis",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/wykswr/breathXplorer",
    license="MIT",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
)
