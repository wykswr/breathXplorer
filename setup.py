from setuptools import setup, find_packages

with open("./requirements.txt") as f:
    required_packages = f.read().splitlines()

setup(
    name="breathXplorer",
    version="0.1.0",
    packages=find_packages(),
    install_requires=required_packages,
    author="wykswr",
    author_email="bifocal.above.0y@icloud.com",
    description="A short description of your package",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/wykswr/breathXplorer-lite",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)
