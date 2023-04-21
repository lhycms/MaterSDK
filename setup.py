from setuptools import setup, find_packages


setup(
    name="matersdk",
    version="v1.0",
    author="Liu Hanyu  &&  LONXUN QUANTUM",
    author_email="domainofbuaa@gmail.com",
    url="https://github.com/lhycms/MaterSDK",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
            "pymatgen>=2022.11.7",
            "numpy>=1.23.5",
            "prettytable>=3.5.0",
            "dpdata>=0.2.13",
            "click>=8.1.3",
            "joblib>=1.2.0",
    ]
)