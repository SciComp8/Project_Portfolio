from setuptools import setup, find_packages

setup(
    author="Anni Liu",
    description="A package for predicting perturbation response in the single cell RNA space",
    name="predict_perturbation",
    packages=find_packages(include=["predict_perturbation", "predict_perturbation.*"]),
    version="0.1.0",
    install_requires=['adjustText', 'matplotlib', 'numpy', 'anndata', 'pertpy', 'scanpy', 'scgen', 'scipy']
)
