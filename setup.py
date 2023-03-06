from setuptools import setup

setup(
    name="GTFProcessing",
    version="0.0.1",
    py_modules=["GTFProcessing", "guide_table_processing", "PerturbLoader_generation", "preprocessing", "runSceptre", "test_multiple_files"],
    install_requires=["pandas", "gtfparse", "tqdm"]
)
