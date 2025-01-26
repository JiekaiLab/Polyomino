from setuptools import setup
import setuptools
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.rst").read_text()
setup(
    name='Zmap-sc',
    version='1.0.0',
    description='An algorithm framework employing multi-layered regional constraints to accurately assign cell locations, enhancing spatial accuracy and resilience to noise.',
    url='https://github.com/caiquanyou/Zmap',
    author='Cai Quanyou',
    author_email='cai_quanyou@gibh.ac.cn',
    python_requires='>=3.8.0',
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=[
        "lap>=0.4.0",
        "scikit-learn>=1.1.3",
        "numpydoc>=1.1",
        "anndata>=0.8.0",
        "scanpy>=1.9.1",
        "scipy>=1.9.3",
        "numba>=0.14.0",
        "tqdm>=4.66.5",
        "igraph>=0.11.6",
        "leidenalg>=0.10.2"
    ],
    entry_points={
        'console_scripts': [
            'zmap=Zmap.Zmap:main_Zmap',
        ]},
    classifiers=["Programming Language :: Python :: 3.8", "Operating System :: MacOS",],
    include_package_data=False
)