# coding utf8
import setuptools
from rrnafinder.versions import get_versions

with open('README.md') as f:
    LONG_DESCRIPTION = f.read()

setuptools.setup(
    name="rRNAFinder",
    version=get_versions(),
    author="Yuxing Xu",
    author_email="xuyuxing@mail.kib.ac.cn",
    description="A python tool for identify rRNA gene in plant genome.",
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    url="https://github.com/SouthernCD/rRNAFinder",
    include_package_data = True,

    entry_points={
        "console_scripts": ["rRNAFinder = rrnafinder.cli:main"]
    },

    packages=setuptools.find_packages(),

    install_requires=[
        "toolbiox>=0.0.40",
        "hugep2g>=1.0.1",
        "bcbio-gff>=0.6.6",
        "biopython<=1.80",
    ],

    python_requires='>=3.5',
)