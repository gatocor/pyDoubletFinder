import sys

try:
    from skbuild import setup
except ImportError:
    print('Please update pip, you need pip 10 or greater,\n'
          ' or you need to install the PEP 518 requirements in pyproject.toml yourself', file=sys.stderr)
    raise
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="pyDoubletFinder",
    version="0.1",
    description="Reimplementation in ppython of DoubletFinder in R.",
    author='Gabriel Torregrosa Cortes',
    author_email="g.torregrosa@outlook.com",
    url="https://github.com/gatocor/pyDoubletFinder",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    packages=['pyDoubletFinder'],
    package_dir={'': 'src'},
    cmake_install_dir='src/scBayesDeconv',
    python_requires = ">=3.5",
    install_requires = ["numpy>=1.17.5","scipy>=1.5.0","dynesty>=0.9.7","pandas>=0.25.0"]
)