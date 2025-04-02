import io
import os
import sys
from shutil import rmtree
from setuptools import find_packages, setup, Command

# Package meta-data.
NAME = 'prism_grn'
DESCRIPTION = 'A Probabilistic Model for Recovering GRNs based on Multi-Omics Data'
URL = 'https://github.com/Ying-Lab/PRISM/'
EMAIL = 'zscotty@stu.xmu.edu.cn'
AUTHOR = 'Wenhao Zhang'
REQUIRES_PYTHON = '>=3.8.10'
VERSION = '1.2.0'

# Required packages
REQUIRED = [
    'torch==1.11.0',
    'pyro-ppl==1.7.0',
    'networkx==3.1',
    'numpy==1.24.4',
    'scipy==1.10.1',
    'pandas==2.0.3'
]

# Optional dependencies
EXTRAS = {
    'dev': ['pytest', 'black', 'flake8'],
}

here = os.path.abspath(os.path.dirname(__file__))

# Read README.md for long_description
try:
    with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

# Version handling
about = {'__version__': VERSION} if VERSION else {}
if not VERSION:
    with open(os.path.join(here, NAME.replace("-", "_"), '__version__.py')) as f:
        exec(f.read(), about)

class UploadCommand(Command):
    """Support setup.py upload."""
    description = 'Build and publish the package.'
    user_options = []

    @staticmethod
    def status(s):
        """Prints messages in bold."""
        print(f'\033[1m{s}\033[0m')

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            self.status('Removing previous builds…')
            rmtree(os.path.join(here, 'dist'))
        except OSError:
            pass

        self.status('Building Source and Wheel distribution…')
        os.system(f'{sys.executable} setup.py sdist bdist_wheel')

        self.status('Uploading the package to PyPI via Twine…')
        os.system('twine upload dist/*')

        self.status('Pushing git tags…')
        os.system(f'git tag v{about["__version__"]}')
        os.system('git push --tags')

        sys.exit()

# Package setup
setup(
    name=NAME,
    version=about['__version__'],
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=find_packages(include=['prism', 'prism.*']),
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    include_package_data=True,
    license="MIT",
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
    cmdclass={
        'upload': UploadCommand,
    },
)