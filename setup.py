from pathlib import Path
from setuptools import setup, find_packages
from versioneer import get_version, get_cmdclass

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(name='ccinput',
        version=get_version(),
        cmdclass=get_cmdclass(),
        description='Computational Chemistry Input Generator',
        long_description=long_description,
        long_description_content_type='text/markdown',
        url='http://github.com/cyllab/ccinput',
        author='Raphaël Robidas',
        author_email='Raphael.Robidas@USherbrooke.ca',
        license='BSD 3-Clause',
        classifiers=[
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: BSD License",
            "Operating System :: OS Independent",
            "Programming Language :: Python :: 3",
            "Topic :: Scientific/Engineering :: Chemistry",
        ],
        packages=find_packages(),
        entry_points={
            'console_scripts': [
                'ccinput = ccinput.wrapper:cmd',
            ],
        },
        install_requires=[
            'periodictable',
            'basis_set_exchange',
            'numpy',
            'versioneer',
        ],
        test_suite='nose.collector',
        tests_require=[
            'nose',
            ],
        python_requires=">=3.6",
        zip_safe=False)
