from setuptools import setup, find_packages

setup(name='ccinput',
        version='1.0.0',
        description='Computational Chemistry Input Generator',
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
        ],
        test_suite='nose.collector',
        tests_require=[
            'nose',
            ],
        python_requires=">=3.6",
        zip_safe=False)
