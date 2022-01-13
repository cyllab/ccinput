from setuptools import setup, find_packages

setup(name='ccinput',
        version='0.1',
        description='Computational Chemistry Input Generator',
        url='http://github.com/cyllab/ccinput',
        author='Raphaël Robidas',
        author_email='Raphael.Robidas@USherbrooke.ca',
        license='BSD 3-Clause',
        packages=find_packages(),
        install_requires=[
            'periodictable',
            'basis_set_exchange',
            'numpy',
        ],
        test_suite='nose.collector',
        tests_require=[
            'nose',
            ],
        zip_safe=False)
