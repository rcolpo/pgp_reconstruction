from setuptools import setup, find_packages


with open('README.md') as readme_file:
    readme = readme_file.read()

# Setting up
setup(
    name="pgp_reconstruction",
    version='0.0.1',
    author="rcolpo (Rodrigo Amarante Colpo)",
    author_email="<rodrigo.colpo-amarante@ufz.de>",
    description='Pathway-Guided Pruning Reconstruction of constraint-based metabolic models',
    long_description_content_type="text/markdown",
    long_description = readme,
    packages=find_packages(),
    include_package_data=True,
    install_requires=["reframed>=1.2.1", "pandas>=1.0", "cobra"],
	license="Apache Software License 2.0",
	keywords=['python', 'carveme', 'PGP Reconstruction', 'gapseq'],
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console', 
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: Apache Software License',
    ],
)