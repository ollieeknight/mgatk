"""
mgatk: a mitochondrial genome analysis toolkit
"""
from setuptools import find_packages, setup

dependencies = [
    'click', 
    'pysam', 
    'pytest', 
    'snakemake', 
    'biopython', 
    'numpy', 
    'pandas', 
    'optparse-pretty', 
    'regex', 
    'ruamel.yaml',
    'pulp<=2.7.0',  # Add pulp with specific version constraint
    'matplotlib'    # Add matplotlib
]

setup(
    name='mgatk',
    version='0.7.0',
    url='https://github.com/caleblareau/mgatk',
    license='MIT',
    author='Caleb Lareau',
    author_email='caleb.lareau@gmail.com',
    description='Mitochondrial genome analysis toolkit.',
    long_description=__doc__,
    packages=find_packages(exclude=['tests']),
    include_package_data=True,
    zip_safe=False,
    platforms='any',
    install_requires=dependencies,
    entry_points={
        'console_scripts': [
            'mgatk = mgatk.cli:main',
            'mgatk-del-find = mgatk.deletioncalling.clifind:main',
            'mgatk-del = mgatk.deletioncalling.clidel:main'
        ],
    },
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]
)