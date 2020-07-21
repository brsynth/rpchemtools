from setuptools import setup
import versioneer

requirements = [
    'rdkit==2018.09.1.0'
]

setup(
    name='rpchemtools',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='Minimalist toolbox to deal with chemicals',
    license='MIT',
    author='Thomas Duigou, Baudoin Delépine',
    author_email='thomas.duigou@inrae.fr',
    url='https://github.com/tduigou/rpchemtools',
    packages=['rpchemtools'],
    install_requires=requirements,
    keywords='rpchemtools',
    classifiers=[
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ]
)
