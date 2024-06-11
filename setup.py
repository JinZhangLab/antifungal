from setuptools import setup, find_packages

setup(
    name='antifungal',
    version='0.1.3',
    packages=find_packages(), 
    include_package_data=True, 
    author='Jin Zhang',
    author_email='jzhang@chemoinfolab.com',
    description='A python package for activity prediction and in rational design of antifungal peitdes.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/JinZhangLab/antifungal',
    install_requires=[
        'numpy>=1.13.3',
        'pandas>=1.2.4',
        'modlamp>=4.3.0',
        'scipy<=1.7.0',
        'scikit-learn==0.24.2',
        'tqdm',
        'matplotlib',
    ],
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    python_requires = '>=3.6, <3.9',
)