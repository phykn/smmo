import setuptools


with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="smmo",
    version="0.0.4",
    author="Kwangnam Yu",
    author_email="phykn.kr@gmail.com",
    description="Scattering-Matrix method for multilayer optics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/phykn/smmo",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        "numpy"
    ]
)