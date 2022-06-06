import setuptools

setuptools.setup(
    name="smmo",
    version="0.0.1",
    license='MIT',
    author="Kwangnam Yu",
    author_email="phykn.kr@gmail.com",
    description="Scattering-Matrix method for multilayer optics",
    long_description=open("README.md").read(),
    url="https://github.com/phykn/smmo",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ],
)