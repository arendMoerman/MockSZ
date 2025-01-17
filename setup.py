import os
import pathlib

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext as build_ext_orig

__version__ = "0.2.4"

class CMakeExtension(Extension):

    def __init__(self, name):
        # don't invoke the original build_ext for this special extension
        super().__init__(name, sources=[])


class build_ext(build_ext_orig):

    def run(self):
        for ext in self.extensions:
            self.build_cmake(ext)
        #super().run()

    def build_cmake(self, ext):
        cwd = pathlib.Path().absolute()

        build_temp = pathlib.Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)
        extdir = pathlib.Path(self.get_ext_fullpath(ext.name))

        config = 'Release ' + '--j4'
        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + str(extdir.parent.absolute()),
            '-DCMAKE_BUILD_TYPE=' + config
        ]

        build_args = [
            '--config', config
        ]

        os.chdir(str(build_temp))
        self.spawn(['cmake', str(cwd)] + cmake_args)
        if not self.dry_run:
            if os.name != "nt":
                self.spawn(['cmake', '--build', '.'] + build_args)
            else:
                self.spawn(['cmake', '--build', '.'])
        os.chdir(str(cwd))

setup(
    name='MockSZ',
    license="MIT",
    version=__version__,
    author="Arend Moerman",
    install_requires = ["numpy", "wheel", "cmake", "scipy"],
    package_dir = {'': 'src'},
    packages=['MockSZ'],
    ext_modules=[CMakeExtension(os.path.join("MockSZ", "libs"))],
    cmdclass={'build_ext': build_ext},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: C",
        "Programming Language :: C++",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8'
)

