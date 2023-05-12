t8wyo
=========
t8wyo is an example software framework built on top of t8code.

---
# 1. Obtaining the code
    git clone https://github.com/dg4est/t8wyo.git

---
# 2. Compilation
    To compile t8wyo and the 3rd party libraries (t8code, metis):
    ./build.sh <options>

**Usage**: ./build.sh `{OPTIONS}`...`{COMPILER OPTIONS}`...`{3PL OPTIONS}`

| OPTIONS:             | Shortcut    | Description                                         |
|:---------------------|-------------|:----------------------------------------------------|
| --3pl                | -3pl        | build the 3rd party libraries: metis, t8code        |
| --t8wyo              | -t8wyo      | build the t8wyo solver                              |
|                      |             |                                                     |
| --help               | -h          | displays this help message                          |
| --clean              | -c          | removes local build directories                     |
| --distclean          | -dc         | removes builds and install directories              |
| --release            | -opt        | compile the project in optimized mode               |
| --debug              | -deb        | compile the project in debug mode                   |
|                      |             |                                                     |
| **COMPILER OPTIONS**:|             |                                                     |
| CC=`<arg>`           | cc=`<arg>`  | set the C compiler                                  |
| CXX=`<arg>`          | cxx=`<arg>` | set the C++ compiler                                |
| FC=`<arg>`           | fc=`<arg>`  | set the Fortran compiler                            |
|                      |             |                                                     |
| **3PL OPTIONS**:     |             |                                                     |
| --ALL3P              | -all3p      | compile all 3rd party libraries                     |
| --T8CODE             | -t8code     | compile t8code                                      |
| --METIS              | -metis      | compile metis                                       |

# Common Build Options:
**Default (-go)**: Sets CC=mpicc CXX=mpicxx FC=mpif90
    ./build.sh -go

**Intel MPI (-impi)**: Sets CC=mpiicc CXX=mpiicpc FC=mpiifort
    ./build.sh -impi

---
# 3. Dependencies
    Software:
        CMake (Version 2.8 or newer): Also recommend ccmake or cmake-gui.
        MPI: optional MPI-3 one-sided functions available
        libtool: apt install libtool
        automake: apt install automake
        zlib: apt install zlib1g zlibc

---
# 4. License

The MIT License (MIT)

Copyright (c) 2023 Andrew C. Kirby, Dimitri J. Mavriplis

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
