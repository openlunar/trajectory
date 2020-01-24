# trajectory

Trajectory optimization algorithms for the earth--moon system.

* [Open Lunar Foundation](https://www.openlunar.org/)

## Description

This package hopes to become a repository for Python tools for
optimizing in the Earth--Moon system. Currently it includes only one
method, which is a planar patched conic strategy optimized via
sequential gradient restoration algorithm. A circular restricted
three-body problem solution is also in work.

## Requirements

* Python 3.x
* Numpy
* SciPy
* Matplotlib

## Installation

Clone this repository:

    git clone https://github.com/openlunar/trajectory.git

Currently, this tool cannot be installed as a package. You must run it
out of the repository directory. Forward work includes converting it
to a package with optimization functions that can be called from other
scripts.

## Usage

Most scripts are in the project root directory. These can be run from
a shell, e.g.

    ./patched_conic.py

## Developers

If you find a bug or wish to make a contribution, use the project's
[github issue tracker](https://github.com/openlunar/trajectory/issues).

## License

Copyright (c) 2019--2020, Open Lunar Foundation.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.