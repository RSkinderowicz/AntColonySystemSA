# AntColonySystemSA

This is an implementation of the Ant Colony System combined with Simulated
Annealing algorithm for solving the Sequential Ordering Problem.
The implementation was used in computational experiments presented in the
article:

_Rafał Skinderowicz, "An improved Ant Colony System for the Sequential Ordering
Problem", Computers & Operations Research, DOI: 10.1016/j.cor.2017.04.012_

# Building

This software is intended to compile and run on Linux.
It was tested with GCC v5.4 but should work with an earlier version, as long as
it supports C++11

To compile & build run:

    make

It could take a while. If everything goes OK, "acs" executable should be
created in the current directory.

# Running

The framework takes a lot of arguments, most of which have some default values.
To see the default help page run:
    ./acs --help

An examplary run of the Enhanced ACS (EACS) algorithm with the
SOP-3-exchange-SA local search heuristic for the R.700.100.15 SOP instance from
the SOPLIB06 repository can be initiated with:

    ./acs --alg=eacs --ants=10 --phmem=std --ls=sop3exchange-sa --beta=2.0
    --test=R.700.100.15.sop
    --q0=10 --phi=0.01 --rho=0.1
    --sa-cooling-ratio=0.9999 --sa-init-accept-prob=0.1
    --sop-ls-sa-cooling-ratio=0.99 --sop-ls-sa-init-accept-prob=0.1
    --timeout=10 --trials=1 

where --timeout=10 sets the timelimit to 10 seconds.

After approx. 10 sec. the program should terminate and the results should be
saved to a \*.js file in "results/" folder. The results are saved in JSON
format which can be easily parsed in almost any programming language.

# License

The source code is licensed under the [MIT
License](http://opensource.org/licenses/MIT):

Copyright 2017 Rafał Skinderowicz (rafal.skinderowicz@us.edu.pl)

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
