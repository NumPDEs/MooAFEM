# Matlab Object Oriented FEM

```
 ________________________________________
/                                        \
| Welcome to MooAFEM!                    |
|                                        |
| A Matlab Object Oriented AFEM package. |
\                                        /
 ----------------------------------------
        \   ^__^
         \  (oo)\_______
            (__)\       )\/\
                ||----w |
                ||     ||
```
The ASCII art is taken from the [cowsay](https://github.com/tnalpgge/rank-amateur-cowsay) package by Tony Monroe.

## Short description

An object oriented Matlab library for (adaptive) finite element analysis.
Advantages:
- Easy to use and modify
- Covers a wide range of equations
- Efficient implementation of general polynomial orders
- A lot of built-in convenience functions

## Installation

To get started, run `setup.m` in the root folder. This adds everything to the
path and compiles .mex files if necessary.
Note that, at least, Matlab version R2020b is required.

## First steps

There are some examples that show the features of this package in the
`examples/` directory. These are simple scripts, so they can be executed either
from the command window or by running it from the editor.

## Modifications

To allow for reliable modifications, a comprehensive test suite is shipped
with the package. These tests use the matlab.unittest framework and can be run
by one of the following commands:
```
% automatically run all unit-tests
runtests('tests/')

% automatically run all unit-tests and all examples
runtests('tests/', 'IncludeSubfolders', true)

% run tests and perform code coverage analysis
runtests('tests/', ['IncludeSubfolders', true,] 'ReportCoverageFor', 'lib/')
