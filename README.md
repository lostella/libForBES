# ForBES

**ForBES** (standing for **For**ward-**B**ackward **E**nvelope **S**olver) is a MATLAB solver for
nonsmooth convex optimization problems.

It is generic in the sense that the user can customize the problem to solve in an easy and flexible way.
It is efficient since it features very efficient algorithms, suited for large scale applications.

For the documentation refer to the [ForBES webpage](http://lostella.github.io/ForBES/).

## Installation

Simply clone the git repository, or click on this [link](https://github.com/lostella/ForBES/archive/master.zip)
to download it as a zip archive and decompress the archive. Then move with the MATLAB command line to
the directory of ForBES, and execute the following command:

```
> setup_forbes
```

This will compile all the necessary source files and install the directory into MATLAB's path.

## How to use it

ForBES consists mainly of one MATLAB routine, `forbes`. In order to use it one
must provide a description of the problem and (optionally) a set of options:

```
out = forbes(f, g, init, aff, constr, opt);
```

Full documentation of ForBES, explaining how to specify these arguments, can be
found at the [ForBES webpage](http://lostella.github.io/ForBES/).

Examples on how to use `forbes` can be found in the [demos folder](https://github.com/lostella/ForBES/tree/master/demos).
Furthermore, you can access the help file of the solvers directly from MATLAB with

```
> help forbes
```

## Credits

ForBES is developed by Lorenzo Stella [`lorenzo.stella-at-imtlucca.it`] and Panos Patrinos [`panagiotis.patrinos-at-imtlucca.it`]
at [IMT Lucca](http://www.imtlucca.it).
Any feedback, bug report or suggestion for future improvements is more than welcome.
We recommend using the [issue tracker](https://github.com/lostella/ForBES/issues) to report bugs.
