<h1 align="center">py2sterimol</h1> 
<h3 align="center">Thin client Python interface to the original Fortran implementation of Sterimol parameters.</h3>

<p align="center">  
  <img alt="py2sterimollogo" src="https://github.com/JacksonBurns/py2sterimol/blob/main/py2sterimol_logo.png">
</p> 
<p align="center">
  <img alt="GitHub Repo Stars" src="https://img.shields.io/github/stars/JacksonBurns/py2sterimol?style=social">
  <img alt="PyPI - Downloads" src="https://img.shields.io/pypi/dm/py2sterimol">
  <img alt="PyPI" src="https://img.shields.io/pypi/v/py2sterimol">
  <img alt="PyPI - License" src="https://img.shields.io/github/license/JacksonBurns/py2sterimol">
</p>

## ...what?
Sterimol parameters are of immense use for a number of subfields in chemistry. Unfortunately, the original implementation of the sterimol algorithm was written in free-form Fortran and left mostly untouched for a few decades. This prompted others to build new implementations using Python, see [wSterimol by Paton _et. al_](https://github.com/bobbypaton/wSterimol) and [morfeus by Jorner _et. al_](https://github.com/kjelljorner/morfeus/blob/main/morfeus/sterimol.py).

This is an excellent solution - Python is ubiquitous, there's no compilation required, and the source code is _readable_.

There are two issues which arise that probably make basically no difference to 99% of users: (1) algorithmic faith and (2) speed.
 1) Values generated with the original implementation and competing versions may yield slightly different (though likely chemically insignificant) results.
 2) Python is of course much slower than Fortran, so virtual screening may struggle.

Above all else, it's cool so why not?

## Compiling and Using Sterimol.f
Buckle up.

Because this file is _super old_ it can be difficult to compile. `gfortran` has not worked at all for me, regardless of specifying the standard, whether it was free or fixed form, etc. the [Intel Fortran Compiler](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran), however, seems to have no issues with the features and extensions that this program relies on.

It should look something like this: `ifort sterimol.f`

And you can run it like this: `type example.inp | .\sterimol.exe`

For unix-based systems, you need the same basic command for compilation, but then calling is a bit easier and actually follows the original documentation: `sterimol.exe < example.inp`

## Relevant Links
 - A like minded individual providing [the original code](http://www.ccl.net/cca/software/SOURCES/FORTRAN/STERIMOL/) as a matter of completionism.
 - Jackson and Paton's [sterimol](https://github.com/ipendlet/Sterimol), the first attempt at implementing sterimol in Python which conveniently includes some input files from which to infer the usage rules.

## License
Good lord, I have basically no idea.

The test input files are retrieved from sterimol, as referenced above. These are made available under the CC-BY license.

The original Fortran code - who knows? It was never released online because it predated the internet and the notion of open source altogether. I'm going to include it in this repository and use it to my heart's content, no modifications. As far as the license for this code - `py2sterimol` is licensed under the MIT license.

## sterimolarchive
All of the relevant source code files, tutorial documents, and whatever else I can find are included in the `sterimolarchive` directory.

## Online Documentation
[Click here to read the documentation](https://JacksonBurns.github.io/py2sterimol/)
