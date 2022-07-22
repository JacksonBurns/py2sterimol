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

## Online Documentation
[Click here to read the documentation](https://JacksonBurns.github.io/py2sterimol/)
