# Square Root Information Filter

## Files and Folders
- ./book: translated Chinese version of the book "Factorization Methods for Discrete Sequential Estimation"
- ./doc: documentation of the ESL
- ./esl: the Estimation Subroutine Library (ESL)
- ./test: a Visual Studio project for testing the ESL
- ./gitignore
- ./README.md: the current file
- ./SRIF.sln: the Visual Studio solution that hosts the test project of the ESL

## Estimation Subroutine Library (ESL)

This Estimation Subroutine Library (ESL) provides implementations of the
algorithms described in the book "Factorization Methods for Discrete Sequential
Estimation" (see below). Primarily UDU**T and SRIF based algorithms for solving
"Kalman Filter" type problems ... both for estimation and for smoothing.

Typical applications range from navigation (e.g. GPS), radar, computer vision,
and even advanced sprinkler controls.

Keywords: UDU, kalman filter, covariance filter, square root information
filter, SRIF, navigation, smoothing

To paint with a broad and simplistic brush, if one is memory constrained or
needs the updates one point at a time, the UDU**T formulation family of
routines is probably preferable. If one's problem is naturally constructed as a
series of batches, the SRIF formulation has some advantages (in addition,
there's a fairly obvious way of exploiting parallel processing for the SRIF
case; but it's not coded in the library per se, but one can use the routines on
multiple parts of a cluster in the obvious fashion ;>).

There's a Creative Commons license and some historical overview (from
20+ years ago) in the file brtsgn.f.


---
Keith H. Bierman   khbkhb@gmail.com
5430 Nassau Circle East
Cherry Hills Village, CO 80113

##