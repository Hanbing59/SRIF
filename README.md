# SRIF

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


ISBN 720 876 1333
http://store.doverpublications.com/0486449815.htm
http://www.amazon.com/Factorization-Methods-Discrete-Sequential-Estimation/dp/0486449815/ref=pd_bbs_sr_1?ie=UTF8&s=books&qid=1204522450&sr=8-1
http://search.barnesandnoble.com/Factorization-Methods-for-Discrete-Sequential-Estimation/Gerald-J-Bierman/e/9780486449814/?itm=1


--
Keith H. Bierman   khbkhb@gmail.com
5430 Nassau Circle East
Cherry Hills Village, CO 80113