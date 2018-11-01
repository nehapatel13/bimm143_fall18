Untitled
================
Neha Patel
10/24/2018

OUR FUNCTION
============

Kinase with drug
================

library(bio3d) x &lt;- "4AKE" pdb &lt;- read.pdb(x) chainA &lt;- trim.pdb(pdb, chain = "A", elety = "CA") sb1 &lt;- chainA*a**t**o**m*b plotb3(sb1, sse=chainA, typ="l", ylab="Bfactor")

Kinase without Drug
===================

library(bio3d) x &lt;- "1AKE" pdb &lt;- read.pdb(x) chainA &lt;- trim.pdb(pdb, chain = "A", elety = "CA") sb2 &lt;- chainA*a**t**o**m*b plotb3(sb2, sse=chainA, typ="l", ylab="Bfactor")

Kinase with Drug
================

library(bio3d) x &lt;- "1E4Y" pdb &lt;- read.pdb(x) chainA &lt;- trim.pdb(pdb, chain = "A", elety = "CA") sb3 &lt;- chainA*a**t**o**m*b plotb3(sb3, sse=chainA, typ="l", ylab="Bfactor")

Plot all together
=================

FUNCTION FOR ASSIGNMENT
=======================

library(bio3d) my\_function &lt;- function(x) { pdb &lt;- read.pdb(x) chainA &lt;- trim.pdb(pdb, chain = "A", elety = "CA") b &lt;- chainA*a**t**o**m*b plotb3(b, sse=chainA, typ="l", ylab="Bfactor") }
