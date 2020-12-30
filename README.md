# Path phasing

Author: Masutani Bansho
Email: ban-m@g.ecc.u-tokyo.ac.jp


# TL;DR

This repository contains 1. implementation of path phasing algorithm and 2. it's command line interface.

To use CLI, you only need a utf-8 encoded tsv file FILE, which contents RECORDS describes paths as follows:

RECORDS = (<RECORD>\n)^*

RECORD = (ID)\t(NODE:CLUSTER\t)^+

Here, ID is a string(,or integer) without tab character, and NODE and CLUSTER are integer.
For example,

RECORD = Test Read\t10:0\t11:0\t2341:13213\t343:324

As you can see, the number of node shouldn't be consective, so, feel easy. The graph containing all of reads would be created on the fly.

Output is tsv, consisting RECORD\tCUSTER. CLUSTER is either 0 or 1.


If you supply dataset consisting of two connected component, then these information would be LOST.
Also, when the occurence of a node is high(>15), we drop some of the reads and predict their phase afterwards.
So I DON'T advocate that this CLI solves the phasing problem exactly or output optimal solution.


# Contents


## What's path phasing?


It is phasing on a graph, not a "reference coordinate".