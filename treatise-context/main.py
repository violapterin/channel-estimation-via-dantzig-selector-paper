#! /usr/bin/env python3

import os

import functions as fct
import classes as cls

# example
NUM_REP_FIG = 6 # example
for i in range (NUM_REP_FIG):
    print (
        "~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\n",
        "Figure", i,
        "~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\n"
    )
    idx_hold = i + 1
    # example
    ver = cls.Version (
        cls.Focus.ASSORTED,
        cls.Size.BIG,
        cls.Ratio.WIDE,
        str (idx_hold))

    fct.execute (ver)
