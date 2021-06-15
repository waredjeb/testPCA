#!/usr/bin/python3
import matplotlib.pyplot as plt
import pickle as pl
import numpy as np
import sys
import mplhep as hep
plt.style.use([hep.style.ROOT, hep.style.firamath])

hep.rcParams.label.data = True
hep.rcParams.label.paper = False

# print (sys.argv)

# Load figure from disk and display
figurename = sys.argv[1]
print(figurename)
fig_handle = pl.load(open(figurename,'rb'))
fig_handle.show()
a = input("A")