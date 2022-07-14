
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import syndat
import scipy.stats as stat
import matplotlib as mpl





def sample_poisson(m):
    return np.random.default_rng().poisson(lam=m)

def sample_gaussian(m):
    return np.random.default_rng().normal(loc=m,scale=np.sqrt(m))

m=1e3

test = sample_gaussian(m)
print(test)
print(round(test))


s = 1e4

for i in s