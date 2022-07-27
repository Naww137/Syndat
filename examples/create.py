import numpy as np
from matplotlib.pyplot import *
import syndat
import pandas as pd
import os

sammy_directory =  os.path.realpath('../synthetic_data/Ta181')

opendata = os.path.join(sammy_directory,'rpi-open-ta181.csv')
sammy_xs = os.path.join(sammy_directory,'SAMMY.LST')


# initialize true underylying reduction parameters
gen = syndat.generation(True, opendata, sammy_xs)

red = syndat.reduction(gen)


#%%
figure()

plot(red.sdat.tof,red.sdat.c, label='cts in'); 
plot(gen.odat.tof,gen.odat.c, label='cts out')
plot(gen.odat.tof, gen.redpar.val.ks*gen.Bi+gen.redpar.val.b0s, label='bkg_i')
plot(gen.odat.tof, gen.redpar.val.ko*gen.Bi+gen.redpar.val.b0o, label='bkg_o')
xscale('log'); yscale('log')
legend()

