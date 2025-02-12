
#***********************************************************************
#DD package, data collection and analysis of 2D electronic spectra
#Copyright (C) 2016, 2017  Jan Alster (Charles Univesity, Prague)
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>
#***********************************************************************

# -*- coding: utf-8 -*-
"""
Example of how to make use of the python wrappers
"""

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import libPyPopot
import math
import time

#Seed initialization
libPyPopot.seed((int)((time.time()%100)*10000));

#Set the dimension of our problem
dimension = 10

def lbound(index):
    return -30

def ubound(index):
    return 30

def stop(fitness, epoch):
    return (fitness <= 100) | (epoch >= 10000)
    
def evaluate(param):
    fit = 0.0
    for i in range(len(param) - 1):
        fit += 100 * (param[i+1] - param[i]**2)**2 + (param[i] - 1.0)**2
    return fit

#Initialize libPyPopot
algo = libPyPopot.Stochastic_PSO_2006(dimension, lbound, ubound, stop, evaluate)

#Run until the stop criterion is met
algo.run(0)
print("Stop after {epoch} steps".format(epoch = algo.getEpoch()))
print("Best fitness : {fitness}".format(fitness = algo.bestFitness()))
print("Best solution : {sol}".format(sol = algo.bestParticule()))

# #Restart the PSO
# algo = libPyPopot.Stochastic_PSO_2006(dimension, lbound, ubound, stop, evaluate)

# #A step by step optimization
# for i in range(100):
#     algo.step()
#     print("Epoch {epoch} : best fitness {fitness}".format(epoch = algo.getEpoch(
# ), fitness = algo.bestFitness()))
