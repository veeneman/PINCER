#!/usr/bin/env python2.7
import os
import sys
import azimuth.model_comparison
import numpy as np
model = os.path.join(os.path.dirname(azimuth.__file__), 'saved_models', 'V3_model_nopos.pickle')

g = np.loadtxt(sys.argv[1],dtype = [('id', object),('seq','S30')])
az = azimuth.model_comparison.predict(g['seq'], None, None, None, model)
if np.size(g,0) != np.size(az):
 sys.exit("Error: Efficacy score not calculated for all guides")
np.savetxt(sys.argv[2], np.c_[g['id'], az], fmt = ['%s', '%.18e'], delimiter = "\t")
