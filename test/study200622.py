#!/usr/bin/python3

import os
os.system("build/salpa -i /data/dw/tmp/littlechunk.dat -P /data/dw/tmp/little_chunk_events.out -o /data/dw/littlesalpa.dat -c 1 -C 384 -b 1 -Z -x 6 -f 5.1")

import numpy as np
with open('/data/dw/littlesalpa.dat', mode='rb') as f:
    dat = f.read()
    dat = np.frombuffer(dat, dtype=np.int16)
    L = len(dat)
    dat = np.reshape(dat, [L//384, 384])
    T, C = dat.shape
    
import pyqplot as qp

tt = np.arange(T) / 30000.0

qp.figure('/tmp/s1')
i0 = 28000
i1 = 32000
qp.plot(tt[i0:i1], dat[i0:i1,0])
