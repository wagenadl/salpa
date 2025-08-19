#!/usr/bin/python3

import salpa_cppcore
import numpy as np
import os

class Salpa:
    def __init__(self, data, tau, rail1=-np.inf, rail2=np.inf, thresh=np.inf,
                 t_blankdepeg=5, t_ahead=5, t_chi2=15):
        self.data = data.astype(np.float32)
        N = len(data)
        self.out = np.zeros(N, np.float32)
        self.csalpa = salpa_cppcore.csalpa(self.data,
                                           self.out,
                                           thresh, tau)
        self.csalpa.set_t_blankdepeg(t_blankdepeg)
        self.csalpa.set_t_ahead(t_ahead)
        self.csalpa.set_t_chi2(t_chi2)
        self.csalpa.setrail(rail1, rail2)
    def complete(self):
        self.partial(len(self.data))
        return self.out
    def partial(self, t_to):
        self.csalpa.process(t_to)
    def forcepeg(self, t_from, t_to):
        self.csalpa.forcepeg(t_from, t_to)
    def result(self):
        return self.out

    
def salpa(data, tau, rail1=-np.inf, rail2=np.inf, thresh=np.inf,
          t_blankdepeg=5, t_ahead=5, t_chi2=15, tt_stimuli=None, t_forcepeg=10):
    '''SALPA - Suppresion of Artifacts by Local Polynomial Approximation
    y = SALPA(x, tau) performs SALPA on the data X, which must be a 1D
    numpy array. Note that THRESH is absolute. You need to estimate the noise
    level with an external tool.
    See Wagenaar and Potter (2001).'''
    slp = Salpa(data, tau, rail1, rail2, thresh, t_blankdepeg, t_ahead, t_chi2)
    if tt_stimuli is not None:
        for t in tt_stimuli:
            slp.forcepeg(t, t+t_forcepeg)
    return slp.complete()


if __name__ == '__main__':
    import qplot as qp
    ## Prepare some fake data
    tt = 4*np.pi*np.arange(1e4)/1e4
    xx = np.sin(tt) + np.random.randn(len(tt))#sin(30*tt)
    A=100
    xx[np.logical_and(tt>4, tt<4.5)] = A
    iart = tt>=4.5
    tart = tt[iart] - 4.5
    xx[iart] += A*np.cos(np.sqrt(tart/.004))*np.exp(-tart/1)

    def show(tt, xx, yy, name):
        qp.figure(f'/tmp/{name}',7,2)
        qp.pen('469', 1, join='round')
        qp.plot(tt, xx)
        qp.pen('r')
        qp.plot(tt,yy)
        qp.shrink()
    
    ## Demonstrate full processing
    yy = salpa(xx, tau=75, thresh=3, rail1=-A*.99, rail2=A*.99)
    print("Hello")
    print(yy.shape)
    show(tt, xx, yy, 'salpa')
    
    ## Demonstrate without rails
    yy = salpa(xx, tau=75, thresh=3)
    show(tt, xx, yy, 'salpa0')
    
    ## Demonstrate partial processing
    print(xx[0], xx[3173])
    slp = Salpa(xx, tau=75, thresh=3)
    istart = np.argmax(xx>.99*A)
    iend = istart + np.argmax(xx[istart:]<.99*A)
    slp.forcepeg(istart-10, iend)
    yy = slp.complete()
    yy[np.isnan(yy)]=0
    show(tt, xx, yy, 'salpa1')

    print("hello world")
