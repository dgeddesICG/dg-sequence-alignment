## Imports
import numpy as np

from pprint import pprint
import time
import math
from copy import copy

def matchFrames(seq1,seq2,drift):
    # user must provide drift (see period.txt files)
    # Drift provided should be in the order (dx, dy).
    # A positive value moves the object imaged in seq2 up/left to correct for drift
    dx = drift[0]
    dy = drift[1]

    #apply shifts
    rectF = [0,seq1[0].shape[0],0,seq1[0].shape[1]]#X1,X2,Y1,Y2
    rect = [0,seq2[0].shape[0],0,seq2[0].shape[1]]#X1,X2,Y1,Y2

    if dy<=0:
        rectF[0] = -dy
        rect[1] = rect[1]+dy
    else:
        rectF[1] = rectF[1]-dy
        rect[0] = dy
    if dx<=0:
        rectF[2] = -dx
        rect[3] = rect[3]+dx
    else:
        rectF[3] = rectF[3]-dx
        rect[2] = +dx

    seq1 = seq1[:,rectF[0]:rectF[1],rectF[2]:rectF[3]]
    seq2 = seq2[:,rect[0]:rect[1],rect[2]:rect[3]]

    return seq1,seq2
