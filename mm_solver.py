# ******************************************************************************
# Name:    Mystery Matrix Solver
# Author:  Vadim Korotkikh
# Email:   va.korotki@gmail.com
# Date:	   February 2017
#
# Description: Find the R matrix that satisfies certain initial conditions

from numpy.linalg import inv
from numpy.linalg import det
from sympy import *
import numpy as np

import itertools
import time
import sys


def vectors():

	a2	= Matrix([ -sqrt(2)/3, -sqrt(6)/3, -1/3 ])
	a3  = Matrix([ 2*sqrt(2)/3, 0, -1/3 ])
	a4	= Matrix([ -sqrt(2)/3, sqrt(6)/3, -1/3 ])

	sqc3	= sqrt(3)/3

	a2 	= np.asarray(a2)
	a3	= np.asarray(a3)
	a4	= np.asarray(a4)
	# b2	= sqc3 * Matrix([ -1, -1, 1])
	# b3 	= sqc3 * Matrix([ -1, 1, -1])
	# b4	= sqc3 * Matrix([ 1, -1, -1])

	b1	= sqc3	*	np.array([[1], [1], [1]])

	b2 	= sqc3	*	np.array([[-1], [-1], [1]])
	b3 	= sqc3	*	np.array([[-1], [1], [-1]])
	b4 	= sqc3	*	np.array([[1], [-1], [-1]])

	rmat	= Matrix([[ -sqrt(6)/6, sqrt(6)/3, -sqrt(6)/6], [ sqrt(2)/2, 0, -sqrt(2)/2], [sqrt(3)/3, sqrt(3)/3, sqrt(3)/3]])

	nprmat  = np.asmatrix(rmat)
	nprmata = np.asarray(rmat)
	trrmat	= np.transpose(nprmat)

	veccheck( [a2, a3, a4], [b2,b3,b4], nprmat)

	xvec	= np.dot(nprmat, b1)
	ymat	= np.dot(nprmat, trrmat)
	zmat	= np.dot(trrmat, nprmat)

	print("R matrix", det(rmat))
	print(nprmat)
	print("")
	print(xvec)
	print("")
	print(ymat)
	print("")
	print(zmat)

def veccheck( avecs, bvecs, rmat):
	print("")
	print("Checking R matrix")
	for b in bvecs:
		tvec = np.dot(rmat, b)
		print(tvec)
		print("")

vectors()
