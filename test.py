#!/usr/bin/python
# coding=utf-8
# Numerical methods for computer engineering
# Author: Miguel Blanco Godón, Computer Engineering, 2021

import colours
from mni import Spline
from mni import Bezier
from mni import BSpline
import morse
import matplotlib.pyplot as pp

if __name__ == "__main__":
	# testing LU factorization
	for i in range(1):
#print(">>> Optimal gradient test:")
		try:
#		a = BSpline.createUniformNodeSet(False, 1, 2, 10)
	 #		print(a)
			cx = [0.0, 1.0, 1.0, 1.0, 2.0]
			cy = [0.0, 0.0, 1.0, 2.0, 2.0]
			BSpline.plotBSpline(cx, cy)
			cx = [0.0, 0.0, 1.5, 3.0, 3.0]
			cy = [0.0, 2.0, 3.0, 2.0, 0.0]
			BSpline.plotBSpline(cx,cy)
			cx = [0.0, 1.0, 1.0, 1.0, 2.0]
			cy = [0.0, 0.0, 1.0, 2.0, 2.0]
			BSpline.plotBSpline(cx,cy)
			BSpline.plotBSplineBases(5,0)
			BSpline.plotBSplineBases(8,1)
			BSpline.plotBSplineBases(5,2)
			BSpline.plotBSplineBases(5,3)
			BSpline.plotBSplineBases(5,4)
			BSpline.plotBSplineBases(14, 3)
			x = [1, 0, 1, 2, 2]
			y = [0, 1, 2, 2, 0]
			Bezier.interpolate(x,y)
			x = [-4,-3,-2,-1,0,1,2,3,4]
			w = [16,9,4,1,0,1,4,9,16]
			Bezier.interpolate(x,w)
			Spline.plotCubicSpline(x, w, 4) 
			x = [-1, 0, 1, 3, 4]
			w = [6, 3, 6, 38, 77]
			Spline.plotCubicSpline(x, w, 4) 
			Bezier.interpolate(x,w)
			cx = [0.0, 0.0, 1.0]
			cy = [1.0, 0.0, 0.0]
			Bezier.plotBezier(cx, cy)
			cx = [0.0, 0.0, 1.5, 3.0, 3.0]
			cy = [0.0, 2.0, 3.0, 2.0, 0.0]
			Bezier.plotBezier(cx, cy)
			cx = [0.0, -2.0, -3.0, -2.0, 0.0]
			cy = [0.0, 0.0, 1.5, 3.0, 3.0]
			Bezier.plotBezier(cx, cy)
			cx = [0.0, 1.0, 1.0, 1.0, 2.0]
			cy = [0.0, 0.0, 1.0, 2.0, 2.0]
			Bezier.plotBezier(cx, cy)
			Bezier.plotRationalBezier(cx, cy, [1, 1.8, 1, 1.8, 1])
			cx = [0.0, 1.0, 2.0, 3.0] 
			cy = [0.0, 2.0, 2.0, 0.0] 
			Bezier.plotBezier(cx, cy)
			cx = [ 3.0, 5.0, 6.0]
			cy = [ 0.0, 0.0, 2.0]
			Bezier.plotBezier(cx, cy)

		except Exception as e:
			print(colours.Colours.red(str(e)))
