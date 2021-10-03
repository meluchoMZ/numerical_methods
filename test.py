#!/usr/bin/python
# coding=utf-8
# Numerical methods for computer engineering
# Author: Miguel Blanco Godón, Computer Engineering, 2021

import colours
import iterative_methods
import lu
import util

if __name__ == "__main__":
	# testing LU factorization
	for i in range(1):
		print(">>> Testing LU factorization in data set '_test" + str(i + 1) + ".dat':")
		try:
			file = open("_test" + str(i + 1) + ".dat", 'r')
			a = [[float(n) for n in line.split()] for line in file]
			file.close()
			"""
			l, u = lu.lu(a)
			print("L:")
			util.Util.print_matrix(l)
			print("\nU:")
			util.Util.print_matrix(u)
			print("\n")
			print("\nLU:")
			"""
			# b = [2, 3, 4, 1]
			# b = [8, 7, -8, 9]
			# b = [5, -1, 7, -5, 12]
			"""
			_lu = lu.lu_same(a)
			
			util.Util.print_matrix(_lu)
			print("\n")
			x = lu.solve_lu(_lu,b)
			print("X:")
			print(x)
			if util.Util.equals_lu(l, u, _lu):
				print(colours.Colours.green("Test passed!"))
			else:
				print(colours.Colours.red("Test failed!"))
				"""
			"""
			#c = [[4, -1, 2], [1, 4, 0], [-1, 1, 3]]
			c = [[10,-1,2,0],[-1,11,-1,3],[2,-1,10,-1],[0,3,-1,8]]
			print("A:")
			util.Util.print_matrix(c)
			print("B:")
			#b = [1,3,2]
			b = [6,25,-11,15]
			print(b)
			x = iterative_methods.IterativeMethods.jacobi(c, b, [0, 0, 0, 0], 10, 0.001)
			print("X:")
			print(x)
			xgauss = iterative_methods.IterativeMethods.gauss_seidel(c, b, [0,0,0,0], 10, 0.001)
			print("X GAUSS-SEIDEL")
			print(xgauss)

			"""
			"""			import mni
			c = [[4,3,0],[3,4,-1],[0,-1,4]]
			print("A:")
			util.Util.print_matrix(c)
			print("B:")
			#b = [1,3,2]
			b = [24,30,-24]
			print(b)
			x,k = mni.IterativeMethods.jacobi(c, b, [1,1,1], 100, 0.001)
			print("X:"+str(k)+" iterations")
			print(x)
			xgauss, k = mni.IterativeMethods.relaxation(c, b, [1,1,1], 100, 0.001, 1)
			print("X GAUSS-SEIDEL" + str(k) + " iterations")
			print(xgauss)
			xrel, k = mni.IterativeMethods.relaxation(c,b,[1,1,1],100,0.001,1.25)
			print("X REL:" + str(k) + " iterations")
			print(xrel)
			print("check error condition")
			print(mni.IterativeMethods.check_error_condition([1.0001,1.9998,-0.9998,0.9998],[0.9997,2.0004,-1.0004,1.0006],0.001))
			"""
			#c = [[4,-1,1],[-1,4.25,2.75],[1,2.75,4]]
			#c = [[1,1,0,3],[2,1,-1,1],[3,-1,-1,2],[-1,2,3,-1]]
			#import mni
			#c = mni.Factorizations.cholesky(c)
			#c = lu.lu_same(c)
			#util.Util.print_matrix(c)
			#b = [8,7,14,-7]
			#b = [1,2,3]
			#x = mni.Factorizations.solve_cholesky(c, b)
			#print(x)
			import mni
			a = [[4,3,0],[3,4,-1],[0,-1,4]]
			b = [24,30,-24]
			x,k = mni.IterativeMethods.relaxation(a,b,[1,1,1],2000,0.0000001, 1)
			print(x)
			print(k, "iters")
			x,k = mni.IterativeMethods.relaxation(a,b,[1,1,1],2000,0.0000001, 1.25)
			print(x)
			print(k, "iters")

		except Exception as e:
			print(colours.Colours.red(str(e)))
