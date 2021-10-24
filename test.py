#!/usr/bin/python
# coding=utf-8
# Numerical methods for computer engineering
# Author: Miguel Blanco Godón, Computer Engineering, 2021

import colours
import mni
import morse

if __name__ == "__main__":
	# testing LU factorization
	for i in range(1):
#print(">>> Optimal gradient test:")
		try:
			coef = [0.2,0.1,4,1,-1,60,1,1,8,-1,-2,4,700]
			col = [0,0,1,0,1,2,0,1,3,1,2,3,4]
			diag = [0,2,5,8,12]
			a = morse.Morse(coef, col, diag)
			x = [1,2,3,4,5]
			print(">>> Gradient methods test: ")
			b = mni.GradientMethods.const_step(a,x,[1,1,1,1,1], 0.0001, 0.00005, 10000000)
			print("------ CONSTANT STEP: ", b)
			b = mni.GradientMethods.optimal_step(a,x,[1,1,1,1,1], 0.0001, 1000000)
			print("------ OPTIMAL STEP: ", b)
			b = mni.GradientMethods.conjugated(a,x,[1,1,1,1,1], 0.0001, 10)
			print("------ CONJUGATE: ", b)

#a = morse.Morse([4,-1,3,1,-2,3],[0,0,1,0,1,2],[0,2,5])
#a = morse.Morse([4,-5,13,-1,2],[0,0,1,0,2],[0,2,4])
#b = [1,1,1]
			b = [1,1,1,1,1]
			print(">>> Eigenvalues calculation: ")
			x, eing, k = mni.Eigenvalues.inverse_power(a,b,0.0001, 100)
			print(str(x)+", "+str(eing) + " in "  + str(k) + " iterations")
			"""
			f = open("./_test1.dat", 'r')
			coef = [[float(n) for n in line.split()] for line in f][0]
			f.close()
			f = open("./_test2.dat", 'r')
			col = [[int(n) for n in line.split()] for line in f][0]
			f.close()
			f = open("./_test3.dat", 'r')
			diag = [[int(n) for n in line.split()] for line in f][0]
			f.close()
			a = morse.Morse(coef, col, diag)
			a.print_vectorized()
			print("MATRIX FORM:")
			a.print_matrix()
			x = [-1, 2, 3, -4]
			b = morse.Morse.matrix_vector(a, x)
			print("A*x = b; b: ")
			print(b)
	"""

		except Exception as e:
			print(colours.Colours.red(str(e)))
