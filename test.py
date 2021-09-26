#!/usr/bin/python
# coding=utf-8
# Numerical methods for computer engineering
# Author: Miguel Blanco Godón, Computer Engineering, 2021

import colours
import lu
import util

if __name__ == "__main__":
	# testing LU factorization
	for i in range(3):
		print(">>> Testing LU factorization in data set '_test"+str(i+1)+".dat':")
		try:
			file = open("_test"+str(i+1)+".dat", 'r')
			a = [[float(n) for n in line.split()] for line in file]
			file.close()
			l, u = lu.lu(a)
			_lu = lu.lu_same(a)
			if util.Util.equals_lu(l, u, _lu):
				print(colours.Colours.green("Test passed!"))
			else:
				print(colours.Colours.red("Test failed!"))

		except Exception as e:
			print(colours.Colours.red(str(e)))