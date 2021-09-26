# coding=utf-8
# Numerical methods for computer engineering
# Author: Miguel Blanco Godón, Computer Engineering, 2021

import colours
import lu
import util


def handle():
	valid_inputs = ['0', '1', '2', '3']
	while (True):
		print("IMPLEMENTED METHODS:")
		print("    (1): LU factorization")
		print("    (2): LD(L^t) factorization")
		print("    (3): Cholesky's factorization")
		print("    (0):" + colours.Colours.red(" Quit"))
		x = input("Please select an operation: ")
		while (x not in valid_inputs):
			x = input(colours.Colours.red("Error. Please select a valid operation: "))
		if x == '0':
			print(colours.Colours.green("Bye!"))
			break
		if x == '1':
			file_name = input("Please enter the file path containing the matrix to factorize: ")
			output_file_name = input("Please enter the path of the output file: ")
			try:
				file = open(file_name, 'r')
				a = [[float(n) for n in line.split()] for line in file]
				_a = lu.lu_same(a)
				file.close()
				file = open(output_file_name, 'w+')
				util.Util.write_matrix_to_file(file, _a, "LU factorization:")
				file.close()
				print(colours.Colours.green("Finished!"))
			except Exception as e:
				print(colours.Colours.red(str(e)))
		elif x == '2':
			print("LDLt")
		else:
			print("LLt")

