# coding=utf-8
# Numerical methods for computer engineering
# Author: Miguel Blanco Godón, Computer Engineering, 2021

import colours
import mni
import util
import time
import morse


def handle():
	valid_inputs = ['0', '1', '2', '3', '4', '5', '6', '7', '8']
	while (True):
		print("IMPLEMENTED METHODS:")
		print("    (1): LU factorization")
		print("    (2): Cholesky's factorization")
		print("    (3): Jacobi's iterative method")
		print("    (4): Relaxation iterative method")
		print("    (5): Constant step gradient method")
		print("    (6): Optimal step gradient method")
		print("    (7): Conjugated gradient method")
		print("    (8): Inverse power method")
		print("    (0):" + colours.Colours.red(" Quit"))
		x = input("Please select an operation: ")
		while (x not in valid_inputs):
			x = input(colours.Colours.red("Error. Please select a valid operation: "))
		if x == '0':
			print(colours.Colours.green("Bye!"))
			break
		if x == '1':
			file_name = input("Please enter the file path containing the matrix to factorize: ")
			file_name2 = input("Please enter the file path containing the 'b' vector (Ax=b): ")
			output_file_name = input("Please enter the path of the output file: ")
			try:
				file = open(file_name, 'r')
				a = [[float(n) for n in line.split()] for line in file]
				file.close()
				t = time.time()
				_a = mni.Factorizations.lu(a)
				t = time.time() - t
				file = open(file_name2, 'r')
				b = [[float(n) for n in line.split()] for line in file][0]
				file.close()
				t2 = time.time()
				x = mni.Factorizations.solve_lu(_a,b)
				t2 = time.time() - t2
				file = open(output_file_name, 'w+')
				util.Util.write_matrix_to_file(file, _a, "LU factorization:")
				util.Util.write_vector_to_file(file, x, "Solution to system:")
				file.close()
				print(colours.Colours.green("Finished factorization in " +str(t) +" seconds"))
				print(colours.Colours.green("System solved in "+str(t2) +" seconds"))
				print(colours.Colours.green("Total elapsed time: "+str(t+t2)+ " seconds"))
			except Exception as e:
				print(colours.Colours.red(str(e)))
		elif x == '2':
			file_name = input("Please enter the file path containing the matrix to factorize: ")
			file_name2 = input("Please enter the file path containing the 'b' vector (Ax=b): ")
			output_file_name = input("Please enter the path of the output file: ")
			try:
				file = open(file_name, 'r')
				a = [[float(n) for n in line.split()] for line in file]
				file.close()
				t = time.time()
				_a = mni.Factorizations.cholesky(a)
				t = time.time() - t
				file = open(file_name2, 'r')
				b = [[float(n) for n in line.split()] for line in file][0]
				file.close()
				t2 = time.time()
				x = mni.Factorizations.solve_cholesky(_a,b)
				t2 = time.time() - t2
				file = open(output_file_name, 'w+')
				util.Util.write_matrix_to_file(file, _a, "Cholesky factorization:")
				util.Util.write_vector_to_file(file, x, "Solution to system:")
				file.close()
				print(colours.Colours.green("Finished factorization in " +str(t) +" seconds"))
				print(colours.Colours.green("System solved in "+str(t2) +" seconds"))
				print(colours.Colours.green("Total elapsed time: "+str(t+t2)+ " seconds"))
			except Exception as e:
				print(colours.Colours.red(str(e)))
		elif x == '3':
			file_name = input("Please enter the file path containing the system matrix: ")
			file_name2 = input("Please enter the file path containing the 'b' vector (Ax=b): ")
			file_name3 = input("Please enter the file paht containing the vector of initial values: ")
			maxiter = int(input("Please enter the maximum number of iterations: "))
			error = float(input("Please enter the error constraint: "))
			output_file_name = input("Please enter the path of the output file: ")
			try:
				file = open(file_name, 'r')
				a = [[float(n) for n in line.split()] for line in file]
				file.close()
				file = open(file_name2, 'r')
				b = [[float(n) for n in line.split()] for line in file][0]
				file.close()
				file = open(file_name3, 'r')
				x0 = [[float(n) for n in line.split()] for line in file][0]
				file.close()
				t2 = time.time()
				x,iters = mni.IterativeMethods.jacobi(a,b,x0,maxiter,error)
				t2 = time.time() - t2
				file = open(output_file_name, 'w+')
				util.Util.write_vector_to_file(file, x, "Solution to system:")
				file.close()
				print(colours.Colours.green("System solved in "+str(iters) +" iterations"))
				print(colours.Colours.green("Total elapsed time: "+str(t2)+ " seconds"))
			except Exception as e:
				print(colours.Colours.red(str(e)))
		elif x == '4':
			file_name = input("Please enter the file path containing the system matrix: ")
			file_name2 = input("Please enter the file path containing the 'b' vector (Ax=b): ")
			file_name3 = input("Please enter the file paht containing the vector of initial values: ")
			maxiter = int(input("Please enter the maximum number of iterations: "))
			error = float(input("Please enter the error constraint: "))
			w = float(input("Please enter the 'Omega' step parameter: " ))
			output_file_name = input("Please enter the path of the output file: ")
			try:
				file = open(file_name, 'r')
				a = [[float(n) for n in line.split()] for line in file]
				file.close()
				file = open(file_name2, 'r')
				b = [[float(n) for n in line.split()] for line in file][0]
				file.close()
				file = open(file_name3, 'r')
				x0 = [[float(n) for n in line.split()] for line in file][0]
				file.close()
				t2 = time.time()
				x,iters = mni.IterativeMethods.relaxation(a,b,x0,maxiter,error,w)
				t2 = time.time() - t2
				file = open(output_file_name, 'w+')
				util.Util.write_vector_to_file(file, x, "Solution to system:")
				file.close()
				print(colours.Colours.green("System solved in "+str(iters) +" iterations"))
				print(colours.Colours.green("Total elapsed time: "+str(t2)+ " seconds"))
			except Exception as e:
				print(colours.Colours.red(str(e)))
			except ValueError as v:
				print(colours.Colours.red(str(e)))
		elif x == '5' or x == '6' or x == '7':
			try:
				file_name1 = input("Please enter the file path containing coefficients vector: ")
				file_name2 = input("Please enter the file path containing columns indices vector: ")
				file_name3 = input("Please enter the file path containing diagonal indices vector: ")
				file_name4 = input("Please enter the file path containing second member of the system: ")
				file_name5 = input("Please enter the file path containing initial vector: ")
				it = int(input("Please enter max number of iterations: "))
				error = float(input("Please enter max error: "))
				if x == '5': 
					step = float(input("Please enter the constant step: "))
				file_name6 = input("Please enter the file path to write the output vector: ")
			
				f = open(file_name1,'r')
				coef = [float(n) for n in f]
				f.close()
				f = open(file_name2,'r')
				col = [int(n) for n in f]
				f.close()
				f = open(file_name3,'r')
				diag = [int(n) for n in f]
				f.close()
				f = open(file_name4,'r')
				b = [float(n) for n in f]
				f.close()
				f = open(file_name5,'r')
				x0 = [float(n) for n in f]
				f.close()
				a = morse.Morse(coef, col, diag)
				if x == '5':
					t = time.time()
					x, k = mni.GradientMethods.const_step(a,b,x0,error,step,it)
					t = time.time() - t
				if x == '6':
					t = time.time()
					x, k = mni.GradientMethods.optimal_step(a,b,x0,error,it)
					t = time.time() - t
				if x == '7':
					t = time.time()
					x, k = mni.GradientMethods.conjugated(a,b,x0,error,it)
					t = time.time() - t

				f = open(file_name6, 'w+')

				util.Util.write_vector_to_file(f, x, "Solution to system:")
				f.close()
				print("System solved in ", k ," iterations. Output saved on file: ", file_name6)
				print("Elapsed time: ", t, " seconds")

			except Exception as e:
				print(colours.Colours.red(str(e)))

		elif x == '8':
			try:
				file_name1 = input("Please enter the file path containing coefficients vector: ")
				file_name2 = input("Please enter the file path containing columns indices vector: ")
				file_name3 = input("Please enter the file path containing diagonal indices vector: ")
				file_name4 = input("Please enter the file path containing initial vector: ")
				it = int(input("Please enter max number of iterations: "))
				error = float(input("Please enter max error: "))
				file_name5 = input("Please enter the file path to write the output vector: ")
				f = open(file_name1, 'r')
				coef = [float(n) for n in f]
				f.close()
				f = open(file_name2, 'r')
				col = [int(n) for n in f]
				f.close()
				f = open(file_name3, 'r')
				diag = [int(n) for n in f]
				f.close()
				f = open(file_name4, 'r')
				init = [float(n) for n in f]
				f.close()
				a = morse.Morse(coef, col, diag)
				t = time.time()
				x, k = mni.Eigenvalues.inverse_power(a, init, error, it) 
				t = time.time() - t
				f = open(file_name5, 'w+')
				util.Util.write_vector_to_file(f, x, "Solution to system:")
				f.close()
				print("System solved in ", k ," iterations. Output saved on file: ", file_name5)
				print("Elapsed time: ", t, " seconds")
			except Exception as e:
				print(colours.Colours.red(str(e)))
		else:
		  	print("404 not found")

