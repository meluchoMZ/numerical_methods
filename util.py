# coding=utf-8
# Numerical methods for computer engineering
# Author: Miguel Blanco Godón, Computer Engineering, 2021

class Util:

	@staticmethod
	def print_matrix(m):
		for k in range(0, len(m)):
			for kk in range(0, len(m)):
				print("%5s" % (str(m[k][kk])) , end=" ")
			print("")

	@staticmethod
	def multiply_matrices(m,n):
		o = []
		for k in range(0, len(m)):
			o.append([])
			for kk in range(0, len(m)):
				o[k].append(0)
		for k in range(0,len(m)):
			for kk in range(0, len(n)):
				for kkk in range(0, len(n)):
					o[k][kk] += m[k][kkk]*n[kkk][kk]
		return o

	@staticmethod
	def matrix_vector(a,b):
		o = []
		for k in range(len(b)):
			o.append(0)
			for kk in range(len(a)):
				o[k] += a[k][kk]*b[kk]
		return o

	@staticmethod
	def sumV(a,b):
		o = []
		for k in range(len(a)):
			o.append(a[k]+b[k])
		return o

	@staticmethod
	def substractV(a,b):
		o = []
		for k in range(len(a)):
			o.append(a[k]-b[k])
		return o

	@staticmethod
	def inner(a,b):
		inn = 0;
		for k in range(len(a)):
			inn += a[k]*b[k]
		return inn

	@staticmethod
	def equals(m,n):
		if len(m) != len(n):
			return False
		for k in range(len(m)):
			for kk in range(len(n)):
				if len(m) != len(n):
					return False
		return True

	@staticmethod
	def equals_lu(l,u,lu_same):
		if len(l) != len(u) != len(lu_same):
			return False
		for k in range(len(l)):
			for kk in range(len(u)):
				if kk < k and l[k][kk] != lu_same[k][kk]:
					return False
				if kk >= k and u[k][kk] != lu_same[k][kk]:
					return False
		return True

	@staticmethod
	def write_matrix_to_file(file, m, message):
		file.write(message+"\n")
		for i in range(len(m)):
			for j in range(len(m[i])):
				file.write("a["+str(i)+"]["+str(j)+"]="+str(m[i][j]) + " ")
			file.write("\n")

	@staticmethod
	def write_vector_to_file(file, m, message):
		file.write(message+"\n")
		for i in range(len(m)):
			file.write("x["+str(i)+"]="+str(m[i]) + " ")
		file.write("\n")
