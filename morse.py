class Morse:
	def __init__(self, coeff, col, diag):
		self.__coeff = coeff
		self.__col = col
		self.__diag = diag
	def getCoeff(self):
		return self.__coeff 
	def setCoeff(self, newCoeff):
		self.__coeff = newCoeff
	def getCol(self):
		return self.__col
	def setCol(self, newCol):
		self.__col = newCol
	def getDiag(self):
		return self.__diag
	def setDiag(self, newDiag):
		self.__diag = newDiag

	def print_vectorized(self):
		print("COEFFICIENTS: ")
		print(self.__coeff)
		print("COLUMNS: ")
		print(self.__col)
		print("DIAGONALS: ")
		print(self.__diag)
	
	def print_matrix(self):
		c = 0;
		for k in range(len(self.__diag)):
			# 0
			# 0, 1 
			# 0, 1, 2...
			for kk in range(k+1):
				if kk == self.__col[c]:
					print(self.__coeff[c], end=" ")
					c = c+1
				else:
					print(0, end=" ")
			print("")

	@staticmethod
	def matrix_vector(a, x):
		coef = a.getCoeff()
		col = a.getCol()
		diag = a.getDiag()
		b = [0 for n in x]

		for k in range(len(x)):
			k1 = diag[k-1]+1
			k2 = diag[k]
			for kk in range(k1,k2):
				b[k] += coef[kk]*x[col[kk]]
				b[col[kk]] += coef[kk]*x[k]
			b[k] += coef[k2]*x[k]

		return b
