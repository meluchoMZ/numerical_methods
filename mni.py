# Numerical methods for computer engineering
# Author: Miguel Blanco Godón, Computer Engineering, 2021

import math

class Factorizations:
	# Factorizations

	@staticmethod
	def lu(a):
		n = len(a)
		for k in range(1, n):
			a[k][0] = a[k][0] / a[0][0]
		for c in range(1, n):
			for r in range(1, n):
				psum = 0
				if c == r:
					for k in range(0, c):
						psum += a[r][k] * a[k][c]
					a[r][c] -= psum
				else:
					if c < r:
						# compute l
						for k in range(0, c):
							psum += a[r][k] * a[k][c]
						a[r][c] -= psum
						a[r][c] /= a[c][c]
					else:
						# compute u
						for k in range(0, r):
							psum += a[r][k] * a[k][c]
						a[r][c] -= psum
		return a

	@staticmethod
	def solve_lu(_lu, b):
		y = [b[0]]
		#decline
		for k in range(1, len(b)):
			psum = 0
			for kk in range(k):
				psum += _lu[k][kk] * y[kk]
			y.append(b[k] - psum)
		l = len(y) - 1
		#trace back
		x = [y[l] / _lu[l][l]]
		for k in range(l - 1, -1, -1):
			psum = 0
			for kk in range(l, k, -1):
				psum += _lu[k][kk] * x[l - kk]
			x.append((y[k] - psum) / _lu[k][k])
		x.reverse()
		return x

	@staticmethod
	def cholesky(a):
		for j in range(len(a)):
			for i in range(len(a)):
				if i==j:
					# diagonal calculation
					qsum = 0
					for k in range(i):
						qsum += a[i][k] * a[i][k]
					a[i][i] = math.sqrt(a[i][i] - qsum)
				else:
					if j < i:
						# l calculation
						psum = 0
						for k in range(j):
							psum += a[i][k] * a[j][k]
						a[i][j] = (a[i][j]-psum)/a[j][j]
						a[j][i] = a[i][j]
		return a

	@staticmethod
	def solve_cholesky(a, b):
		y = []
		y.append(b[0]/a[0][0])
		#decline
		for k in range(1, len(b)):
			psum = 0
			for kk in range(k):
				psum += a[k][kk] * y[kk]
			y.append((b[k] - psum)/a[k][k])
		print("Y:")
		print(y)
		#trace back
		l = len(y)-1
		x = [y[l] / a[l][l]]
		for k in range(l - 1, -1, -1):
			psum = 0
			for kk in range(l, k, -1):
				psum += a[k][kk] * x[l - kk]
			x.append((y[k] - psum) / a[k][k])
		x.reverse()
		return x



class IterativeMethods:
	@staticmethod
	def jacobi(a, b, x0, max_iter, error_quota):
		n = len(b)
		k = 0
		x = x0[:]
		while k < max_iter:
			for i in range(n):
				psum = 0
				for j in range(i):
					psum += (-a[i][j]*x0[j])
				for j in range(i+1, n):
					psum += (-a[i][j]*x0[j])
				x[i] = (1/a[i][i])*(psum+b[i])
			if IterativeMethods.check_error_condition(x, x0, error_quota):
				break
			k += 1
			x0 = x[:]
		return x,k

	@staticmethod
	def check_error_condition(x, y, epsilon):
		num = []
		den = []
		for i in range(len(x)):
			num.append(x[i]-y[i])
			num[i] = num[i]*num[i]
			den.append(x[i]*x[i])
		num = math.sqrt(sum(num))
		den = math.sqrt(sum(den))
		return (num/den) < epsilon

	@staticmethod
	def relaxation(a, b, x0, max_iter, error_quota, w):
		n = len(b)
		k = 0
		x = x0[:]
		while k < max_iter:
			for i in range(n):
				psum = 0
				for j in range(i):
					psum += (-a[i][j]*x[j])
				for j in range(i+1, n):
					psum += (-a[i][j]*x0[j])
				x[i] = (1-w)*x0[i] + (w/a[i][i])*(b[i]+psum)
			if IterativeMethods.check_error_condition(x,x0,error_quota):
				break
			k += 1
			x0 = x[:]
		return x, k
