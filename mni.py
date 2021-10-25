# coding=utf-8
# Numerical methods for computer engineering
# Author: Miguel Blanco Godón, Computer Engineering, 2021

import math
import util
import morse

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



class GradientMethods:

	@staticmethod
	def const_step(a, b, x0, min_error, alpha, max_iterations):
		axo = morse.Morse.matrix_vector(a,x0)
		r = util.Util.substractV(b,axo)
		p = r
		mb = abs(min_error*max(b))
		x = x0
		k = 0
		mr = abs(max(r))
		while (mr >= mb and k < max_iterations):
			ar = morse.Morse.matrix_vector(a, r)
			x = util.Util.sumV(x, [alpha*n for n in r])
			r = util.Util.substractV(r,[alpha*n for n in ar])
			mr2 = abs(max(r))
			beta = (mr2*mr2)/(mr*mr)
			k = k + 1
			mr = mr2
		return x,k
	
	@staticmethod
	def optimal_step(a, b, x0, min_error, max_iterations):
		axo = morse.Morse.matrix_vector(a,x0)
		r = util.Util.substractV(b, axo)
		p = r
		mb = abs(min_error*max(b))
		x = x0
		k = 0
		mr = abs(max(r))
		while (mr >= mb and k < max_iterations):
			ar = morse.Morse.matrix_vector(a, r)
			alpha = util.Util.inner(r, r)/util.Util.inner(r, ar)
			x = util.Util.sumV(x, [alpha*n for n in r])
			r = util.Util.substractV(r,[alpha*n for n in ar])
			mr2 = abs(max(r))
			beta = (mr2*mr2)/(mr*mr)
			k = k + 1
			mr = mr2
		return x,k
	
	@staticmethod
	def conjugated(a, b, x0, min_error, max_iterations):
		axo = morse.Morse.matrix_vector(a,x0)
		r = util.Util.substractV(b, axo)
		p = r
		mb = abs(min_error*max([abs(n) for n in b]))
		x = x0
		k = 0
		mr = max([abs(n) for n in r])
		while (mr >= mb and k < max_iterations):
			ap = morse.Morse.matrix_vector(a, p)
			alpha = util.Util.inner(r, r)/util.Util.inner(p, ap)
			x = util.Util.sumV(x, [alpha*n for n in p])
			r1 = util.Util.substractV(r,[alpha*n for n in ap])
			mr2 = max([abs(n) for n in r])
			beta = util.Util.inner(r1, r1)/util.Util.inner(r, r)
			p = util.Util.sumV(r1, [beta*n for n in p])
			mr = mr2
			r = r1
			k = k + 1
		return x,k

class Eigenvalues:

	@staticmethod
	def inverse_power(m, q, epsilon, max_iterations):
		k = 0
		while (k < max_iterations):
			x, p = GradientMethods.conjugated(m,q,[1 for n in q], epsilon, max_iterations)
			mx = max(x)
			mn = min(x)
			beta = mx if abs(mx) > abs(mn) else mn
			q1 = [n/beta for n in x]
			e = max([abs(n) for n in util.Util.substractV(q1, q)])/max([abs(n) for n in q1])
			if e < epsilon:
				break
			q = q1
			k = k + 1
		return 1/beta, k

