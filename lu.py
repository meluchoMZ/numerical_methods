# LU factorization + systems solver
# coding=utf-8
# Numerical methods for computer engineering
# Author: Miguel Blanco Godón, Computer Engineering, 2021
import util


# Computes LU factorization and returns both matrices in one
def lu_same(a):
	n = len(a)
	# l(k,k) = 1 for all k in [1,n]
	# u(1,k) = a(1, k) for all k in [1, n]
	# l(k,1) = a(k,1)/u(1,1)
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
					#for k in range(r+1, c):
					#	psum += a[r][k] * a[k][c]
					a[r][c] -= psum
	return a


# Computes LU factorization and returns the tuple [l, u]
def lu(a):
	u = []
	l = []
	n = len(a)
	for k in range(0, n):
		u.append([])
		l.append([])
		for kk in range(0, n):
			u[k].append(0)
			l[k].append(0)

	# l(k,k) = 1 for all k in [1,n]
	# u(1,k) = a(1, k) for all k in [1, n]
	# l(k,1) = a(k,1)/u(1,1)
	for k in range(0, n):
		l[k][k] = 1
		u[0][k] = a[0][k]
		l[k][0] = a[k][0] / u[0][0]

	for c in range(1, n):
		for r in range(1, n):
			if c == r:
				u[r][c] = 0
				for k in range(0, c):
					u[r][c] += l[r][k] * u[k][c]
				u[r][c] = a[r][c] - u[r][c]
			else:
				if c < r:
					# compute l
					l[r][c] = 0
					for k in range(0, c):
						l[r][c] += l[r][k] * u[k][c]
					l[r][c] = a[r][c] - l[r][c]
					l[r][c] /= u[c][c]
				else:
					# compute u
					u[r][c] = 0
					for k in range(0, r):
						u[r][c] += l[r][k] * u[k][c]
					for k in range(r+1, c):
						u[r][c] += l[r][k] * u[k][c]
					u[r][c] = a[r][c] - u[r][c]
					u[r][c] /= l[r][r]
	return l, u

def solve_l_system(_lu, b):
	y = [b[0]]
	for k in range(1, len(b)):
		psum = 0
		for kk in range(k):
			psum += _lu[k][kk] * y[kk]
		y.append(b[k] - psum)
	return y

def solve_u_system(_lu, y):
	l = len(y)-1
	x = [y[l]/_lu[l][l]]
	for k in range(l-1, -1, -1):
		psum = 0
		for kk in range(l, k, -1):
			psum += _lu[k][kk] * x[l-kk]
		x.append((y[k]-psum)/_lu[k][k])
	x.reverse()
	return x

def solve_lu(_lu, b):
	y = solve_l_system(_lu,b)
	x = solve_u_system(_lu,y)
	return x
