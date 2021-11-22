# coding=utf-8
# Numerical methods for computer engineering
# Author: Miguel Blanco Godón, Computer Engineering, 2021

import math
import util
import morse
import matplotlib.pyplot as pp
import thomas
import numpy as np

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


class Spline:

	@staticmethod
	def plotCubicSpline(x, w, cn):
		h = []
		n = len(x)-1
		for i in range(n):
			h.append(x[i + 1] - x[i])
		h.append(0)
		diag = [1.0]
		infDiag = [0]
		supDiag = [-2.0]
		wt = [0]
		for i in range(1, n):
			diag.append(2 * (h[i-1] + h[i]))
			infDiag.append(h[i-1])
			supDiag.append(h[i])
			wt.append((6*(w[i+1]-w[i])/h[i]) - (6*(w[i]-w[i-1])/h[i-1]))
		diag.append(1)
		infDiag.append(0) 
		wt.append(cn)
		c = thomas.thomas02(diag, infDiag, supDiag, wt)
		t = x[0]
		i = 0
		s = []
		delta = 1/1000
		T = []
		while t <= x[-1]:
			if t > x[i+1]:
				i+=1
			s.append((c[i]*(((x[i+1]-t)**3)/(6*h[i])))+(c[i+1]*(((t-x[i])**3)/(6*h[i])))+((w[i]/h[i]-c[i]*h[i]/6)*(x[i+1]-t))+((w[i+1]/h[i]-c[i+1]*h[i]/6))*(t-x[i]))
			t += delta
			T.append(t)
		pp.scatter(x,w, color='black')
		pp.plot(T, s, 'r-')
		pp.show()


class Bezier:
	@staticmethod
	def __deCasteljau(t, px, py):
		n = len(px)
		for m in range(1, n):
			for k in range(0, n - m):
				px[k] = (1 - t) * px[k] + t * px[k + 1]
				py[k] = (1 - t) * py[k] + t * py[k + 1]
		return px[0], py[0]

	@staticmethod
	def __deCasteljauRecursive(t, i, n, ck):
		if n == 0:
			return ck[i]
		else:
			return (1-t)*Bezier.deCasteljauRecursive(t, i, n-1, ck) + t*Bezier.deCasteljauRecursive(t, i+1, n-1, ck)
		
	@staticmethod
	def plotBezier(cx, cy):
		n = len(cx)
		bx = []
		by = []
		t = 0
		tf = 1 
		delta = 1/1000
		while (t < tf):
			qx = [float(ck) for ck in cx]
			qy = [float(ck) for ck in cy]
			points = Bezier.__deCasteljau(t, qx, qy)
			bx.append(points[0])
			by.append(points[1])
			t += delta
		# plotting
		for i in range(n):
			pp.scatter(cx[i], cy[i], color="black") # control points
		pp.plot(cx, cy, 'b--') # control polygon
		pp.plot(bx,by, 'r-') # Bezier curve
		pp.show()

	@staticmethod
	def __rationalCasteljau(t, px, py, w):
		n = len(px)
		for m in range(1, n):
			for k in range(0, n - m):
				den = (1-t)*w[k] + t*w[k+1]
				px[k] = ((w[k]*(1-t))/den)*px[k] + ((w[k+1]*t)/den)*px[k+1]
				py[k] = ((w[k]*(1-t))/den)*py[k] + ((w[k+1]*t)/den)*py[k+1]
		return px[0], py[0]

	@staticmethod
	def plotRationalBezier(cx, cy, w):
		n = len(cx)
		bx = []
		by = []
		t = 0
		tf = 1 
		delta = 1/1000
		while (t < tf):
			qx = [float(ck) for ck in cx]
			qy = [float(ck) for ck in cy]
			wk = [float(w_k) for w_k in w]
			points = Bezier.__rationalCasteljau(t, qx, qy, wk)
			bx.append(points[0])
			by.append(points[1])
			t += delta
		# plotting
		for i in range(n):
			pp.scatter(cx[i], cy[i], color="black") # control points
		pp.plot(cx, cy, 'b--') # control polygon
		pp.plot(bx,by, 'r-') # Bezier curve
		pp.show()
	
	@staticmethod
	def __computeBezier(cx, cy):
		n = len(cx)
		bx = []
		by = []
		t = 0
		tf = 1 
		delta = 1/1000
		while (t < tf):
			qx = [float(ck) for ck in cx]
			qy = [float(ck) for ck in cy]
			points = Bezier.__deCasteljau(t, qx, qy)
			bx.append(points[0])
			by.append(points[1])
			t += delta
		return bx, by
	
	@staticmethod
	def interpolate(px, py):
		n = len(px)
		m = []
		delta = 1/(n-1)
		t = [0]
		for i in range(n-1):
			t.append(t[-1]+delta)

		for i in range(n):
			a = []
			for k in range(n):
				a.append(math.comb(n-1,k)*((1-t[i])**(n-k-1))*(t[i]**k))
			m.append(a)

		cx = np.linalg.solve(m, px)
		cy = np.linalg.solve(m, py)
		bx, by = Bezier.__computeBezier(cx, cy)
		# plotting
		pp.scatter(px, py, marker='+', color='blue')
		pp.scatter(cx, cy, color='black')
		pp.plot(cx, cy, '--', color = 'grey')
		pp.plot(bx, by, 'r-')
		pp.show()



class BSpline:

	@staticmethod
	def __createOpenNodeSet(n, p):
		o = []
		k = 0
		while k <= n:
			if k <= p:
				o.append(0)
			else:
				if k <= n:
					o.append(k-p)
				else:
					o.append(n+1-p)
			k += 1 
		return o

	@staticmethod
	def __createClosedNodeSet(start, end, n):
		delta = (end - start) / n
		o = []
		i = start
		while i <= end:
			o.append(i)
			i += delta
		o.append(i)
		return o

	@staticmethod
	def createUniformNodeSet(openSet, start, end, n, p):
		if (openSet):
			return BSpline.__createOpenNodeSet(n, p)
		else:
			return BSpline.__createClosedNodeSet(start, end, n)
	
	
	@staticmethod
	def __deBoor(u, px, py, t, p):
		n = len(px)
		for k in range(0, n):
			for r in range(1, p):
				for j in range(p, r, p-1):
					alpha = (t-u[j+k-p])/(u[j+1+k-r] - u[j+k-p])
					px[j] = (1 - alpha) * px[j-1] + alpha * px[j]
					py[j] = (1 - alpha) * py[j-1] + alpha * py[j]
		return px[0], py[0]


	@staticmethod
	def plotBSplineBases(n, p):
		u = BSpline.createUniformNodeSet(False, 0, 1, n, 0)
		delta = 1/1000
		T = []
		Nk = []
		for k in range(len(u)-1):
			t = u[0]
			N0k = []
			T = []
			while t < u[-1]:
				if (u[k] <= t and t < u[k+1]):
					N0k.append(1)
				else:
					N0k.append(0)
				T.append(t)
				t += delta
			Nk.append(N0k)
		for i in range(p):
			for x in range(len(T)):
				for k in range(u):
					Nk[k][x] = ((T[x]-u[k])/(u[k+p]-u[k]))*Nk[k][x] + ((u[k+p+1]-T[x])/(u[k+p+1]-u[k+1]))*Nk[k+1][x]

		for i in range(len(Nk)):
			pp.plot(T, Nk[i])
		pp.show()



	@staticmethod
	def plotBSpline(cx, cy, u, t):
		n = len(cx)
		bx = []
		by = []
		t = 0
		tf = 1
		delta = 1/1000
		while (t < tf):
			qx = [float(ck) for ck in cx]
			qy = [float(ck) for ck in cy]
			points = BSpline.__deBoor(u, qx, qy, t, 0)

