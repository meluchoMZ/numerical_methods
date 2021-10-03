import math
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
		return x

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
	def gauss_seidel(a, b, x0, max_iter, error_quota):
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
				x[i] = (1/a[i][i])*(b[i]+psum)
				if IterativeMethods.check_error_condition(x,x0,error_quota):
					break;
			k += 1
			x0 = x[:]
		return x
