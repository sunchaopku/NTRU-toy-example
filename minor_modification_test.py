#q: the modulus
#n: degree of the number field
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler


def poly_to_anti(n, f):
	g = []
	g.append(-f[ n - 1])
	for i in range(n - 1):
		g.append(f[i])
	return g
def poly_to_matrix(n, f):
	M = matrix(ZZ, n, n)
	h = []
	for i in range(n):
		h.append(f[i])
	for i in range(n):
		M[i] = copy(h)
		h = poly_to_anti(n, h)
	return M

def matrix_to_poly(n, M):
	R = PolynomialRing(ZZ, 'x')
	x = R.gen()
	P = R.quotient(x^n + 1, 'x')
	temp = []
	for i in range(n):
		temp.append(M[0,i])
	res =P(temp)
	return res
def zpoly_to_qpoly(n, f):
	R1 = PolynomialRing(QQ, 'x')
	z = R1.gen()
	P1 = R1.quotient(z^n + 1, 'x')
	temp = []
	for i in range(n):
		temp.append(f[i])
	return P1(temp)

def zpoly_to_qpoly_no_mod(n, f):
	R1 = PolynomialRing(QQ, 'x')
	temp = []
	for i in range(n):
		temp.append(f[i])
	return R1(temp)

def P_to_R(n, f):
	R1 = PolynomialRing(QQ, 'x')
	temp = []
	for i in range(n):
		temp.append(f[i])
	return R1(temp)

def lcmf(n, Rf_temp, RN_temp):
	res = 1
	for i in range(n):
		res = lcm(res, Rf_temp[i].denominator())
		res = lcm(res, RN_temp[i].denominator())
	return res

def round_f(n, f):
	R1 = PolynomialRing(QQ, 'x')
	temp = []
	for i in range(n):
		t = floor(f[i])
		if f[i] - t < 0.5:
			temp.append(t)
		else:
			temp.append(t + 1)
	return R1(temp)


from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
def DLPkeygen(alpha):
	n = 512
	q = 12289
	cnt = 0
	#def DLPkeygen(n, q):
	xx = uniform(0,1)
	xx = xx^(1/n)
	print("xx is ", xx)
	sigma = alpha * sqrt(q/(2 * n))
	D = DiscreteGaussianDistributionIntegerSampler(sigma, c = 0)
	radius = alpha * sqrt(q)

	R = PolynomialRing(QQ, 'x')
	x = R.gen()
	P = R.quotient(x^n + 1, 'y')
	while 1:
			cnt = cnt + 1
			#print("current repetition number is ", cnt) 
			direction_f = []
			direction_g = []
			T = RealDistribution('gaussian', 1)
			for i in range(n):
				direction_f.append( T.get_random_element())
				direction_g.append( T.get_random_element())
			len_fg = 0 
	
			for i in range(n):
				len_fg = len_fg + direction_f[i]^2
				len_fg = len_fg + direction_g[i]^2
			len_fg = sqrt(len_fg)
			for i in range(n):
				direction_f[i] = direction_f[i] / len_fg
				direction_g[i] = direction_g[i] / len_fg
		#	print(direction_f)
			f_vec = []
			g_vec = []
			for i in range(n):
				f_vec.append(0)
				g_vec.append(0)
			for i in range(n):
				f_vec[i] = round(RR(xx * radius * direction_f[i]))
				g_vec[i] = round(RR( xx * radius * direction_g[i]))
		#	print(f_vec)
			f = zpoly_to_qpoly(n, P(f_vec))
			g = zpoly_to_qpoly(n, P(g_vec))
			Mf = poly_to_matrix(n, f)
			Mg = poly_to_matrix(n, g)
			len1 = 0
			for i in range(n):
				len1 = len1 + f[i] * f[i] + g[i] * g[i]
			len1 = sqrt(len1)
			f_ = zpoly_to_qpoly(n, matrix_to_poly(n, Mf.T))
			g_ = zpoly_to_qpoly(n, matrix_to_poly(n, Mg.T))
			T = f * f_ + g * g_ 
		#	print(T)
			t1 = q * f_  / T 
			t2 = q * g_ / T
			print("f G - gF", f * t1 + g * t2)
			len2 = 0
			for i in range(n):
				len2 = len2 + t1[i]^2 + t2[i]^2
			len2 = sqrt(len2)
			print(RR(len1 * len2))
			if len1 > alpha * sqrt(q) or len2 > alpha * sqrt(q):
				continue
			break
	print("alpha is ", alpha)
	print("total repetition number is ", cnt) 

alpha = 1.17
for i in range(10):
	for j in range(10):
		DLPkeygen(alpha)

	alpha = alpha - 0.01





