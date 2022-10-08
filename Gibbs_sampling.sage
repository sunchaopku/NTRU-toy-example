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

def coef_test(alpha, n_, f_, g_ ):
	n = n_
	q = 12289
	cnt = 0

	R = PolynomialRing(QQ, 'x')
	x = R.gen()
	P = R.quotient(x^n + 1, 'y')
	while 1:
			cnt = cnt + 1
			#print("current repetition number is ", cnt) 
			f_vec = []
			g_vec = []
			for i in range(n):
				f_vec.append(round(f_[i, 0]))
				g_vec.append(round(g_[i, 0]))
			f = zpoly_to_qpoly(n, P(f_vec))
			g = zpoly_to_qpoly(n, P(g_vec))
		#	Mf = poly_to_matrix(n, f)
		#	Mg = poly_to_matrix(n, g)
			len1 = 0
			for i in range(n):
				len1 = len1 + f[i] * f[i] + g[i] * g[i]
			len1 = sqrt(len1)
			f_ad_vec = []
			g_ad_vec = []
			for i in range(n):
				if i == 0:
					f_ad_vec.append(f_vec[i])
					g_ad_vec.append(g_vec[i])
				else:
					f_ad_vec.append(- f_vec[n - i] )
					g_ad_vec.append(- g_vec[n - i] )

			f_ = zpoly_to_qpoly(n, P(f_ad_vec))
			g_ = zpoly_to_qpoly(n, P(g_ad_vec))
			T = f * f_ + g * g_ 
			t1 = q * f_  / T 
			t2 = q * g_ / T
	#		print("f G - gF", f * t1 + g * t2)
			len2 = 0
			for i in range(n):
				len2 = len2 + t1[i]^2 + t2[i]^2
			len2 = sqrt(len2)
		#	print(RR(len1 * len2))
		#	if len1 > alpha * sqrt(q):
		#		print("unnecessary sampling")
		#	else:
		#		print("necessary sampling")
			if len1 > alpha * sqrt(q) or len2 > alpha * sqrt(q):
				return False
			else:
				return True


def DLPkeygen(alpha,n_):
	n = n_
	q = 12289
	cnt = 0
	#def DLPkeygen(n, q):
#	sigma = alpha * sqrt(q/(2 * n))
#	D = DiscreteGaussianDistributionIntegerSampler(sigma, c = 0)


#	R = PolynomialRing(QQ, 'x')
#	x = R.gen()
#	P = R.quotient(x^n + 1, 'y')
	f_vec = []
	for i in range(2 * n):
		f_vec.append(alpha^2 * q / 2)
	while 1:
			cnt = cnt + 1
			#print("current repetition number is ", cnt) 
			
			#g_vec = []
			#for i in range(2 * n):
			#	f_vec.append(alpha^2 * q / 2)
			#	g_vec.append(0)
		#	len1 = 0
		#	for i in range(n):
		#		len1 = len1 + f[i] * f[i] + g[i] * g[i]
		#	len1 = sqrt(len1)
			for i in range( 2 * n):
				#print("current i is ", i)
				temp_sum_1 = 0
				for j in range(2 * n):
					if j != i:
						temp_sum_1 = temp_sum_1 + f_vec[j]
				upper_bound = n * alpha * alpha * q - temp_sum_1
				temp_sum_2 = 0
				for j in range(n):
					temp_sum_2 = temp_sum_2 + q^2 / (  f_vec[2 * j] + f_vec[2 * j + 1] )
				lower_bound = 0
				if i % 2  == 0:
					temp_sum_2 = temp_sum_2 - q^2 /(f_vec[i] + f_vec[i + 1] )
					lower_bound = q^2 / ( n * alpha * alpha * q - temp_sum_2 ) - f_vec[i + 1]
				else:
					temp_sum_2 = temp_sum_2 - q^2 /(f_vec[i] + f_vec[ i - 1] )
					lower_bound = q^2 / ( n * alpha * alpha * q - temp_sum_2 ) - f_vec[i - 1]
				lower_bound = max(lower_bound, 0)
				#if lower_bound < 0:
				#	print(lower_bound)
				#	print("lower_bound below 0")
				#T = RealDistribution('uniform', [lower_bound, upper_bound])
			#	f_vec[i] = ZZ.random_element( ceil(lower_bound), floor(upper_bound))
				temp_val = 0
				if upper_bound <= lower_bound: 
					#f_vec[i] = 0
					continue
					print("upper_bound", upper_bound)
					print("lower_bound", lower_bound)
					print("WWWWWWWW")
				while 1:
					temp_val = Rq(ZZ.random_element( ceil(lower_bound), floor(upper_bound)))
					if temp_val.is_square() == True:
						break
				f_vec[i] =  ZZ(temp_val)


	#		print("f G - gF", f * t1 + g * t2)
			len2 = 0
			for i in range(2 * n):
				len2 = len2 + f_vec[i]
			len1 = 0
			for i in range(n):
				len1 = len1 + q^2/ (f_vec[2 * i] + f_vec[2 * i + 1] )
		#	print(RR(len1 * len2))
			if len2 > n * alpha * alpha * q or len1 > n * alpha * alpha * q :
				print("unnecessary sampling")
			else:
				print("necessary sampling")
			#print(f_vec)
			print("current repetition number is ", cnt)
			if len2 > n * alpha * alpha * q or len1 > n * alpha * alpha * q :
				continue
			break
			print("current repetition number is ", cnt)
	print("alpha is ", alpha)
	print("total repetition number is ", cnt) 
	return f_vec
q = 12289
n = 4
Rq = GF(q)
alpha = 1.1
output_fg = DLPkeygen(alpha, n )
delta = Rq.multiplicative_generator()  
f = matrix(Rq, n, 1)
g = matrix(Rq, n, 1)
for i in range(2 * n):
	temp_val = Rq(output_fg[i])
	print(temp_val)
	temp_val = temp_val.square_root()
	print(temp_val)
	if i % 2 == 0:
		f[i/2, 0] = temp_val
	else:
		g[(i-1)/2, 0] = temp_val
w = delta^( (q - 1) / (2 * n) )
MF = matrix(Rq, n, 1)
MG = matrix(Rq, n, 1)
Trans = matrix(Rq, n, n )
for i in range(n):
	temp_w = w^(2 * i + 1)
	for j in range(n):
		Trans[i, j] = temp_w^j
Inv_Trans = Trans.inverse()
f_coef = Inv_Trans * f 
g_coef = Inv_Trans * g
print(coef_test(alpha, n, f_coef, g_coef))
#for i in range(100):
#	if is_prime(i * 1024 +1):
#		print(i*1024 + 1)
#		break

#for i in range(10):
#	for j in range(10):
#		DLPkeygen(alpha)

#	alpha = alpha - 0.01


#from sage.rings.polynomial.convolution import _negaconvolution_fft

#n = 3 # degree 1024
#Rq = GF(12289)

#a = FastFourierTransform(2 * n)

#R.<X> = PolynomialRing(Rq)
#R.<X> = PolynomialRing(RR)
#S.<x> = R.quotient_ring(X^(2^n) + 1)
#u = S.random_element()
#v = S.random_element()
#TT  = _negaconvolution_fft(list(u), list(v), n)
#w_1 = S(TT)
#w_2 = u * v
#assert(w_1 == w_2)




#a = FastFourierTransform(4)



