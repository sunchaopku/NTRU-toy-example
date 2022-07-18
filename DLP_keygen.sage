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



q = 2
n = 4
m = 8
while True:
        q = next_prime(q)
        if (q % m) == 1:
            break
#K = GF(q)
#w = K(1).nth_root(m)


R = PolynomialRing(ZZ, 'x')
x = R.gen()
P = R.quotient(x^n + 1, 'y')
f = zpoly_to_qpoly(n, (P.random_element()))
g = zpoly_to_qpoly(n, P.random_element())
Mf = poly_to_matrix(n, f)
Mg = poly_to_matrix(n, g)
f1 = zpoly_to_qpoly(n, matrix_to_poly(n, Mf))
f_ = zpoly_to_qpoly(n, matrix_to_poly(n, Mf.T))
g_ = zpoly_to_qpoly(n, matrix_to_poly(n, Mg.T))
fq = zpoly_to_qpoly(n, f)

T = (f * f_ + g * g_)
t1 = q * f_/T 

xgcd(P_to_R(n, f),P_to_R(n, g))


def DLPkeygen(n, q):
	sigma = 1.17 * sqrt(q/(2 * n))
	D = DiscreteGaussianDistributionIntegerSampler(sigma, c = 0)

	P = PolynomialRing(QQ, 'x')
	x = R.gen()
	P = R.quotient(x^n + 1, 'y')
	while 1:
		f_vec = []
		g_vec = []
		for i in range(n):
			f_vec.append(D())
			g_vec.append(D())
		f = zpoly_to_qpoly(n, P(f_vec))
		g = zpoly_to_qpoly(n, P(g_vec))
		len1 = 0
		for i in range(n):
			len1 = len1 + f[i] * f[i] + g[i] * g[i]
		len1 = sqrt(len1)
		Mf = poly_to_matrix(n, f)
		Mg = poly_to_matrix(n, g)
		f_ = zpoly_to_qpoly(n, matrix_to_poly(n, Mf.T))
		g_ = zpoly_to_qpoly(n, matrix_to_poly(n, Mg.T))
		T = f * f_ + g * g_ 
		t1 = q * f_  / T 
		t2 = q * g / T
		len2 = 0
		for i in range(n):
			len2 = len2 + t1[i]^2 + t2[i]^2
		len2 = sqrt(len2)
		if len1 > 1.17 * sqrt(q) or len2 > 1.17 * sqrt(q):
			continue
		tempf = xgcd(P_to_R(n,f), x^n + 1)
		
		Rf_temp = tempf[1]
		RN_temp = tempf[2]
		Rf = lcmf(n, Rf_temp, RN_temp)
		rau_f = Rf_temp * Rf 

		tempg = xgcd(P_to_R(n,g), x^n + 1)
		
		Rg_temp = tempg[1]
		RN_temp = tempg[2]
		Rg = lcmf(n, Rg_temp, RN_temp)
		rau_g = Rg_temp * Rg 

		print("rau_f is ", rau_f)
		print("rau_f * f is ", rau_f * zpoly_to_qpoly(n, f))
		print("Rf is ", Rf)

		print("rau_g is ", rau_g)
		print("rau_g * g is ", rau_g * zpoly_to_qpoly(n ,g))
		print("Rg is ", Rg)
		if gcd(Rf, Rg) != 1 or gcd( Rf, q) != 1:
			continue
		u = xgcd(Rf, Rg)[1]
		v = xgcd(Rf, Rg)[2]
		print("u *Rf + v* Rg = ", u *Rf + v* Rg)
		F = q * v * rau_g
		G = -q * u * rau_f
		print("f * G - g * F = ", f * G - g * F )
		print("q is ", q)
		k = (F * f_ + G * g_) / T 
		print("k is ", k)
		k = round_f(n, k)
		print("after rounding k is ", k)
		F = F - k * f 
		G = G - k * g
		print("f * G - g * F = ", f * G - g * F )
		print("q is ", q)
		break


q = 2
n = 4
m = 8
while True:
        q = next_prime(q)
        if (q % m) == 1:
            break

DLPkeygen(n, q )







