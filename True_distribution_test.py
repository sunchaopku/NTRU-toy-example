#q: the modulus
#n: degree of the number field
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
import matplotlib.pyplot as plt
import sys

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
def DLPkeygen(alpha, n_, q_, upper_bound, fg_vec, full_idx):
	#global t1,t2,t3,t4,t5
	n = n_
	q = q_
	cnt = 0
	#def DLPkeygen(n, q):
	sigma = alpha * sqrt(q/(2 * n))

	round_upper_bound = floor(sqrt(upper_bound))

	if full_idx > 0:
		for xxxx in range(0, round_upper_bound):
			fg_vec.append(xxx)
			DLPkeygen(alpha, n_, q_, upper_bound - xxxx^2, fg_vec, full_idx - 1)



	D = DiscreteGaussianDistributionIntegerSampler(sigma, c = 0)
	x_axis = []
	y_axis = []
	for kkk in range(12000, 16823):
		x_axis.append(kkk)
		y_axis.append(0)

	R = PolynomialRing(QQ, 'x')
	x = R.gen()
	P = R.quotient(x^n + 1, 'y')
	#while 1:
	min_quality = 2
	#upper_bound = round(alpha * sqrt(q)) + 1

	while 1:
		#	cnt = cnt + 1
			#print("current repetition number is ", cnt) 
			f_vec = []
			g_vec = []
			for i in range(n):
				f_vec.append(D())
				g_vec.append(D())
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
				continue
			else:
				return len1^2
				#y_axis[len1^2 - 12000] = y_axis[len1^2- 12000] + 1
				#print(len1^2 - 12000)
			#	len1_quality = RR(len1/sqrt(q))
			#	min_quality = min(min_quality, quality)
			#	print("alpha is ", alpha)
			#	print("total repetition number is ", cnt) 
			#	print("quality is ", quality)
			#	print("current min quality is ", min_quality)
				#return

	#list_plot(y_axis, x_axis)
	#(list_plot(list(zip(x_axis,y_axis)))).show()
	#print("alpha is ", alpha)
	#print("total repetition number is ", cnt) 
#	plt.rcParams["text.usetex"] =True
#	plt.rcParams["mathtext.fontset"] = "cm"
# plotting the points 
#	plt.plot(x_axis, y_axis, label = 'test')
#	plt.plot(x1, y2, label = 'experimental')
	# naming the x axis
#	plt.xlabel(r'rounding error: $\epsilon$')
	# naming the y axis
#	plt.ylabel('success rate')
	  
	# giving a title to my graph
#	plt.title(r'$q = 12289, \alpha = 1.17, d = 512$ ') #d = 512 alpha = 1.17
	#plt.title(r'$q = 12289, \alpha = 1.63, d = 1024$ ')
	 #pad_inches = 0 
	#plt.axis('off')
#	plt.legend()
#	plt.savefig('test.pdf')

	# function to show the plot
	#plt.show()
global t1,t2,t3,t4,t5
t1 = 0 #   13000
t2 = 0 # 13000 - 14000
t3 = 0 # 14000 - 15000
t4 = 0 # 15000 - 16000
t5 = 0 # 16000 - 16823
q = 11

for xxx in range(10000):
	print("Falcon round ", xxx)
	res = DLPkeygen(alpha = 1.17, n_ = 32, q_ = 11)
	if res < 13000:
		t1 = t1 + 1
	elif res < 14000:
		t2 = t2 + 1
	elif res < 15000:
		t3 = t3 + 1
	elif res < 16000:
		t4 = t4 + 1
	else:
		t5 = t5 + 1
print("t1 ", t1)
print("t2 ", t2)
print("t3 ", t3)
print("t4 ", t4)
print("t5 ", t5)



