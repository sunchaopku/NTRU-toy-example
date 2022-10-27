from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler

import numpy as np
import itertools
from datetime import datetime

params = {  512: {'mK': 1024, 'thr': 2.04, 'q': 12289},
            648: {'mK': 1944, 'thr': 2.13, 'q': 3889 },
            768: {'mK': 2304, 'thr': 2.20, 'q': 18433},
            864: {'mK': 2592, 'thr': 2.25, 'q': 10369},
            972: {'mK': 2916, 'thr': 2.30, 'q': 17497},
           1024: {'mK': 2048, 'thr': 2.33, 'q': 12289} }

param = params[512]

mK  = param['mK']
thr0= param['thr']
q   = param['q']

d = euler_phi(mK)
zmstar = [a for a in range(mK) if gcd(a,mK)==1]
K.<z> = CyclotomicField(mK) # WARNING: degree(K) = d

def permmK(a):
    def rank(xs):
        # double argsort trick
        return xs.argsort().argsort()
    assert(gcd(a,mK)==1)
    return rank(np.array([(a*x)%mK for x in zmstar]))

gc = np.vstack([permmK(a) for a in zmstar[:d//2]])

def selecttest(m,thr=thr0,sigratio=1.17,trials=50):
    l=[]
    for i in range(trials):
        now = datetime.now()
        print("%s: test %d/%d" % (now.strftime("%H:%M:%S"),i+1,trials))
        l.append(selectfg(sigratio=sigratio,m=m,thr=thr)[0])
    return sorted(l)

def Kfft(x):
    return np.fft.fft(vector(RDF,x.list()),mK)[zmstar]

def qual(f,g):
    z = np.absolute(Kfft(f))^2 + np.absolute(Kfft(g))^2
    return sqrt(max(np.max(z)/q, q/np.min(z)))

def qualfft(f,g):
    z = np.absolute(f)^2 + np.absolute(g)^2
    return sqrt(max(np.max(z)/q, q/np.min(z)))

def ihft(f_hft):
    f_fft = np.hstack([f_hft, np.conjugate(f_hft)[::-1]])

    f_fftext = np.zeros(mK, dtype=complex)
    f_fftext[zmstar] = f_fft

    f = 2*np.real(np.fft.ifft(f_fftext)[:d])
    return K(list(np.round(f)))

def randround(a):
    a0 = np.round(a)
    rs = np.random.binomial(n=1, p=np.absolute(a-a0))
    r0 = a0 + rs*np.sign(a-a0)
    return r0

def ihftrand(f_hft):
    f_fft = np.hstack([f_hft, np.conjugate(f_hft)[::-1]])

    f_fftext = np.zeros(mK, dtype=complex)
    f_fftext[zmstar] = f_fft

    f = 2*np.real(np.fft.ifft(f_fftext)[:d])
    return K(list(randround(f)))

def ihftgauss(f_hft):
    f_fft = np.hstack([f_hft, np.conjugate(f_hft)[::-1]])

    f_fftext = np.zeros(mK, dtype=complex)
    f_fftext[zmstar] = f_fft

    f = 2*np.real(np.fft.ifft(f_fftext)[:d])

    eps = 1/(4*sqrt(256*2^64))
    smoothing = RDF(1/pi * sqrt(log(2*d*(1+1/eps))/2))
    fgauss = [DiscreteGaussianDistributionIntegerSampler(sigma=smoothing,
        c=x)() for x in list(f)]
    print(np.array(fgauss[:10])-f[:10])
    return K(fgauss)


def distihft(r):
    fhat_hft_rho2  = np.random.uniform(high=q*r^2,   size=d//2)
    fhat_hft_theta = np.random.uniform(high=2*np.pi, size=d//2)
    fhat_hft = np.sqrt(fhat_hft_rho2) * np.exp(1.0j*fhat_hft_theta)
    f = ihft(fhat_hft)

    f_fft = Kfft(f)
    #u_hft = (np.absolute(f_fft)^2)[:d//2]
    #return np.sqrt(np.max(u_hft/q))
    
    d_hft = np.absolute(f_fft[:d//2]-fhat_hft)^2
    return np.sum(d_hft)*2

def selectfgunif(sigratio=1.17,maxtrials=100,thr=thr0,af=None,ag=None,verbose=False):
    trial = 1

    if ag is None:
        ag = thr
    if af is None:
        af = ag

    while True:
        if verbose:
            print("Trial %d" % trial)
            if trial > maxtrials:
                print("Failure")
                return None

        fhat_hft_rho2  = np.random.uniform(high=q*af^2,  size=d//2)
        fhat_hft_theta = np.random.uniform(high=2*np.pi, size=d//2)
        fhat_hft = np.sqrt(fhat_hft_rho2) * np.exp(1.0j*fhat_hft_theta)

        f = ihft(fhat_hft)

        f_fft = Kfft(f)
        u_hft = (np.absolute(f_fft)^2)[:d//2]

        if np.max(u_hft) >= q*ag^2:
            trial = trial + 1
            continue

        a = np.maximum(q/ag^2 - u_hft, 0)
        b = q*ag^2 - u_hft
        
        ghat_hft_rho2  = np.random.uniform(low=a, high=b)
        ghat_hft_theta = np.random.uniform(high=2*np.pi, size=d//2)
        ghat_hft = np.sqrt(ghat_hft_rho2) * np.exp(1.0j*ghat_hft_theta)

        g = ihft(ghat_hft)
        qfg = qual(f,g)
        print("Quality (f,ghat) = %f" % qualfft(f_fft[:d//2],ghat_hft))
        print("Quality (f,g)    = %f" % qfg)

        if qfg <= thr:
            break
        trial = trial + 1

    if verbose:
        print("Success!")
        print("Achieved quality = %f" % qual(f,g))
    return (f,g)

def selectfgunifsame(maxtrials=100,thr=thr0,errsz=None,errfact=2,verbose=False):
    trial = 1
    suc_cnt = 0
    for xxx in range(1000):
       # if verbose:
        #    print("Trial %d" % trial)
         #   if trial > maxtrials:
          #      print("Failure")
           #     return None

        fghat_hft_rho2 = np.random.uniform(low=(sqrt(q)/thr + errsz)^2,
                                          high=(sqrt(q)*thr - errsz)^2,
                                          size=d//2)
        fghat_hft_theta= np.random.uniform(high=2*np.pi, size=d//2)

        fhat_hft_rho   = np.sqrt(fghat_hft_rho2) * np.cos(fghat_hft_theta)
        ghat_hft_rho   = np.sqrt(fghat_hft_rho2) * np.sin(fghat_hft_theta)

        fhat_hft_theta = np.random.uniform(high=2*np.pi, size=d//2)
        ghat_hft_theta = np.random.uniform(high=2*np.pi, size=d//2)

        fhat_hft = fhat_hft_rho * np.exp(1.0j*fhat_hft_theta)
        ghat_hft = ghat_hft_rho * np.exp(1.0j*ghat_hft_theta)

        f = ihft(fhat_hft)
      #  print("f ", f)
        g = ihft(ghat_hft)

        f_fft = Kfft(f)
        g_fft = Kfft(g)

        qfg = qual(f,g)
    #    print("Quality (fhat,ghat) = %f" % qualfft(fhat_hft,ghat_hft))
     #   print("Quality (f,g)       = %f" % qfg)

        if qfg <= thr:
            suc_cnt = suc_cnt + 1
           # print("total trial ", trial)
           # break
       # trial = trial + 1
        #if trial >= 1000:
         #   print("Failure")
          #  return
    return suc_cnt
    #if verbose:
     #   print("Success!")
      #  print("Achieved quality = %f" % qual(f,g))
   # return (f,g)
def selectfg(sigratio=1.17,maxtrials=5000,thr=thr0,af=None,ag=None,verbose=False):
    sigma0 = RDF(sigratio*sqrt(q/(2*d)))

    eps = 1/(4*sqrt(256*2^64))
    smoothing = 1/pi * sqrt(log(2*d*(1+1/eps))/2)
    if verbose:
        print("sigma_0 = %f" % sigma0)
        print("eta_eps = %f" % smoothing)

    #D0 = DiscreteGaussianDistributionIntegerSampler(sigma=sigma0)

    trial = 1

    if ag is None:
        ag = thr
    if af is None:
        af = ag

    while True:
        if verbose:
            print("Trial %d" % trial)
            if trial > maxtrials:
                print("Failure")
                return None

        fhat_hft = np.zeros(d//2, dtype=complex)
        w = np.arange(d//2),
        while w[0].size > 0:
            fhat_hft[w] = np.random.normal(scale=sqrt(q)/2,size=(w[0].size,2)).view(complex).reshape(w[0].size)
            u_hft = np.absolute(fhat_hft)^2
            w = np.where(u_hft >= q*af^2)

        f = ihft(fhat_hft)

        f_fft = Kfft(f)
        u_hft = (np.absolute(f_fft)^2)[:d//2]

        if np.max(u_hft) >= q*ag^2:
            trial = trial + 1
            continue

        a = q/ag^2 - u_hft
        b = q*ag^2 - u_hft
        aa= np.where(a>0, a, b/2)
        s = np.where(a>0, np.sqrt((b-aa)/np.log(b/aa)), 0)
        s = np.stack([s/np.sqrt(2),s/np.sqrt(2)]).T

        ghat_hft = np.zeros(d//2, dtype=complex)
        w = np.where(a>0)
        while w[0].size > 0:
            #print(w[0].size)
            ghat_hft[w] = np.random.normal(scale=s[w]).view(complex).reshape(w[0].size)
            z_hft = u_hft + np.absolute(ghat_hft)^2
            w     = np.where((z_hft <= q/ag^2) | (z_hft >= q*ag^2))

        g = ihft(ghat_hft)
        qfg = qual(f,g)
        print("Quality (f,ghat) = %f" % qualfft(f_fft[:d//2],ghat_hft))
        print("Quality (f,g)    = %f" % qfg)

        if qfg <= thr:
            break
        trial = trial + 1

    if verbose:
        print("Success!")
        print("Achieved quality = %f" % qual(f,g))
    return (f,g)


def computeFG(f,g):
    # outputs F,G such that fG-gF = q and (f,g), (F,G) generates the
    #  orthogonal lattice of (f/g mod q, -1) as long as g is invertible mod q
    Nf = f.norm()
    Ng = g.norm()
    s = Nf/f
    t = Ng/g

    d,u,v = xgcd(Nf,Ng) # u*Nf + v*Nd = d in ZZ
    F_big = -q*t*v/d
    G_big = q*s*u/d

    assert f*G_big-g*F_big == q
    ktemp = (f.conjugate()*F_big+g.conjugate()*G_big)/(f*f.conjugate() + g*g.conjugate())
    k = K( [ e.round() for e in ktemp.list() ] )

    F = F_big - k*f
    G = G_big - k*g

    assert f*G-g*F == q
    return (F,G)

def keygen(sigratio=1.17,m=16,ng=d//2,thr=thr0,verbose=False,early=False,alg='qf'):
    itfg = 0
    if verbose:
        print("Generating (f,g) pair.")
    while True:
        itfg += 1
        if verbose:
            print("...iteration %d" % itfg)
        try:
            s,f,g = selectfg(sigratio,m,ng,thr,False,early,alg)
            if s <= thr:
                if verbose:
                    print("...success!")
                break
            else:
                if verbose:
                    print("...failed (s = %f > %f)" % (s,thr))
        except RuntimeError:
            if verbose:
                print("...failed (not found)")
            pass

    if verbose:
        print("Generating (F,G) pair.")
    F,G = computeFG(f,g)
    if verbose:
        print("Key generation complete!")
    return f,g,F,G

factor = 0.001
for i in range(1000):
    print("factor is ", factor)
    print(selectfgunifsame(thr=1.17,errsz= factor * sqrt(12289)))
    factor = factor + 0.001
# vim: ft=python ts=4 expandtab
