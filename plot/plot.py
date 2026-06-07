#!/usr/bin/env python

from scipy.special import expi
from numpy import exp, log, linspace, euler_gamma, insert, zeros_like, arange, ceil, where, append
from matplotlib import pyplot as plt

plt.rcParams.update({
	'text.usetex': True,
	'font.size': 10,
	'font.family': 'lmodern',
	'text.latex.preamble': r'\usepackage{lmodern}'
})

def li(x):
	return expi(log(x))

def repow(x, k):
	return where(x>0, x**k, 0.0)

p = linspace(0.001, 0.999, 999)
mu1 = (li(p) - log(-log(p)) - euler_gamma) / log(p)
p = [0, *p, 1]
mu1 = [0, *mu1, 1]

plt.figure(figsize=(6,2))
plt.plot(p, mu1)
plt.ylim(0, 1)
plt.xlabel('$P$')
plt.ylabel(r'$\mu_1(P)$')
plt.savefig('mu1.pdf', bbox_inches='tight', transparent=True, format='pdf')

def x_and_F(P, upper_x = 1.3):
	lnP = log(P)
	step = 0.005
	hard_zero = 0.02 # F(P,x<hard_zero) is hardcoded to zero
	hard_zero_i = int(hard_zero / step)
	x = arange(0, 1, step)
	if upper_x > 1:
		x = append(x, [upper_x])
	F = zeros_like(x)
	s = 1
	s_fac = 1
	while True:
		x1 = arange(hard_zero, 1/s, step)
		if len(x1) == 0:
			break
		s_fac *= s
		F[hard_zero_i:hard_zero_i+len(x1)] += P**(s*x1) * ((1-s*x1)*lnP)**s / s_fac / s
		s += 1
	F -= li(P) - log(-lnP) - euler_gamma
	F /= lnP
	F += x
	F[:hard_zero_i] = 0
	return x, F

def plot1(p):
	plt.plot(*x_and_F(p), label='$P=%.2f$'%p)

plt.figure(figsize=(4.8,3.6))
plot1(0.01)
plot1(0.03)
plot1(0.09)
plot1(0.3)
plot1(0.6)
plot1(0.9)
xmin, xmax = plt.xlim()
plt.xlim(xmin, 1.2)
plt.ylim(0, 0.6)
plt.legend()
plt.xlabel('$x$')
plt.ylabel('$F(P,x)$')
plt.savefig('F.pdf', bbox_inches='tight', transparent=True, format='pdf')

def plot_f_k(k):
	if k == 0:
		plt.plot([0, 1, 1], [0, 0, 5], label='$k=0$')
		return
	if k == 1:
		plt.plot([0, 0.5, 0.5, 1], [0, 0, 2, 2], label='$k=1$')
		return
	x = linspace(0, 1, 201)
	f = 0
	fac = k*(k+1)
	for s in range(1, k+2):
		f += fac * repow(1-s*x, k-1)
		fac *= (k+1-s) / -s
	plt.plot(x, f, label=f'$k={k}$')

plt.figure(figsize=(2.8,2.1))
plot_f_k(0)
plot_f_k(1)
plot_f_k(2)
plot_f_k(3)
plot_f_k(4)
plt.ylim(0, 4)
plt.xlabel('$x$')
plt.ylabel('$f_k(x)$')
plt.savefig('f_k.pdf', bbox_inches='tight', transparent=True, format='pdf')

def x_and_F_k(k, upper_x = 1.3):
	x = linspace(0, 1, 201)
	if upper_x > 1:
		x = append(x, [upper_x])
	F = zeros_like(x)
	fac = 1/(k+1)
	for s in range(1, k+2):
		fac *= (k+2-s) / -s
		F += fac/s * (1 - repow(1-s*x, k+1))
	F += x
	return x, F

def plot_F_k(k):
	plt.plot(*x_and_F_k(k), label=f'$k={k}$')

plt.figure(figsize=(2.8,2.1))
plot_F_k(0)
plot_F_k(1)
plot_F_k(2)
plot_F_k(3)
plot_F_k(4)
xmin, xmax = plt.xlim()
plt.xlim(xmin, 1.2)
plt.ylim(0, 0.6)
plt.legend()
plt.xlabel('$x$')
plt.ylabel('$F_k(x)$')
plt.savefig('BigF_k.pdf', bbox_inches='tight', transparent=True, format='pdf')

def plot2(P):
	lmd = -log(P)
	x, F = x_and_F(P, upper_x = (2 + log(lmd))/lmd)
	y = x * lmd - log(lmd)
	plt.plot(y, F * lmd, label=r'$\mathrm e^{-\lambda}=%.4f$' % P)

plt.figure(figsize=(2.8,2.1))
plot2(0.0001)
plot2(0.001)
plot2(0.01)
plot2(0.1)
y = linspace(-2, 2, 201)
plt.plot(y, -expi(-exp(-y)), label=f'Gumbel')
plt.xlim(-2, 2)
plt.ylim(0, 1)
plt.legend()
plt.xlabel('$y$')
plt.ylabel(r'$\lambda F(\mathrm e^{-\lambda},(y+\ln\lambda)/\lambda)$')
plt.savefig('Gumbel.pdf', bbox_inches='tight', transparent=True, format='pdf')

def plot_F_k_2(k):
	x, F = x_and_F_k(k, upper_x = (2 + log(k))/k)
	y = x * k - log(k)
	plt.plot(y, F * k, label=r'$k=%i$' % k)

plt.figure(figsize=(2.8,2.1))
plot_F_k_2(5)
plot_F_k_2(20)
plot_F_k_2(35)
plot_F_k_2(50)
y = linspace(-2, 2, 201)
plt.plot(y, -expi(-exp(-y)), label=f'Gumbel')
plt.xlim(-2, 2)
plt.ylim(0, 1)
plt.legend()
plt.xlabel('$y$')
plt.ylabel(r'$kF_k((y+\ln k)/k)$')
plt.savefig('Gumbel_F_k.pdf', bbox_inches='tight', transparent=True, format='pdf')
