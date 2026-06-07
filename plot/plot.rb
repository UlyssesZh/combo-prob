#!/usr/bin/env ruby

ENV['PYTHONPATH'] = Dir.glob(File.expand_path '../../lib/python3.*/site-packages', `readlink -f #{`which python`.strip}`.strip).first

require 'matplotlib/pyplot'

PLT = Matplotlib::Pyplot
PLT.rcParams.update({
	'text.usetex' => true,
	'font.size' => 10,
	'font.family' => 'lmodern',
	'text.latex.preamble' => '\usepackage{lmodern}'
})

def combo_pc m
	(1..m).each_with_object [[[1]]] do |n, dp|
		dp[n] = Array.new(n+1) { Array.new n+1, 0 }

		# dp[n][0] = (1-Y)**n
		0.upto(n-1) { dp[n][0][_1] = dp[n-1][0][_1] } # will be multiplied by 1-Y later

		1.upto n/2-1 do |k|
			# dp[n][k] = (1-Y) * (Y**k * (0..k).sum { |j| dp[n-k-1][j] } + (0..k-1).sum { |r| Y**r * dp[n-r-1][k] })
			0.upto(k) { |j| 0.upto(n-k-1) { dp[n][k][_1+k] += dp[n-k-1][j][_1] } }
			0.upto(k-1) { |r| 0.upto(n-r-1) { dp[n][k][_1+r] += dp[n-r-1][k][_1] } }
		end

		if n % 2 == 1
			k = n/2
			# dp[n][k] = (1-Y) * (Y**k + (0..k-1).sum { |r| Y**r * dp[n-r-1][k] })
			dp[n][k][k] = 1
			0.upto(k-1) { |r| 0.upto(n-r-1) { dp[n][k][_1+r] += dp[n-r-1][k][_1] } }
		end

		((n+1)/2).upto n-1 do |k|
			# dp[n][k] = (1-Y) * (Y**k + (0..n-k-1).sum { |r| Y**r * dp[n-r-1][k] })
			dp[n][k][k] = 1
			0.upto(n-k-1) { |r| 0.upto(n-r-1) { dp[n][k][_1+r] += dp[n-r-1][k][_1] } }
		end

		0.upto(n-1) { |k| n.downto(1) { dp[n][k][_1] -= dp[n][k][_1-1] } } # multiply by 1-Y
		
		# dp[n][n] = Y**n
		dp[n][n][n] = 1
	end.last
end

def horner coef, x
	coef.reverse.reduce(0) { |r, c| r * x + c }
end

def inv_coef a
	Array.new a.size do |j|
		binom = nil
		(j..N).sum do |k|
			a[k] * (binom = binom ? binom * k / (k-j) : 1)
		end * (-1)**j
	end
end

N = 50
COEF = combo_pc N
COEF_I = COEF.map { |a| inv_coef a }

def val k, y
	y < 0.5 ? horner(COEF[k], y) : horner(COEF_I[k], 1-y)
end

def plot y, label = '$p=%.3f$'%y
	xs = (0..N).to_a
	ys = xs.map { |k| val k, y }
	PLT.plot xs, ys, label:
end

def plot_finite
	PLT.figure figsize: [4.8, 3.6]
	plot 0.9
	plot 0.92
	plot 0.94
	plot 0.96
	plot 0.97
	plot 0.98
	plot 0.99
	plot 0.995
	plot 0.999
	PLT.ylim 0, 0.07
	PLT.legend
	PLT.xlabel '$k$'
	PLT.ylabel '$P_{%d,k}(p)$'%N
	PLT.savefig 'n50.pdf', transparent: true, format: :pdf, bbox_inches: :tight
end
plot_finite

include Math

class Integer
	def fac
		(1..self).reduce 1, :*
	end
end

def inf_first y, x
	-y**x*log(y)*(2+log(y**(x-1)))
end

def inf_rest y, x
	smax = x == 0.5 ? 3 : (1/x).ceil
	(2...smax).sum do |s|
		logysx1 = log y**(s*x-1)
		y**(s*x) * logysx1**(s-2) * s*(-1)**s/(s-2).fac * (
			1 + logysx1*2/(s-1) + logysx1**2/s/(s-1)
		)
	end * log(y)
end

def plot2 y, label = '$P=%.2f$'%y
	xs = Array.new(101) { |i| i/200.0 } + Array.new(101) { |i| 0.5+i/200.0 } + [1]
	ys1 = Array.new 81 do |i|
		x = i/200.0 + 0.1
		inf_first(y,x) + inf_rest(y,x)
	end
	ys2 = Array.new 101 do |i|
		x = i/200.0 + 0.5
		inf_first(y,x)
	end
	ys = [0]*20 + ys1 + ys2 + [5]
	PLT.plot xs, ys, label:
end

def plot_inf
	PLT.figure figsize: [4.8, 3.6]
	plot2 0.01
	plot2 0.03
	plot2 0.09
	plot2 0.3
	plot2 0.6
	plot2 0.9
	PLT.xlabel '$x$'
	PLT.ylabel '$f(P,x)$'
	PLT.ylim 0, 3
	PLT.legend
	PLT.savefig 'inf.pdf', transparent: true, format: :pdf, bbox_inches: :tight
end
plot_inf

def plot3 y
	y **= 1.0/N
	xs = (0..N).map { |k| k/N.to_f }
	ys = (0..N).map { |k| val(k, y)*N }
	PLT.plot xs, ys
end

def plot_compare
	PLT.figure figsize: [4.8, 3.6]
	plot3 0.01
	plot2 0.01
	plot3 0.2
	plot2 0.2
	plot3 0.5
	plot2 0.5
	plot3 0.9
	plot2 0.9
	PLT.xlabel '$x$'
	PLT.ylabel '$f(P,x)$ or $nP_{N,xn}(P^{1/n})$'
	PLT.ylim 0, 3
	PLT.legend
	PLT.title '$n=%d$'%N
	PLT.savefig 'compare.pdf', transparent: true, format: :pdf, bbox_inches: :tight
end
plot_compare
