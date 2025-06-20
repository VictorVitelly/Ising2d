set terminal qt size 850,550
set lmargin 24
set xlabel 't' font ',22'
set ylabel 'C_H (t)' font ',22' rotate by 0 offset -10,0
set title 'Autocorrelation, L=32' font ',24'
set key font ',18'
set xtics font ", 18"
set ytics font ", 18"
set logscale y
set rmargin 5
set grid x,y
set format y "10^{%L}"

array column1[5]
array column2[5]
array TT[5]
array fitfrom[5]
array fitto[5]

do for [i=1:5] {
set style line i pt 1
column1[i]=2*i
column2[i]=2*i+1
TT[i]=-2.0-0.2*i
}

#8
#fitfrom[1]=5
#fitfrom[2]=4
#fitfrom[3]=5
#fitfrom[4]=3
#fitfrom[5]=4
#fitto[1]=17
#fitto[2]=18
#fitto[3]=14
#fitto[4]=14
#fitto[5]=14
#16
#fitfrom[1]=8
#fitfrom[2]=8
#fitfrom[3]=7
#fitfrom[4]=4
#fitfrom[5]=4
#fitto[1]=30
#fitto[2]=30
#fitto[3]=30
#fitto[4]=14
#fitto[5]=14
#32
fitfrom[1]=16
fitfrom[2]=16
fitfrom[3]=16
fitfrom[4]=3
fitfrom[5]=3
fitto[1]=50
fitto[2]=50
fitto[3]=50
fitto[4]=8
fitto[5]=8
#64

f(x,A,tau)=A*exp(-x/tau)
#g(x)=C*exp(-x/tau2)
array tA[5]
array tAerr[5]
array ttau[5]
array ttauerr[5]
array chi2[5]
tau=2
A=0.01

do for [i=1:5] {
#fit f(x,A,tau) '../data/autocorr.dat' u 1:2:3 every ::20::139 via A, tau
fit f(x,A,tau) '../data/autocorr32.dat' u 1:column1[i]:column2[i] every ::fitfrom[i]::fitto[i] via A, tau
tA[i]=A
tAerr[i]=A_err
ttau[i]=tau
ttauerr[i]=tau_err
chi2[i]=(FIT_STDFIT*FIT_STDFIT)
}

#set xrange[0:20]
set yrange[0.00000001:10]

plot for [i=1:5] '../data/autocorr32.dat' u 1:column1[i]:column2[i] w errorbars notitle linestyle i lw 2, for [i=1:5] f(x,tA[i],ttau[i]) title sprintf('T=%.2f, τ=%.2f±%.2f, χ^2/dof=%.2f',TT[i],ttau[i],ttauerr[i],chi2[i]) linestyle i

pause -1
