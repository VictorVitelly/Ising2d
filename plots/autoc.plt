set terminal qt size 850,550
set lmargin 24
set xlabel 't' font ',22'
set ylabel 'C_H (t)' font ',22' rotate by 0 offset -10,0
set title 'Autocorrelation, L=64' font ',24'
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
TT[i]=2.0+0.2*i
}

#8 Metropolis
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
#8 Cluster
#fitfrom[1]=2
#fitfrom[2]=3
#fitfrom[3]=6
#fitfrom[4]=3
#fitto[1]=10
#fitto[2]=16
#fitto[3]=16
#fitto[4]=8
#16 Metropolis
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
#16 Cluster
#fitfrom[1]=3
#fitfrom[2]=4
#fitfrom[3]=3
#fitfrom[4]=2
#fitto[1]=13
#fitto[2]=16
#fitto[3]=12
#fitto[4]=18
#32
#fitfrom[1]=16
#fitfrom[2]=16
#fitfrom[3]=16
#fitfrom[4]=3
#fitfrom[5]=3
#fitto[1]=50
#fitto[2]=50
#fitto[3]=50
#fitto[4]=8
#fitto[5]=8
#32 cluster
#fitfrom[1]=2
#fitfrom[2]=3
#fitfrom[3]=1
#fitfrom[4]=1
#fitto[1]=23
#fitto[2]=21
#fitto[3]=10
#fitto[4]=29

#64
#fitfrom[1]=15
#fitfrom[2]=15
#fitfrom[3]=4
#fitfrom[4]=3
#fitfrom[5]=3
#fitto[1]=50
#fitto[2]=50
#fitto[3]=14
#fitto[4]=10
#fitto[5]=10
#64 Cluster
fitfrom[1]=1
fitfrom[2]=1
fitfrom[3]=0
fitfrom[4]=0
fitto[1]=15
fitto[2]=14
fitto[3]=25
fitto[4]=30


f(x,A,tau)=A*exp(-x/tau)
#g(x)=C*exp(-x/tau2)
array tA[5]
array tAerr[5]
array ttau[5]
array ttauerr[5]
array chi2[5]
tau=2
A=0.01

do for [i=1:4] {
fit f(x,A,tau) '../data/autocC/autocorr64.dat' u 1:column1[i]:column2[i] every ::fitfrom[i]::fitto[i] via A, tau
tA[i]=A
tAerr[i]=A_err
ttau[i]=tau
ttauerr[i]=tau_err
chi2[i]=(FIT_STDFIT*FIT_STDFIT)
}

set xrange[0:20]
set yrange[0.000000001:1]

plot for [i=1:4] '../data/autocC/autocorr64.dat' u 1:column1[i]:column2[i] w errorbars notitle linestyle i lw 2, for [i=1:4] f(x,tA[i],ttau[i]) title sprintf('T=%.1f, τ=%.3f±%.3f, χ^2/dof=%.2f',TT[i],ttau[i],ttauerr[i],chi2[i]) linestyle i

pause -1
