scaled_exp(x,m,beta) = exp( (m-x)/beta)
gumbel(x,m,beta) = scaled_exp(x,m,beta)* exp(-scaled_exp(x,m,beta)) /beta
lognormal(x,mu,sigma) = exp(-(log(x)-mu)**2/(2*sigma**2))/(x*sigma*sqrt(2*pi))

gumbel_P(x,m,beta) = 1-exp(-scaled_exp(x,m,beta))
# Note: gumbel_P may have rounding error problems as exp(-epsilon) is
#	very close to 1.  
# We can approximate it in these cases better using exp(-epsilon) approx 1-epsilon.
# (with double precision arithmetic, we have enough accuracy it
# doesn't matter which way we do it for this problem)
gumbel_P_approx(x,m,beta)= scaled_exp(x,m,beta)


lognormal_P(x,mu,sigma) = 0.5*erfc( (log(x)-mu)/(sigma*sqrt(2)))

set ylabel "probability"
set xlabel "ORF length (codons)"


set title "random codon model"
num_prot = 2994

unset logscale x
set logscale y
m=50; beta=30;
# fit  gumbel(x,m,beta) 'hist1' using 1:3 via m,beta
# fit  gumbel_P(x,m,beta) 'hist1' using 1:4 via m,beta
fit  log(gumbel_P(x,m,beta)) 'hist1' using 1:(log($4)) via m,beta
plot [10:400] gumbel(x,m,beta), 'hist1' using 1:3

print gumbel_P(388, m, beta), gumbel_P_approx(388, m, beta)
print "E-value=", num_prot*gumbel_P(388, m, beta)


pause -1 "press return to continue"

set ylabel "E-value"
plot [10:400] gumbel_P(x,m,beta), 'hist1' using 1:4

pause -1 "press return to continue"


mu=4; sigma=0.5;
# fitting to the log of the P-value seems most robust
# fit lognormal(x,mu,sigma) 'hist1' using 1:3 via mu,sigma
# fit lognormal_P(x,mu,sigma) 'hist1' using 1:4 via mu,sigma
fit log(lognormal_P(x,mu,sigma)) 'hist1' using 1:(log($4)) via mu,sigma

set logscale xy
set ylabel "probability"
plot [10:400] lognormal(x,mu,sigma), 'hist1' using 1:3

print "p-value=",lognormal_P(388, mu, sigma), "  E-value=", num_prot*lognormal_P(388, mu, sigma)

pause -1 "press return to continue"

set ylabel "E-value"
plot [10:400] lognormal_P(x,mu,sigma), 'hist1' using 1:4



