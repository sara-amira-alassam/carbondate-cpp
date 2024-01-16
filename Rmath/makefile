all: ../Rmath.a

# Give the names of the required files
OBJ := bd0.o bessel_y.o cospi.o dexp.o dnbeta.o dnt.o fmin2.o lgamma.o pbinom.o pgamma.o pnchisq.o pt.o qbinom.o qgamma.o qnorm.o qunif.o rexp.o rmultinom.o rt.o stirlerr.o
OBJ := $(OBJ) bessel_i.o beta.o dbeta.o df.o dnbinom.o dpois.o gamma.o lgammacor.o pchisq.o phyper.o pnorm.o ptukey.o qchisq.o qhyper.o qpois.o rbeta.o rf.o rnchisq.o runif.o sunif.o
OBJ := $(OBJ) bessel_j.o chebyshev.o dbinom.o dgamma.o dnchisq.o dunif.o gammalims.o mlutils.o pexp.o pnbeta.o pnt.o punif.o qexp.o qnbeta.o qt.o rbinom.o rgamma.o rnorm.o sexp.o toms708.o
OBJ := $(OBJ) bessel_k.o choose.o dchisq.o dhyper.o dnorm.o fmax2.o lbeta.o pbeta.o pf.o pnbinom.o ppois.o qbeta.o qf.o qnchisq.o qtukey.o rchisq.o rhyper.o rpois.o snorm.o

.NOTPARALLEL ../Rmath.a: $(OBJ)
	ar rc $@ $(OBJ)

%.o: %.c include/Rmath.h
	gcc -c -O2 $< -Iinclude -Iinclude/R_ext -DMATHLIB_STANDALONE

clean:
	-rm $(OBJ) ../Rmath.a
