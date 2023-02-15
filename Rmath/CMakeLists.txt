add_compile_definitions(MATHLIB_STANDALONE)
add_library(Rmath
        mlutils.c fmax2.c fmin2.c
        chebyshev.c lgammacor.c gammalims.c stirlerr.c
        bd0.c gamma.c lgamma.c beta.c lbeta.c cospi.c
        bessel_i.c bessel_j.c bessel_k.c bessel_y.c
        choose.c snorm.c sexp.c dgamma.c sunif.c
        pgamma.c qgamma.c rgamma.c dbeta.c
        pbeta.c qbeta.c rbeta.c dunif.c
        punif.c qunif.c runif.c dnorm.c pnorm.c qnorm.c rnorm.c
        df.c pf.c qf.c rf.c pt.c qt.c rt.c dnt.c dchisq.c pchisq.c qchisq.c
        rchisq.c rnchisq.c dbinom.c pbinom.c qbinom.c rbinom.c rmultinom.c
        dexp.c pexp.c qexp.c rexp.c dhyper.c phyper.c qhyper.c rhyper.c
        dpois.c ppois.c qpois.c rpois.c
        dnchisq.c pnchisq.c qnchisq.c dnbeta.c pnbeta.c qnbeta.c
        pnt.c ptukey.c qtukey.c toms708.c)
target_include_directories(Rmath PUBLIC include)
