#### Chiang interval 
CIex <- function(age, mx, dth, ns=1000, level=0.95, all.draws = FALSE){
        # number of ages
        m = length(mx)
        # estimated probs
        qx = derive_ex_values(mx, age, full = TRUE)$qx
        # trials for binomial, rounded
        Ntil = round(dth/qx)
        # simulated dth counts
        # from binomial dist
        Y = suppressWarnings(matrix(rbinom(m*ns,
                                           Ntil,
                                           qx),
                                    m,ns))
        QX = Y/Ntil
        
        if (all.draws) {
                exsim = sapply(1:ns, function(i) derive_ex_values_qx(mx, QX[, i], age)[1])
                return(exsim)
        }
        else {
                exsim = sapply(1:ns, function(i) derive_ex_values_qx(mx, QX[, i], age))
                ## confidence interval
                CI = apply(exsim, 1, quantile, probs = c((1-level)/2, 0.5, 1 - (1-level)/2))
                return(CI)
        }
}