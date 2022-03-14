cap mata mata drop ManyIVStats()
cap mata mata drop ManyIVreg_IM()

cap mata mata drop sf_helper_epsilon()
cap mata mata drop sf_helper_sig()
cap mata mata drop sf_helper_annihilator()
cap mata mata drop sf_helper_solve()
cap mata mata drop sf_helper_tsolve()
cap mata mata drop sf_helper_licols()

mata
struct ManyIVStats {
    real scalar F
    real matrix Omega
    real matrix Xi
    real vector Sargan
    real vector CD
}

class ManyIVreg_IM
{
    real rowvector beta
    real matrix se, RFS
    struct ManyIVStats scalar stats
    real scalar n, K, L, F

    real scalar clustered
    real scalar estimatese
    real scalar estimatestats
    real scalar small
    real scalar cons
    real scalar nabsorbed_w
    real scalar nabsorbed_z
    real scalar jive

    void results()
    void print()
    void fit()

    string colvector betaLabels
    string colvector seLabels
}

// y              = `Y'
// T              = `X'
// Z              = `Z'
// W              = `W'
// _estimatese    = `estimatese'
// _estimatestats = `estimatestats'
// _small         = `small'
// _cons          = `cons'
// cluster        = ManyIVreg_Absorb_New(tokens(st_local("cluster")), "`touse'")
// Absorb         = ManyIVreg_Absorb_New(tokens(st_local("absorb")), "`touse'")
// AbsorbIV       = ManyIVreg_Absorb_New(tokens(st_local("absorbiv")), "`touse'")
// AbsorbIV.append(Absorb)

void function ManyIVreg_IM::fit(
    real colvector y,
    real colvector T,
    real matrix Z,
    real matrix W,
    real scalar _estimatese,
    real scalar _estimatestats,
    real scalar _small,
    real scalar _cons,
    class ManyIVreg_Absorb scalar cluster,
    class ManyIVreg_Absorb scalar Absorb,
    class ManyIVreg_Absorb scalar AbsorbIV)
{
    real scalar i, G, qc, skipcons
    real scalar mmin, lamre, Qs, c, Lam11, Lam22, h, Vvalid, Vinvalid
    real matrix Xi
    real vector overid, pvalue, selw
    real colvector sel_i, yp, Tp, yq, Tq
    real colvector wselix, zselix
    real rowvector coll

    real matrix YY, YPY, YMY, ZW, yTZW
    real matrix epsilon, epsilon_i, hatP, hatP_i
    real matrix Wp, Zq, Wq, MW_yT, MWD_yT, MWD_Z, HZ_yT, S, Sp
    real matrix Omre, Omure, Gamma, Sig
    real vector k, ei, sec, DZW, DW, iIDZW, iIDW, a, b
    real vector hatTjive, hatTujive, hatPjive, hatPujive

    estimatese    = _estimatese
    estimatestats = _estimatestats
    small         = _small
    cons          = _cons
    nabsorbed_w   = 0
    nabsorbed_z   = 0
    jive          = ((AbsorbIV.nabsorb <= 2) | AbsorbIV.d_computed) & ((Absorb.nabsorb <= 2) | Absorb.d_computed)

    if ( jive == 0 ) {
        errprintf("jive/ujive will not be computed with more than 2 absorb groups\n")
    }

    if ( (AbsorbIV.nsingledrop < AbsorbIV.nsingletons) | (Absorb.nsingledrop < Absorb.nsingletons) ) {
        errprintf("jive/ujive with singleton absorb groups not implemented; will skip\n")
        jive = 0
    }

    n = rows(Z)
    K = cols(Z)
    L = cols(W)
    betaLabels = ("OLS", "TSLS", "LIML", "MBTSLS", "JIVE", "UJIVE", "RTSLS")'
    seLabels   = ("Homoskedastic", "Heteroscedastic", "Cluster", "ManyIV", "ManyIV")'

    // Note: This is not really necessary bc I check there is at least one
    // instrument and there can be at most one exogenous variable. Take out?
    if ( (K + AbsorbIV.df) < cols(T) ) {
        errprintf("need at least as many instruments (found %g) as endogenous variables (%g)\n",
                  K + AbsorbIV.df, cols(T))
        error(198)
    }

    // 2. Point estimates
    // Note: These apply Frisch–Waugh–Lovell; project the covariate
    // of interest and the instrument onto W's null space then use
    // univariate formulas.

    if ( Absorb.nabsorb ) {
        Absorb.flagredundant()
        nabsorbed_w = Absorb.df + !cons

        yTZW   = Absorb.hdfe((y, T, Z, W))
        zselix = cols(Z)? 3::(2 + cols(Z)): J(0, 1, 0)
        wselix = cols(W)? (3 + cols(Z))::(2 + cols(Z) + cols(W)): J(0, 1, 0)
        coll   = sf_helper_licols(yTZW[., wselix], Absorb.hdfetol / n)

// TODO: xx does this actually merit error? Perfectly projecting into
// covariates might be something to allow. Otherwise don't omit y, T.
//
//         if ( coll[1] == 0 ) {
//             errprintf("dependent variable collinear with absorb groups\n")
//             error(1234)
//         }
//
//         if ( coll[2] == 0 ) {
//             errprintf("endogenous variable of interest collinear with absorb groups\n")
//             error(1234)
//         }

        yp   = yTZW[., 1]
        Tp   = yTZW[., 2]
        Wp   = cols(W)? select(yTZW[., wselix], coll): W
        L    = cols(Wp) + nabsorbed_w
    }
    else {
        yp   = y
        Tp   = T
        Wp   = W
    }

    if ( AbsorbIV.nabsorb ) {
        AbsorbIV.flagredundant()
        nabsorbed_z = AbsorbIV.df

        skipcons = cons & !AbsorbIV.allhaveskip()
        selw   = cols(W)? (J(1, cols(W)-1, 1), !skipcons): J(1, 0, 0)
        yTZW   = AbsorbIV.hdfe((y, T, Z, select(W, selw)))
        zselix = cols(Z)? 3::(3 + cols(Z) - 1): J(0, 1, 0)
        wselix = any(selw)? (3 + cols(Z))::cols(yTZW): J(0, 1, 0)
        coll   = sf_helper_licols(yTZW[., (zselix \ wselix)], AbsorbIV.hdfetol / n)

// TODO: xx ibid.
        yq   = yTZW[., 1]
        Tq   = yTZW[., 2]
        Zq   = cols(Z)?   select(yTZW[., zselix], coll[zselix :- 2]): Z
        Wq   = any(selw)? select(yTZW[., wselix], coll[wselix :- 2]): select(W, selw)
        K    = cols(Zq) + nabsorbed_z
    }
    else {
        yq   = yp
        Tq   = Tp

        // Account for collinear columns | instrument in RF and FS (only
        // needed of no absorb instruments)
        if ( Absorb.nabsorb ) {
            coll = sf_helper_licols(yTZW[., (zselix \ wselix)], Absorb.hdfetol / n)
            Zq   = cols(Z)? select(yTZW[., zselix], coll[zselix :- 2]): Z
            Wq   = cols(W)? select(yTZW[., wselix], coll[wselix :- 2]): Wp
// TODO: xx ibid.
        }
        else {
            Zq = Z
            Wq = Wp
        }
        K = cols(Zq)
    }

    if ( K < cols(T) ) {
        errprintf("instruments collinear with absorb levels (only %g independent; need %g)\n", K, cols(T))
        error(198)
    }

    if ( cols(Zq) < cols(Z) ) {
        printf("(dropped %g collinear instruments)\n", cols(Z) - cols(Zq))
    }

    MW_yT  = sf_helper_annihilator(Wp, (yp, Tp))          // [y_⊥ T_⊥] = M_W [y T]
    MWD_yT = sf_helper_annihilator(Wq, (yq, Tq))          // [y_⊥ T_⊥] = M_W M_D [y T]
    MWD_Z  = sf_helper_annihilator(Wq, Zq)                // Z_⊥ = M_W M_D Z
    YY     = (MW_yT' * MW_yT)                             // [y_⊥ T_⊥]' [y_⊥ T_⊥] = [y T]' M_W [y T]
    RFS    = sf_helper_solve(MWD_Z, MWD_yT)               // [solve(Z_⊥, y_⊥) solve(Z_⊥, T_⊥)] = Reduced form and First stage
    HZ_yT  = (cols(Zq)? MWD_Z * RFS: 0) :+ MW_yT - MWD_yT // H_{Z_⊥} [y_⊥ T_⊥]
    YPY    = MW_yT' * HZ_yT                               // [y_⊥ T_⊥]' H_{Z_⊥} [y_⊥ T_⊥]
    YMY    = YY - YPY                                     // [y_⊥ T_⊥]' M_{Z_⊥} [y_⊥ T_⊥]

    // 2.1 k-class: OLS, TSLS, LIML, MBTLS
    // Note: These are all coded as
    //
    //     (T_⊥' (I - k W_{Z_⊥}) y_⊥) / (T_⊥' (I - k W_{Z_⊥}) T_⊥)
    //
    // So different values of k give different estimands.
    //
    // - k = 0 -> YY[1, 2]  / YY[2, 2]  = (y_⊥' T_⊥) / (T_⊥' T_⊥)
    // - k = 1 -> YPY[1, 2] / YPY[2, 2] = (y_⊥' H_{Z_⊥} T_⊥) / (T_⊥' H_{Z_⊥} T_⊥)
    // - The other two give liml and mbtsls

    eigensystem(invsym(YMY) * YY, ., ei=.)
    k    = (0, 1, min(Re(ei)), (1 - L/n) / (1 - (K - 1) / n - L/n))
    beta = (YY[1, 2] :- k :* YMY[1, 2]) :/ (YY[2, 2] :- k :* YMY[2, 2])

    // 2.2 JIVE, UJIVE
    if ( jive ) {
        ZW  = (Wq, Zq)
        DZW = rowsum((ZW * invsym(ZW' * ZW)) :* ZW) // D_{Z W} = diag(H_{Z W}) as a vector
        DW  = rowsum((Wp * invsym(Wp' * Wp)) :* Wp) // D_W     = diag(H_W) as a vector

        if ( (Absorb.nabsorb > 0) & (AbsorbIV.nabsorb > 0) ) {
            DW  = DW  :+ Absorb.d_projection()
            DZW = DZW :+ AbsorbIV.d_projection()
        }
        else if ( Absorb.nabsorb > 0 ) {
            DW  = DW  :+ Absorb.d_projection()
            DZW = DZW :+ Absorb.d_projection()
        }
        else if ( AbsorbIV.nabsorb > 0 ) {
            DZW = DZW :+ AbsorbIV.d_projection()
        }

        iIDZW  = 1 :/ edittozero(1 :- DZW, n) // (I - D_{Z W})^{-1} as a vector
        iIDW   = 1 :/ edittozero(1 :- DW, n)  // (I - D_W)^{-1} as a vector

        hatTujive = T :- iIDZW :* sf_helper_annihilator(ZW, Tq) //     (I - (I - D_{Z W})^{-1} M_{Z W}) T
        hatPjive  = sf_helper_annihilator(Wp, hatTujive)        // M_W (I - (I - D_{Z W})^{-1} M_{Z W}) T
        Absorb._hdfe(hatPjive)

        hatTjive  = T :- iIDW :* sf_helper_annihilator(Wp, Tp)  // (I - (I - D_W)^{-1} M_W) T
        hatPujive = hatTujive - hatTjive // (I - D_W)^{-1} M_W T - (I - D_{Z W})^{-1} M_{Z W} T
                                         // = (I - D_W)^{-1} (I - H_W) T - (I - D_{Z W})^{-1} (I - H_{Z W}) T
                                         // = ((I - D_{Z W})^{-1} (H_{Z W} - I) - (I - D_W)^{-1} (H_W - I) T) T
                                         // = (
                                         //     (I - D_{Z W})^{-1} (H_{Z W} - D_{Z W} - (I - D_{Z W})) -
                                         //     (I - D_W)^{-1} (H_W - D_W - (I - D_W))
                                         // ) T
                                         // = ((I - D_{Z W})^{-1} (H_{Z W} - D_{Z W}) - (I - D_W)^{-1} (H_W - D_W)) T

        beta =  (
            beta,
            (hatPjive'  * y) / (hatPjive'  * T),
            (hatPujive' * y) / (hatPujive' * T),
            YPY[1, 1] / YPY[1, 2]
        )
    }
    else {
        beta =  beta, J(1, 2, .), YPY[1, 1] / YPY[1, 2]
    }

    // ------------------
    // 5. Standard Errors
    // ------------------

    Sp = YMY/(n-K-L)  // S_{perp}
    se = J(5, 7, .)
    if ( estimatese ) {

        // -----------------
        // 5.1 Homoscedastic
        // -----------------

        if ( jive ) {
            se[1, 1::6] = sqrt((
                sf_helper_sig(MW_yT, beta[1]) / (MW_yT[.,2]'*MW_yT[.,2]),
                (
                    sf_helper_sig(MW_yT, beta[2]), sf_helper_sig(MW_yT, beta[3]), sf_helper_sig(MW_yT, beta[4])
                ) / YPY[2,2],
                sf_helper_sig(MW_yT, beta[5]) * (hatPjive'  * hatPjive)  / (hatPjive'  * T)^2,
                sf_helper_sig(MW_yT, beta[6]) * (hatPujive' * hatPujive) / (hatPujive' * T)^2
            ))
        }
        else {
            se[1, 1::4] = sqrt((
                sf_helper_sig(MW_yT, beta[1]) / (MW_yT[.,2]'*MW_yT[.,2]),
                (
                    sf_helper_sig(MW_yT, beta[2]), sf_helper_sig(MW_yT, beta[3]), sf_helper_sig(MW_yT, beta[4])
                ) / YPY[2,2]
            ))
        }

        // -------------------
        // 5.2 Heteroscedastic
        // -------------------

        // ols, tsls, liml, mbtsls, jive, ujive
        // Note for asymptotics k -> 1 for the various estimators.
        if ( jive ) {
            hatP    = MW_yT[.,2], J(1, 3, HZ_yT[.,2]), hatPjive, hatPujive
            epsilon = sf_helper_epsilon(MW_yT, beta[1::6])
            se[2, 1::6] = sqrt(colsum((epsilon :* hatP):^2)) :/ (T' * hatP)
        }
        else {
            hatP    = MW_yT[.,2], J(1, 3, HZ_yT[.,2])
            epsilon = sf_helper_epsilon(MW_yT, beta[1::4])
            se[2, 1::4] = sqrt(colsum((epsilon :* hatP):^2)) :/ (T' * hatP)
        }

        // -------------------------
        // 5.4 Cluster, if requested
        // -------------------------

        clustered = cluster.nabsorb
        if ( clustered > 1 ) {
            errprintf("multi-way clustering not implemented; will ignore\n")
            clustered = 0
        }
        else if ( clustered ) {
            sec   = J(1, jive? 6: 4, 0)
            G     = cluster.nlevels[1]
            for(i = 1; i <= G; i++) {
                sel_i     = panelsubmatrix(cluster.index(1), i, cluster.info(1))
                hatP_i    = hatP[sel_i, .]
                epsilon_i = epsilon[sel_i, .]
                sec = sec + colsum(epsilon_i :* hatP_i):^2
            }
            // qc = small? sqrt(G / (G - 1)): 1
            qc = small? sqrt(((n - 1) / (n - L - 1)) * (G / (G - 1))): 1
            if ( jive ) {
                se[3, 1::6] = qc :* sqrt(sec) :/ abs(T' * hatP)
            }
            else {
                se[3, 1::4] = qc :* sqrt(sec) :/ abs(T' * hatP)
            }
        }

        // --------------------
        // 5.3 Many instruments
        // --------------------

        // Notation
        S    = YPY/n
        eigensystem(sf_helper_solve(Sp, S), ., ei=.)
        mmin = min(Re(ei))

        // Hessian of random-effects
        lamre = max(Re(ei)) - K/n
        a     = beta[3] \ 1
        b     = 1 \ -beta[3]
        Omre  = (n-K-L) * Sp/(n-L) + n * (S :- lamre * sf_helper_tsolve((a*a'), sf_helper_tsolve(a', Sp) * a)) / (n-L)
        Qs    = (b' * S * b) / (b' * Omre * b)
        c     = lamre * Qs / ((1-L/n) * (K/n+lamre))

        se[4, 3] = sqrt(
            -b'*Omre*b / (n*lamre) * (lamre+K/n) /
            (Qs*Omre[2,2] - S[2,2] + (c/(1-c)) * Qs / (sf_helper_tsolve(a', Omre) * a))
        )

        // mbtsls, using maximum URE likelihood plug-in estimator
        b = 1 \ -beta[4] // b_mbtsls
        Lam11 = max((0, b' * (S - K/n * Sp) * b))

        if ( mmin > K/n ) {
            Lam22 = S[2, 2] - K/n * Sp[2, 2]
            Omure = Sp
        }
        else {
            Lam22 = lamre / (sf_helper_tsolve(a', Omre) * a)
            Omure = Omre
        }

        Gamma = (1, 0) \ (-beta[4], 1)
        Sig   = Gamma' * Omure * Gamma
        h     = ((1-L/n) * (K-1)/n) / (1 - L/n - (K-1)/n)

        Vvalid   = Sig[1,1] / Lam22 + h * (Sig[1,1] * Sig[2,2] + Sig[1,2]^2) / Lam22^2
        Vinvalid = Vvalid + (Lam11 * Omure[2,2] + Lam11 * Lam22 * n/K) / Lam22^2

        se[4::5, 4] = sqrt((Vvalid \ Vinvalid) / n)

        // Small-sample adjustment, if requested
        qc = small? sqrt(n / (n - L - 1)): 1
        se[(1\2\4\5), .] = qc :* se[(1\2\4\5), .]
    }

    F = YPY[2, 2] / (K * Sp[2,2]) // First-stage F
    if ( estimatestats ) {
        Xi = YPY/n - (K/n) * Sp // % Xi

        overid = J(2, 1, .)
        pvalue = J(2, 1, .)
        if ( K > 1 ) {
            overid[1] = n * mmin / (1 - K/n - L/n + mmin) // n* J_sargan
            pvalue[1] = chi2tail(K - 1, overid[1])        // p-value for Sargan

            overid[2] = n * mmin  // Cragg-Donald
            pvalue[2] = 1 - normal(sqrt((n-K-L)/(n-L)) * invnormal(chi2(K - 1, overid[2])))
        }

        stats.F      = F
        stats.Omega  = Sp
        stats.Xi     = Xi
        stats.Sargan = overid[1], pvalue[1]
        stats.CD     = overid[2], pvalue[2]
    }
}

void function ManyIVreg_IM::results(string scalar bname, string scalar sename)
{

    if ( (bname != "") ) {
        st_matrix(bname, editmissing(beta, 0))
        st_matrixcolstripe(bname, (J(rows(betaLabels), 1, ""), betaLabels))
    }

    if ( (sename != "") ) {
        st_matrix(sename, se)
        st_matrixcolstripe(sename, (J(rows(betaLabels), 1, ""), betaLabels))
        st_matrixrowstripe(sename, (J(rows(seLabels),   1, ""), seLabels))
    }
}

void function ManyIVreg_IM::print()
{
    string scalar note_absw, note_absz
    real scalar j, maxl

    printf("\n")
    note_absw  = nabsorbed_w? " (" + strtrim(sprintf("%21.0fc", nabsorbed_w)) + " absorbed)": ""
    note_absz  = nabsorbed_z? " (" + strtrim(sprintf("%21.0fc", nabsorbed_z)) + " absorbed)": ""

    maxl = max(strlen(betaLabels)) + 1
    if ( estimatese ) {
        printf(sprintf("%%%gs %%9s %%11s\n", maxl), "", "Coef.", clustered? "Cluster": "Hetero")
        printf(sprintf("%%%gs %%9s %%11s\n", maxl), maxl * " ", 5 * "-", (clustered? 7: 6) * "-")
        for(j = 1; j <= length(beta); j++) {
            printf(sprintf("%%%gs %%9.4f %%11s\n", maxl),
                   betaLabels[j], beta[j], "(" + strtrim(sprintf("%9.4f", se[2 :+ clustered, j])) + ")")
        }
    }
    else {
        printf(sprintf("%%%gs %%9s\n", maxl), "", "Coef.")
        printf(sprintf("%%%gs %%9s\n", maxl), maxl * " ", 5 * "-")
        for(j = 1; j <= length(beta); j++) {
            printf(sprintf("%%%gs %%9.4f\n", maxl), betaLabels[j], beta[j])
        }
    }
    printf("\n%s observations, ", strtrim(sprintf("%21.0fc", n)))
    printf("%s instrument%s%s, ", strtrim(sprintf("%21.0fc", K)), (K > 1)? "s": "", note_absz)
    printf("%s covariate%s%s, ", strtrim(sprintf("%21.0fc", L)), (L > 1)? "s": "", note_absw)
    printf("first-stage F = %s\n", strtrim(sprintf("%21.3fc", F)))
}

real matrix function sf_helper_epsilon(real matrix Yp, real rowvector beta)
{
    return(Yp[., 1] :- Yp[., 2] * beta)
}

real scalar function sf_helper_sig(real matrix Yp, real scalar beta)
{
    real colvector e
    e = sf_helper_epsilon(Yp, beta)
    return(e' * e / length(e))
}

real matrix function sf_helper_annihilator(real matrix X, real matrix Y)
{
    return((cols(X) & cols(Y))? (Y - X * sf_helper_solve(X, Y)): Y)
}

real matrix function sf_helper_solve(real matrix X, real matrix Y)
{
    return((cols(X) & cols(Y))? (invsym(cross(X, X)) * cross(X, Y)): J(cols(X), cols(Y), 0))
}

real matrix function sf_helper_tsolve(real matrix X, real matrix Y)
{
    return((invsym(cross(Y', Y')) * cross(Y', X'))')
}

real rowvector function sf_helper_licols(real matrix X, | real scalar tol)
{
    real rowvector p
    real matrix R
    real colvector D
    if ( (cols(X) == 0) | (rows(X) == 0) ) return(J(1, 0, 0))
    if ( args() < 2 ) tol = 0
    tol = max((tol, epsilon(rows(X))))

    qrdp(cross(X, X), ., R = ., p = .)
    D = abs(diagonal(R))
    return(rowshape((D[order(p', 1)] :/ max(1 \ D)) :>= tol, 1))
}
end
