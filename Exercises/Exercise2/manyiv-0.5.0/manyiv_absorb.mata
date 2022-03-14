cap mata mata drop ManyIVreg_Absorb()
cap mata mata drop ManyIVreg_Factor()
cap mata mata drop ManyIVreg_Absorb_New()
cap mata mata drop sf_helper_counts()
cap mata mata drop sf_helper_corder()
cap mata mata drop sf_helper_unique()

mata
class ManyIVreg_Factor
{
    real scalar nobs
    real scalar encoded
    real scalar nlevels
    real matrix info
    real colvector index
    real colvector nj
    real colvector skip
    real colvector groupid
}

class ManyIVreg_Absorb
{
    transmorphic  absorbinfo
    string vector absorbvars
    string scalar touse

    real scalar df
    real scalar nobs
    real scalar nabsorb
    real scalar hdfetol
    real scalar maxiter
    real scalar squarem
    real vector nlevels
    real vector total_nlevels
    real vector encoded
    real vector original
    real scalar ncombined
    real scalar nlevels_combined
    real scalar nsingledrop
    real scalar nsingletons
    real scalar d_computed
    real vector d_order
    real colvector d_projection

    void new()
    void init()

    real scalar allhaveskip()
    real scalar onesingleton()
    real scalar makepanel()
    real matrix info()
    real colvector d_projection()
    real colvector index()
    real colvector nj()
    real colvector skip()
    real colvector groupid()
    real colvector groupindex()
    real colvector singletonindex()

    void set_info()
    void set_index()
    void set_nj()
    void set_skip()
    void set_groupid()

    void encode()
    void combine()
    void append()
    void dropsingletons()
    void dropfromindex()
    void skipsingletons()
    void flagredundant()
    void exportc()
    void importc()

    real matrix hdfe()
    void _hdfe()
    void _hdfe_fpi()
    void _hdfe_squarem()
    real scalar _hdfe_halperin()
    real scalar _hdfe_halperin_symm()
    real scalar _demean()
}

class ManyIVreg_Absorb scalar function ManyIVreg_Absorb_New(string vector _absorbvars, string scalar _touse, | real scalar _squarem)
{
    class ManyIVreg_Absorb scalar Absorb
    if ( args() < 3 ) _squarem = 0
    Absorb = ManyIVreg_Absorb()
    Absorb.init(_absorbvars, _touse, _squarem)
    return(Absorb)
}

void function ManyIVreg_Absorb::new()
{
    this.hdfetol   = 1e-8
    this.maxiter   = 1000
    this.encoded   = 0
    this.squarem   = 0
    this.ncombined = 0
    this.nobs      = 0
    this.nlevels_combined = 0
    this.nsingledrop = 0
    this.nsingletons = 0
    this.d_computed  = 0
}

void function ManyIVreg_Absorb::init(string vector _absorbvars, string scalar _touse, | real scalar _squarem)
{
    real scalar j
    absorbvars    = _absorbvars
    touse         = _touse
    nabsorb       = length(absorbvars)
    absorbinfo    = asarray_create()
    total_nlevels = 0
    nlevels       = J(1, nabsorb, .)
    encoded       = J(1, nabsorb, 0)
    original      = J(1, nabsorb, 1)
    for(j = 1; j <= nabsorb; j++) {
        nlevels[j] = makepanel(absorbvars[j])
        total_nlevels = total_nlevels + nlevels[j]
    }
    if ( args() >= 3 ) squarem = _squarem
    // Note: The idea is for df to not be modified by append
    df = total_nlevels - nabsorb
}

real scalar function ManyIVreg_Absorb::makepanel(string scalar var)
{
    real scalar nl
    string colvector svar
    real colvector varindex, rvar
    real matrix varinfo

    if ( strpos(st_vartype(var), "str") ) {
        svar     = st_sdata(., var, touse)
        varindex = order(svar, 1)
        varinfo  = panelsetup(svar[varindex], 1)
    }
    else {
        rvar     = st_data(., var, touse)
        varindex = order(rvar, 1)
        varinfo  = panelsetup(rvar[varindex], 1)
    }

    nl   = rows(varinfo)
    nobs = length(varindex)

    asarray(absorbinfo, var + ".index", varindex)
    asarray(absorbinfo, var + ".info",  varinfo)
    asarray(absorbinfo, var + ".nj",    varinfo[., 2] :- varinfo[., 1] :+ 1)
    asarray(absorbinfo, var + ".skip",  J(nl, 1, 0))

    return(nl)
}

void function ManyIVreg_Absorb::combine()
{
    real colvector sel, gid, order, gencode, gcumsum
    real colvector combinedid, combinedix, combinednj, combinedskip
    real rowvector counts
    real scalar i, j, goffset, ioffset, refj, start
    real matrix ginfo, combinedinfo

    if ( nabsorb < 2 ) return

    encode()
    if ( ncombined ) {
        refj  = 0
    }
    else {
        refj = 1
        nlevels_combined = nlevels[1]
        ncombined = 1
    }
    start = ncombined + 1

    for(j = start; j <= nabsorb; j++) {
        ioffset = 0
        goffset = 0

        combinedid   = J(nobs, 1, .)
        combinedix   = index(refj)
        combinedinfo = J(0, 2, 0)
        combinednj   = J(0, 1, 0)
        combinedskip = J(0, 1, 0)

        for(i = 1; i <= nlevels_combined; i++) {
            sel = panelsubmatrix(index(refj), i, info(refj))
            gid = groupid(j)[sel]

            counts  = colshape(sf_helper_counts(gid, nlevels[j]), 1)
            order   = sf_helper_corder(gid, gcumsum = runningsum(counts))
            gcumsum = select(gcumsum, counts)
            gencode = runningsum(counts :> 0)
            if ( length(gcumsum) > 1 ) {
                ginfo = (1 \ (gcumsum[|1 \ (length(gcumsum) - 1)|] :+ 1)), gcumsum
            }
            else {
                ginfo = 1, gcumsum
            }

            combinedix[|info(refj)[i, 1] \ info(refj)[i, 2]|] = sel[order]
            combinedid[sel] = gencode[gid] :+ goffset
            combinednj      = combinednj \ select(counts, counts)
            combinedinfo    = combinedinfo \ (ginfo :+ ioffset)
            combinedskip    = combinedskip \ (skip(refj)[i]? J(rows(gcumsum), 1, 1): select(skip(j), counts))

            goffset = goffset + gencode[length(gencode)]
            ioffset = ioffset + length(sel)
        }

        set_info(0, combinedinfo)
        set_index(0, combinedix)
        set_nj(0, combinednj)
        set_id(0, combinedid)
        set_skip(0, combinedskip)

        nlevels_combined = rows(combinedinfo)
        refj = 0
        ncombined = j
    }
}

void function ManyIVreg_Absorb::flagredundant(| real scalar skipredundant)
{
    string scalar strgib
    real matrix D12, R
    real scalar i, j, refj, prod, offset, gib, nskip
    real colvector D, sel, gid, counts, selskip, redundant, redundantj, skipj

    if ( args() < 1 ) skipredundant = 0

    // real colvector skipj, nskip, skipix, newskip
    // real colvector sel, gid, counts, gsel
    if ( nabsorb < 2 ) {
        return
    }
    else if ( d_computed ) {
        // if d_computed then this was already done in importc
    }
    else if ( nabsorb > 2 ) {
        errprintf("SE incorrect; use the C++ plugin to compute the dof adjustment with more than 2 levels\n")
        error(1234)
    }
    else {
        encode()
        if ( all(original) ) {
            refj = selectindex(max(nlevels) :== nlevels)[1]
        }
        else {
            refj = selectindex(max(nlevels :* (1 :- original)) :== nlevels)[1]
        }

        prod = 1
        for (j = 1; j <= nabsorb; j++) {
            if ( j == refj ) continue
            prod = prod * nlevels[j]
        }

        gib = ((prod * nlevels[refj] + prod:^2) * 8) / (1024:^3)
        if ( gib > 1 ) {
            strgib = strtrim(sprintf("%15.1fc", gib))
            errprintf("resource warning: %sGiB required for redundant level detection\n", strgib)
        }

        offset = 0
        D12 = J(prod, nlevels[refj], 0)
        for (j = 1; j <= nabsorb; j++) {
            if ( j == refj ) continue
            for (i = 1; i <= nlevels[refj]; i++) {
                if ( skip(refj)[i] ) continue
                sel    = panelsubmatrix(index(refj), i, info(refj))
                gid    = groupid(j)[sel]
                counts = colshape(sf_helper_counts(gid, nlevels[j]), 1)
                D12[|(offset + 1, i) \ (offset + nlevels[j], i)|] = counts
            }
            offset = offset + nlevels[j]
        }

        D = 1 :/ nj(refj)
        if ( any(skip(refj)) ) {
            selskip = selectindex(skip(refj))
            D[selskip] = J(rows(selskip), 1, 0)
        }

        // ------------------------------------------
        // This portion doesn't work with nabsorb > 2
        // ------------------------------------------
        D12 = - D12 * (D12' :* D)
        for (i = 1; i <= nabsorb; i++) {
            if ( i == refj ) continue
            j = i
        }
        skipj = skip(j)
        D = nj(j)
        if ( any(skipj) ) {
            selskip = selectindex(skipj)
            D[selskip] = J(rows(selskip), 1, 0)
        }
        _diag(D12, diagonal(D12) + D)
        // ------------------------------------------
        // This portion doesn't work with nabsorb > 2
        // ------------------------------------------

        qrd(D12, ., R = .)
        D = abs(diagonal(R))

        redundant = (D :/ max(D)) :< hdfetol
        if ( any(redundant) ) {
            if ( any(original) ) {
                df = sum(select(nlevels, original))
            }

            offset = 0
            for (j = 1; j <= nabsorb; j++) {
                skipj = skip(j)
                if ( j != refj ) {
                    redundantj = selectindex(redundant[|(offset + 1) \ (offset + nlevels[j])|])
                    if ( rows(redundantj) ) {
                        skipj[selectindex(redundantj)] = J(rows(redundantj), 1, 1)
                    }
                    offset = offset + nlevels[j]
                }
                if ( original[j] ) {
                    nskip = sum(skipj)
                    df = df - (nskip? nskip: 1)
                }
                if ( skipredundant ) set_skip(j, skipj)
            }
        }
    }
}

real colvector function ManyIVreg_Absorb::singletonindex(real scalar j)
{
    real colvector sel
    sel = selectindex(nj(j) :== 1)
    return(rows(sel)? index(j)[info(j)[sel, 1]]: sel)
}

void function ManyIVreg_Absorb::dropsingletons(|real colvector singleix, real scalar method)
{
    real scalar j
    if ( args() < 2 ) method = 1

    if ( nabsorb == 0 ) {
        singleix = J(0, 1, .)
        return
    }

    if ( method == 1 ) {
        singleix = J(0, 1, .)
        for (j = 1; j <= nabsorb; j++) {
            singleix = singleix \ singletonindex(j)
        }
        singleix = uniqrows(singleix)
        dropfromindex(singleix)

        if ( nobs :== 0 ) {
            errprintf("all singleton groups")
            error(1234)
        }
    }
    else if ( method == 2 ) {
        singleix = J(0, 1, .)
        skipsingletons()
    }
    else if ( method == 3 ) {
        singleix = J(0, 1, .)
        for (j = 1; j <= nabsorb; j++) {
            singleix = J(0, 1, .) \ singletonindex(j)
        }
        singleix = uniqrows(singleix)
    }

    nsingletons = rows(singleix)
    nsingledrop = method == 3? 0: nsingletons
}

void function ManyIVreg_Absorb::skipsingletons()
{
    real colvector selskip
    real scalar j, nskip

    if ( any(original) ) {
        df = sum(select(nlevels, original))
    }

    for (j = (ncombined? 0: 1); j <= nabsorb; j++) {
        selskip = (nj(j) :== 1)
        if ( j ) {
            if ( original[j] ) {
                nskip = sum(selskip)
                df = df - (nskip? nskip: 1)
            }
        }
        set_skip(j, selskip)
    }
}

void function ManyIVreg_Absorb::dropfromindex(real colvector dropindex, | real scalar reference)
{
    real colvector counts, sel, selix, order, singleoffset, index, nj, groupid, skip
    real matrix info
    real scalar j, nl

    if ( rows(dropindex) :== 0 ) return
    if ( args() == 1 ) reference = .

    encode()
    order = J(nobs, 1, 0)
    for (j = (ncombined? 0: 1); j <= nabsorb; j++) {
        nl      = j? nlevels[j]: nlevels_combined
        info    = info(j)
        groupid = groupid(j)
        index   = index(j)
        skip    = skip(j)
        counts  = colshape(sf_helper_counts(groupid[dropindex], nl), 1)
        sel     = selectindex(counts)
        if ( j == reference ) {
            selix  = info[sel, 1]
        }
        else {
            // note: index has a 1 to nobs index so its sort order is unique
            // order = order(index, 1)
            order[index] = 1::nobs
            selix = order[dropindex]
        }

        singleoffset = J(nobs, 1, 0)
        singleoffset[dropindex] = J(rows(dropindex), 1, 1)
        singleoffset = runningsum(singleoffset)

        index = index :- singleoffset[index]
        index[selix] = J(rows(dropindex), 1, 0)
        index = select(index, index)

        groupid[dropindex] = J(rows(dropindex), 1, 0)
        groupid = select(groupid, groupid)

        singleoffset = J(rows(info), 1, 0)
        singleoffset[sel] = select(counts, counts)
        singleoffset = runningsum(singleoffset)

        info[., 2] = info[., 2] :- singleoffset
        nj = (0 \ info[|(1, 2) \ (nl - 1, 2)|])
        info[., 1] = nj :+ 1
        nj = info[., 2] :- nj
        if ( j != reference ) {
            sel = selectindex(nj :== 0)
            if ( rows(sel) ) {
                singleoffset = J(rows(info), 1, 0)
                singleoffset[sel] = J(rows(sel), 1, 1)
                singleoffset = runningsum(singleoffset)
            }
        }
        info = select(info, nj)
        skip = select(skip, nj)
        nj   = select(nj, nj)

        if ( rows(info) < nl ) {
            groupid = groupid :- singleoffset[groupid]
        }

        set_info(j, info)
        set_index(j, index)
        set_nj(j, nj)
        set_groupid(j, groupid)
        set_skip(j, skip)

        if ( j ) {
            nlevels[j] = rows(info)
        }
        else {
            nlevels_combined = rows(info)
        }
    }

    total_nlevels = sum(nlevels)
    nobs = max((0, nobs - rows(dropindex)))
}

// Counting sort order given counts of number of repeated elements (assumes id is encoded)
real colvector sf_helper_corder(real vector id, real colvector cumsum)
{
    real colvector order, ix
    real scalar i, level, nobs
    nobs = length(id)
    ix = 0 \ cumsum[|1 \ length(cumsum) - 1|]
    order = J(nobs, 1, .)
    for (i = 1; i <= nobs; i++) {
        level = id[i]
        order[ix[level] = ix[level] + 1] = i
    }
    return(order)
}

// Count number of times each level repeats given number of levels (assumes id is encoded)
real rowvector sf_helper_counts(real vector id, real scalar nlevels)
{
    real rowvector counts
    real scalar i, j
    counts = J(1, nlevels, 0)
    i = length(id) + 1
    while (--i) {
        j = id[i]
        counts[j] = counts[j] + 1
    }
    return(counts)
}

// Unique elements given number of levels (assumes id is encoded)
real vector sf_helper_unique(real vector id, real scalar nlevels)
{
    real rowvector present
    real scalar i
    if ( rows(id) > cols(id) ) {
        present = J(nlevels, 1, 0)
    }
    else {
        present = J(1, nlevels, 0)
    }
    i = length(id) + 1
    while (--i) {
        present[id[i]] = 1
    }
    return(selectindex(present))
}

void function ManyIVreg_Absorb::encode()
{
    real colvector sel, groupid
    real scalar i, j

    for(j = 1; j <= nabsorb; j++) {
        if ( encoded[j] ) continue
        groupid = J(nobs, 1, .)
        for(i = 1; i <= nlevels[j]; i++) {
            sel = panelsubmatrix(index(j), i, info(j))
            groupid[sel] = J(length(sel), 1, i)
        }
        set_groupid(j, groupid)
        encoded[j] = 1
    }
}

void function ManyIVreg_Absorb::append(class ManyIVreg_Absorb scalar Absorb)
{
    real vector exists, append_index
    real scalar j

    if ( (nabsorb == 0) | (Absorb.nabsorb == 0) ) {
        return
    }

    if ( nobs != Absorb.nobs ) {
        errprintf("cannot append different number of observations (%g vs %g)\n", nobs, Absorb.nobs)
        error(1234)
    }

    exists = J(Absorb.nabsorb, 1, 0)
    for (j = 1; j <= Absorb.nabsorb; j++) {
        if ( any(Absorb.absorbvars[j] :== absorbvars) ) {
            exists[j] = 1
        }
    }

    if ( all(exists) ) {
        printf("all absorb variables already exist; no append\n")
        return
    }
    append_index = selectindex(!exists)

    absorbvars = absorbvars, Absorb.absorbvars[append_index]
    nabsorb = nabsorb + Absorb.nabsorb
    hdfetol = min((hdfetol, Absorb.hdfetol))
    maxiter = max((maxiter, Absorb.maxiter))
    nlevels = nlevels, Absorb.nlevels[append_index]
    total_nlevels = total_nlevels + Absorb.total_nlevels
    encoded  = encoded, Absorb.encoded[append_index]
    original = original, J(1, rows(append_index), 0)

    for (j = 1; j <= Absorb.nabsorb; j++) {
        if ( exists[j] ) continue
        asarray(absorbinfo, Absorb.absorbvars[j] + ".index", Absorb.index(j))
        asarray(absorbinfo, Absorb.absorbvars[j] + ".info",  Absorb.info(j))
        asarray(absorbinfo, Absorb.absorbvars[j] + ".nj",    Absorb.nj(j))
        asarray(absorbinfo, Absorb.absorbvars[j] + ".skip",  Absorb.skip(j))
    }
}

real colvector function ManyIVreg_Absorb::d_projection(| real scalar base)
{
    real scalar i, j, jix, nskip
    real matrix M, S, invS, invSM, invSMM
    real colvector D, sel, gid, selskip
    real rowvector s_i

    if ( args() < 1 ) base = 0

    encode()
    if ( d_computed ) {
        return(d_projection)
    }
    else if ( nabsorb == 0 ) {
        return(0)
    }
    else if ( nabsorb == 1 ) {
        D = (1 :/ nj(1))[groupid(1)]
        if ( base ) {
            D[groupindex(1, base)] = J(nj(1)[base], 1, 0)
        }
        if ( any(skip(1)) ) {
            for(i = 1; i <= nlevels[1]; i++) {
                if ( skip(1)[i] ) D[groupindex(1, i)] = J(nj(1)[i], 1, 0)
            }
        }
        return(D)
    }
    else if ( nabsorb == 2 ) {
        // NOTE: Base is required for at least one group bc of collinearity.
        // However, this assumes skipredundant() has been run, and that should
        // have already flagged collinear columns to be Skipped.

        nskip = 0
        for(j = 1; j <= nabsorb; j++) {
            nskip = nskip + sum(skip(j))
        }
        if ( (args() < 1) & (nskip == 0) ) base = 1

        if ( nlevels[1] > nlevels[2] ) {
            jix = 2
            j   = 1
        }
        else {
            jix = 1
            j   = 2
        }

        D = (1 :/ nj(j))[groupid(j)]
        if ( base ) {
            D[groupindex(j, base)] = J(nj(j)[base], 1, 0)
        }
        if ( any(skip(j)) ) {
            for(i = 1; i <= nlevels[j]; i++) {
                if ( skip(j)[i] ) {
                    D[groupindex(j, i)] = J(nj(j)[i], 1, 0)
                }
            }
        }

        selskip = selectindex(skip(jix))
        nskip   = rows(selskip)

        M = J(nlevels[j], nlevels[jix], 0)
        S = J(nlevels[jix], nlevels[jix], 0)
        for(i = 1; i <= nlevels[j]; i++) {
            sel = panelsubmatrix(index(j), i, info(j))
            gid = groupid(jix)[sel]
            s_i = sf_helper_counts(gid, nlevels[jix])
            if ( nskip ) {
                s_i[selskip] = J(1, nskip, 0)
            }
            if ( (i == base) | skip(j)[i] ) {
                S = S + diag(s_i)
            }
            else {
                M[i, .] = s_i / length(sel)
                S  = S + diag(s_i) - (s_i' * s_i) / length(sel)
            }
        }

        if ( nskip ) {
            selskip = selectindex(!skip(jix))
            invS    = J(rows(S), cols(S), 0)
            invS[selskip, selskip] = invsym(S[selskip, selskip])
        }
        else {
            invS = invsym(S)
        }
        invSM  = M * invS
        invSMM = invSM :* M

        for(i = 1; i <= nlevels[jix]; i++) {
            if ( skip(jix)[i] ) continue
            sel = panelsubmatrix(index(jix), i, info(jix))
            D[sel] = D[sel] + invS[groupid(jix)[sel], i] - invSM[groupid(j)[sel], i]
            D = D + invSMM[groupid(j), i] - invS[groupid(jix), i] :* M[groupid(j), i]
        }

        return(D)
    }
    else {
        errprintf("use the C++ plugin to compute the projection with more than 2 levels\n")
        error(1234)
    }
}

real matrix function ManyIVreg_Absorb::info(real scalar j)
{
    return(asarray(absorbinfo, (j? absorbvars[j]: "[combined]") + ".info"))
}

real colvector function ManyIVreg_Absorb::index(real scalar j)
{
    return(asarray(absorbinfo, (j? absorbvars[j]: "[combined]") + ".index"))
}

real colvector function ManyIVreg_Absorb::nj(real scalar j)
{
    return(asarray(absorbinfo, (j? absorbvars[j]: "[combined]") + ".nj"))
}

real colvector function ManyIVreg_Absorb::skip(real scalar j)
{
    return(asarray(absorbinfo, (j? absorbvars[j]: "[combined]") + ".skip"))
}

real colvector function ManyIVreg_Absorb::groupid(real scalar j)
{
    return(asarray(absorbinfo, (j? absorbvars[j]: "[combined]") + ".groupid"))
}

real colvector function ManyIVreg_Absorb::groupindex(real scalar j, real scalar i)
{
    return(panelsubmatrix(index(j), i, info(j)))
}

void function ManyIVreg_Absorb::set_info(real scalar j, real matrix info)
{
    return(asarray(absorbinfo, (j? absorbvars[j]: "[combined]") + ".info", info))
}

void function ManyIVreg_Absorb::set_index(real scalar j, real colvector index)
{
    return(asarray(absorbinfo, (j? absorbvars[j]: "[combined]") + ".index", index))
}

void function ManyIVreg_Absorb::set_nj(real scalar j, real colvector nj)
{
    return(asarray(absorbinfo, (j? absorbvars[j]: "[combined]") + ".nj", nj))
}

void function ManyIVreg_Absorb::set_skip(real scalar j, real colvector skip)
{
    return(asarray(absorbinfo, (j? absorbvars[j]: "[combined]") + ".skip", skip))
}

void function ManyIVreg_Absorb::set_groupid(real scalar j, real colvector groupid)
{
    return(asarray(absorbinfo, (j? absorbvars[j]: "[combined]") + ".groupid", groupid))
}

real scalar function ManyIVreg_Absorb::allhaveskip()
{
    real scalar j, alls
    alls = 1
    for(j = 1; j <= nabsorb; j++) {
        // alls = alls & (sum(skip(j)) > 1)
        alls = alls & any(skip(j))
    }
    return(alls)
}

real scalar function ManyIVreg_Absorb::onesingleton()
{
    real scalar j, one
    one = 0
    for(j = 1; j <= nabsorb; j++) {
        one = one | (sum(nj(j) :== 1) :== 1)
    }
    return(one)
}

real matrix function ManyIVreg_Absorb::hdfe(real matrix X, | real scalar base)
{
    _hdfe(X, base)
    return(X)
}

void function ManyIVreg_Absorb::_hdfe(real matrix X, | real scalar base)
{
    if ( nabsorb == 1 ) {
        (void) _demean(X, 1, base)
    }
    else if ( nabsorb > 1 ) {
        // Note: base option doesn't work with multipe FE
        if ( squarem ) {
            _hdfe_squarem(X)
        }
        else {
            _hdfe_fpi(X)
        }
    }
}

// Note: In the past you tried breaking the iteration if the tolerance was
//       achieved by any given _pair_ if fixed effects. This doesn't work;
//       the easiest counter example is if fe2 is collinear with fe1 but fe3
//       is not. Then the algorithm would incorrectly stop at fe2.

// Note: This is basically fixed point iteration
void function ManyIVreg_Absorb::_hdfe_fpi(real matrix X)
{
    real scalar i, dev
    i = 0
    while ( i++ < maxiter ) {
        dev = _hdfe_halperin(X)
        if ( dev < hdfetol ) break
    }
    if ( i > maxiter ) {
        errprintf("maximum number of hdfe iterations exceeded (%g)\n", maxiter)
        error(1234)
    }
    else {
        printf("hdfe convergence (fpi) after %g projections (error = %5.3g)\n", i * nabsorb, dev)
    }
}

void function ManyIVreg_Absorb::_hdfe_squarem(real matrix X)
{
    real matrix X1, X2, Q1, Q2
    real scalar i, feval, dev, mstep
    real rowvector alpha, stepmin, stepmax, sr2, sv2, stepsel

    i = 0
    feval = 0
    stepmax = J(1, cols(X), 1)
    stepmin = J(1, cols(X), 1)
    mstep = 4
    while ( i++ < maxiter ) {
        feval++
        dev = _hdfe_halperin(X1 = X)
        if ( dev < hdfetol ) {
            X = X1
            break
        }

        feval++
        Q1  = X1 - X
        dev = _hdfe_halperin(X2 = X1)
        if ( dev < hdfetol ) {
            X = X2
            break
        }

        // Form quadratic terms
        Q2  = X2 - X1
        sr2 = colsum(Q1:^2)
        sv2 = colsum((Q2 - Q1):^2)

        // Get the step-size
        feval++
        alpha = colmax((stepmin \ colmin((stepmax \ sqrt(sr2 :/ sv2)))))
        X     = X + 2 :* alpha :* Q1 :+ alpha:^2 :* (Q2 - Q1)
        dev   = _hdfe_halperin(X)
        if ( dev < hdfetol ) break

        stepsel = alpha :== stepmax
        if ( any(stepsel) ) {
            stepsel = selectindex(stepsel)
            stepmax[stepsel] = mstep * stepmax[stepsel]
        }
    }

    if ( i > maxiter ) {
        errprintf("maximum number of hdfe iterations exceeded (%g)\n", maxiter)
        error(1234)
    }
    else {
        printf("hdfe convergence (squarem) after %g projections (error = %5.3g)\n", feval * nabsorb, dev)
    }
}

real scalar function ManyIVreg_Absorb::_hdfe_halperin(real matrix X)
{
    real scalar j, dev
    dev = 1
    for (j = 1; j <= nabsorb; j++) {
        dev = _demean(X, j)
    }
    return(dev)
}

real scalar function ManyIVreg_Absorb::_hdfe_halperin_symm(real matrix X)
{
    real scalar j, dev
    dev = 1
    for (j = 1; j <= nabsorb; j++) {
        dev = _demean(X, j)
    }
    for (j = nabsorb - 1; j > 1; j--) {
        dev = _demean(X, j)
    }
    return(dev)
}

real scalar function ManyIVreg_Absorb::_demean(real matrix X, real scalar j, | real scalar base)
{
    real colvector avg
    real colvector sel
    real matrix X_i
    real scalar i, dev
    dev = 0
    if ( args() < 3 ) base = 0
    for(i = 1; i <= nlevels[j]; i++) {
        if ( (i == base) | skip(j)[i] ) continue
        sel = panelsubmatrix(index(j), i, info(j))
        X_i = X[sel, .]
        avg = (colsum(X_i) :/ length(sel))
        X[sel, .] = X_i :- avg
        dev = max((abs(avg), dev))
        // dev = dev + sum(avg:^2)
    }
    return(dev)
}

void function ManyIVreg_Absorb::exportc(string scalar fname)
{
    real scalar j, fh
    colvector C
    encode()
    fh = fopen(fname, "rw")
    C = bufio()
    fbufput(C, fh, "%4bu", nabsorb)
    fbufput(C, fh, "%4bu", nobs)

    d_order = order(-colshape(nlevels, 1), 1)
    fbufput(C, fh, "%4bu", nlevels[d_order])
    for(j = 1; j <= nabsorb; j++) {
        fbufput(C, fh, "%4bu", index(d_order[j]) :- 1)
        fbufput(C, fh, "%4bu", groupid(d_order[j]) :- 1)
        fbufput(C, fh, "%4bu", ((info(d_order[j])[., 1] :- 1) \ nobs))
        fbufput(C, fh, "%4bu", skip(d_order[j]))
    }
    // printf("exported %g bytes\n",
    //        4 * (2 + nabsorb + 2 * nabsorb * nobs + sum(nlevels) + nabsorb))
    fclose(fh)
}

void function ManyIVreg_Absorb::importc(string scalar fname, | real scalar skipredundant)
{
    real colvector skipj
    real scalar j, fh, nskip
    colvector C
    if ( args() < 2 ) skipredundant = 0
    if ( any(original) ) df = sum(select(nlevels, original))

    fh = fopen(fname, "r")
    C = bufio()
    d_projection = fbufget(C, fh, "%8z", nobs)'
    for(j = 1; j <= nabsorb; j++) {
        skipj = (fbufget(C, fh, "%4bu", nlevels[d_order[j]]) :== 0)'
        nskip = sum(skipj)
        if ( original[d_order[j]] ) df = df - (nskip? nskip: 1)
        if ( skipredundant ) set_skip(d_order[j], (skip(d_order[j]) :| skipj))
    }
    fclose(fh)

    d_computed = 1
}
end
