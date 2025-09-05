"""Pure-function refactor of bpm.calcBPM logic.
API: calc_bpm_windows(records, window_size=100000, minimum_snps=2, num_pops=2, outname='out')
records: iterable of lists [pop,ploidy,scaff,pos,an,dp,geno1,geno2,...]
returns: list of dicts with keys matching original output columns
"""

from collections import defaultdict


def NestedAnova(locus_list):
    locus = []
    for l in locus_list:  # Remove missing data from the site
        locus.append([x for x in l if x != "-9"])
    r = float(len(locus))  # Number of populations
    p_i = []
    n_i = []
    ploidy_list = []
    Fssg = 0.0
    Fssi = 0.0
    Fssw = 0.0
    Fsst = 0.0
    Rssg = 0.0
    Rssi = 0.0
    Rsst = 0.0
    p_ij = []
    tac = 0
    tan = 0
    ac_ij = []
    df_w = 0.0
    FS_Nij2 = 0.0
    FS_SNij2 = 0.0
    RS_SNij2 = 0.0
    FS_SNij2_over_SNij = 0.0
    for pop_site in locus:
        FSNij2temp = 0.0
        p = []
        ac_i = []
        ploidy = float(pop_site[1])
        nnn = float(len(pop_site[6:]))
        an = ploidy * nnn
        ac = sum([float(geno) for geno in pop_site[6:]])
        n_i.append(nnn)
        p_i.append(ac / an if an != 0 else 0.0)
        ploidy_list.append(ploidy)
        pop_an = 0
        for ind in pop_site[6:]:
            p_ind = float(ind) / ploidy
            p.append(p_ind)
            ac_i.append(float(ind))
            num_alleles = 0
            for c in range(0, int(ploidy) - int(ind)):
                num_alleles += 1
                tan += 1
            for c in range(0, int(ind)):
                num_alleles += 1
                tac += 1
                tan += 1
            FS_Nij2 += ploidy**2
            FSNij2temp += ploidy**2
            pop_an += ploidy
            assert num_alleles == int(ploidy)
            df_w += float(num_alleles) - 1.0
        FS_SNij2_over_SNij += FSNij2temp / pop_an if pop_an != 0 else 0.0
        FS_SNij2 += pop_an**2
        RS_SNij2 += nnn**2
        p_ij.append(p)
        ac_ij.append(ac_i)

    p_bar = float(tac) / float(tan) if tan != 0 else 0.0
    df_g = r - 1.0
    df_i = sum([x - 1.0 for x in n_i]) if n_i else 0.0
    df_t = df_g + df_i + df_w
    fn0bis = (
        (FS_SNij2_over_SNij - (FS_Nij2 / tan)) / df_g
        if (df_g != 0 and tan != 0)
        else 0.0
    )
    fn0 = (tan - FS_SNij2_over_SNij) / df_i if df_i != 0 else 0.0
    fnb0 = (tan - (FS_SNij2 / tan)) / df_g if (df_g != 0 and tan != 0) else 0.0
    rnb0 = (
        (float(sum(n_i)) - (RS_SNij2 / float(sum(n_i)))) / df_g
        if df_g != 0 and sum(n_i) != 0
        else 0.0
    )

    # Protect against degenerate cases
    if df_t != 0 and tan != 0:
        pass

    for i, pop in enumerate(ac_ij):
        for j, ind in enumerate(pop):
            for ref in range(0, int(int(ploidy_list[i]) - int(ind))):
                Fssg += (p_i[i] - p_bar) ** 2
                Fssi += (p_ij[i][j] - p_i[i]) ** 2
                Fssw += (0 - p_ij[i][j]) ** 2
                Fsst += (0 - p_bar) ** 2
            for alt in range(0, int(ind)):
                Fssg += (p_i[i] - p_bar) ** 2
                Fssi += (p_ij[i][j] - p_i[i]) ** 2
                Fssw += (1 - p_ij[i][j]) ** 2
                Fsst += (1 - p_bar) ** 2
            Rssi += (p_ij[i][j] - p_i[i]) ** 2
            Rssg += (p_i[i] - p_bar) ** 2
            Rsst += (p_ij[i][j] - p_bar) ** 2

    # Mean squares and variance components
    try:
        FMS_g = Fssg / df_g
        FMS_i = Fssi / df_i
        FMS_w = Fssw / df_w
    except Exception:
        FMS_g = FMS_i = FMS_w = 0.0
    try:
        RMS_g = Rssg / df_g
        RMS_i = Rssi / df_i
    except Exception:
        RMS_g = RMS_i = 0.0

    try:
        fs2_w = FMS_w
        fs2_i = (FMS_i - fs2_w) / fn0 if fn0 != 0 else 0.0
        fs2_g = (FMS_g - fs2_w - fn0bis * fs2_i) / fnb0 if fnb0 != 0 else 0.0
        rs2_i = RMS_i
        rs2_g = (RMS_g - rs2_i) / rnb0 if rnb0 != 0 else 0.0
    except Exception:
        fs2_w = fs2_i = fs2_g = rs2_i = rs2_g = 0.0

    rnum = rs2_g
    rden = rs2_i + rs2_g
    fnum = fs2_g
    fden = fs2_w + fs2_g + fs2_i

    poly = not (all(p == 0.0 for p in p_i) or all(p == 1.0 for p in p_i))
    return rnum, rden, fnum, fden, poly


def calcDxy(locus_list):
    locus = []
    for l in locus_list:
        locus.append([x for x in l if x != "-9"])
    p1 = p2 = 0.0
    for i, pop_site in enumerate(locus):
        ploidy = float(pop_site[1])
        nnn = float(len(pop_site[6:]))
        an = ploidy * nnn
        ac = sum([float(geno) for geno in pop_site[6:]])
        if i == 0:
            p1 = (ac / an) if an != 0 else 0.0
        if i == 1:
            p2 = (ac / an) if an != 0 else 0.0

    dxy = (p1 * (1.0 - p2)) + (p2 * (1.0 - p1))
    afd = abs(p1 - p2)
    return dxy, afd


def calc_bpm_windows(
    records, window_size=10000, minimum_snps=2, num_pops=2, outname="out"
):
    """Process iterable of records and return list of window dicts and genome metrics"""
    results = []
    snp_count = 0
    site_count = 0
    start = 0.0
    end = float(window_size)
    winexclcount = 0
    num_wind = 0
    Locus = []
    oldscaff = None
    old_pos = None
    fst = [0.0, 0.0]
    rho = [0.0, 0.0]
    Fst = [0.0, 0.0]
    Rho = [0.0, 0.0]
    dxy = 0.0
    Dxy = 0.0
    dn = 0
    Dn = 0
    afd = 0.0
    AFD = 0.0

    # records assumed sorted by scaffold and position
    for i, line in enumerate(records):
        pop, ploidy, scaff, pos, an, dp = line[:6]
        pos = float(pos)
        if i % 100000 == 0:
            pass
        if i == 0:
            old_pos = pos
            Locus = []
            oldscaff = scaff

        if pos > start and pos <= end and scaff == oldscaff:
            if pos == old_pos:
                Locus.append(line)
            elif len(Locus) == num_pops:
                rnum, rden, fnum, fden, poly = NestedAnova(Locus)
                if poly:
                    if num_pops == 2:
                        d, a = calcDxy(Locus)
                        dxy += d
                        Dxy += d
                        afd += a
                        AFD += a
                        if a == 1.0:
                            dn += 1
                            Dn += 1
                    snp_count += 1
                    site_count += 1
                    fst = [sum(x) for x in zip(fst, [fnum, fden])]
                    rho = [sum(x) for x in zip(rho, [rnum, rden])]
                    Fst = [sum(x) for x in zip(Fst, [fnum, fden])]
                    Rho = [sum(x) for x in zip(Rho, [rnum, rden])]
                else:
                    site_count += 1
                Locus = [line]
                old_pos = pos
            else:
                Locus = [line]
                old_pos = pos

        elif pos > end or scaff != oldscaff:
            if len(Locus) == num_pops:
                rnum, rden, fnum, fden, poly = NestedAnova(Locus)
                if poly:
                    snp_count += 1
                    site_count += 1
                    fst = [sum(x) for x in zip(fst, [fnum, fden])]
                    rho = [sum(x) for x in zip(rho, [rnum, rden])]
                    Fst = [sum(x) for x in zip(Fst, [fnum, fden])]
                    Rho = [sum(x) for x in zip(Rho, [rnum, rden])]
                    if num_pops == 2:
                        d, a = calcDxy(Locus)
                        dxy += d
                        Dxy += d
                        afd += a
                        AFD += a
                        if a == 1.0:
                            dn += 1
                            Dn += 1
                else:
                    site_count += 1

            if snp_count >= minimum_snps:
                num_wind += 1
                try:
                    fst_val = fst[0] / fst[1]
                    fac = rho[0] / rho[1]
                    rho_i = fac / (1 + fac)
                    if num_pops == 2:
                        dxy_val = dxy / float(site_count)
                        afd_val = afd / float(site_count)
                    else:
                        dxy_val = -9
                        dn = -9
                        afd_val = -9
                    results.append(
                        {
                            "outname": outname,
                            "scaff": scaff,
                            "start": start,
                            "end": end,
                            "win_size": window_size,
                            "num_sites": site_count,
                            "num_snps": snp_count,
                            "Rho": rho_i,
                            "Fst": fst_val,
                            "dxy": dxy_val,
                            "AFD": afd_val,
                            "FixedDiff": dn,
                        }
                    )
                except ZeroDivisionError:
                    pass
            else:
                winexclcount += 1

            snp_count = 0
            site_count = 0
            fst = [0.0, 0.0]
            rho = [0.0, 0.0]
            dxy = 0.0
            dn = 0
            afd = 0.0

            if pos > end:
                while pos > end:
                    end += window_size
                start = end - window_size
            elif scaff != oldscaff:
                oldscaff = scaff
                start = 0.0
                end = window_size

            if pos > start and pos <= end and scaff == oldscaff:
                Locus = [line]
                old_pos = pos

    # Final locus
    if len(Locus) == num_pops:
        rnum, rden, fnum, fden, poly = NestedAnova(Locus)
        if poly:
            snp_count += 1
            site_count += 1
            fst = [sum(x) for x in zip(fst, [fnum, fden])]
            rho = [sum(x) for x in zip(rho, [rnum, rden])]
            Fst = [sum(x) for x in zip(Fst, [fnum, fden])]
            Rho = [sum(x) for x in zip(Rho, [rnum, rden])]
            if num_pops == 2:
                d, a = calcDxy(Locus)
                dxy += d
                Dxy += d
                afd += a
                AFD += a
                if a == 1.0:
                    dn += 1
                    Dn += 1
        else:
            site_count += 1

    if snp_count >= minimum_snps:
        num_wind += 1
        try:
            fst_val = fst[0] / fst[1]
            fac = rho[0] / rho[1]
            rho_i = fac / (1 + fac)
            if num_pops == 2:
                dxy_val = dxy / float(site_count)
                afd_val = afd / float(site_count)
            else:
                dxy_val = -9
                dn = -9
                afd_val = -9
            results.append(
                {
                    "outname": outname,
                    "scaff": scaff,
                    "start": start,
                    "end": end,
                    "win_size": window_size,
                    "num_sites": site_count,
                    "num_snps": snp_count,
                    "Rho": rho_i,
                    "Fst": fst_val,
                    "dxy": dxy_val,
                    "AFD": afd_val,
                    "FixedDiff": dn,
                }
            )
        except ZeroDivisionError:
            pass

    # genome-wide metrics - return minimal summary too
    FAC = Rho[0] / Rho[1] if Rho[1] != 0 else 0.0
    rho_G = FAC / (1 + FAC) if (1 + FAC) != 0 else 0.0
    Fst_G = Fst[0] / Fst[1] if Fst[1] != 0 else 0.0
    # use Dxy and AFD totals tracked as Dxy and AFD; compute per-site averages if Site_count available
    if (
        num_pops == 2
        and "Site_count" in globals()
        and globals().get("Site_count", 0) != 0
    ):
        Dxy_val = Dxy / float(globals().get("Site_count"))
        AFD_val = AFD / float(globals().get("Site_count"))
    else:
        Dxy_val = -9
        AFD_val = -9

    return results, {
        "num_windows": num_wind,
        "excluded_windows": winexclcount,
        "rho": rho_G,
        "Fst": Fst_G,
        "Dxy": Dxy_val,
        "AFD": AFD_val,
    }
