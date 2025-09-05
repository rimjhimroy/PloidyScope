"""Pure-function refactor of wpm.calcwpm logic.
Expose `calc_wpm_windows(records, sampind, window_size, minimum_snps)`
"""

import numpy as np


def calc_wpm_windows(records, sampind=5, window_size=50000, minimum_snps=2):
    """
    Pure function version of WPM calculation. Input: list of records (pop, ploidy, scaff, pos, an, dp, genotypes...)
    Returns: list of dicts, one per window, with all metrics.
    """
    import numpy as np

    results = []
    cur_scaff = None
    window_start = 0.0
    window_end = window_size
    buffer = []

    # Helper to process a window
    def process_window(loci):
        if not loci:
            return None
        pop, ploidy = loci[0][0], int(loci[0][1])
        AN = sampind * ploidy
        n = float(AN)
        aw = sum(1.0 / j for j in range(1, AN))
        a2 = sum(1.0 / (j**2) for j in range(1, AN))
        bw = a2 + 1.0 / (n**2)
        b1 = (n + 1) / (3 * (n - 1))
        b2 = (2 * (n**2 + n + 3)) / (9 * (n * (n - 1)))
        c1 = b1 - (1 / aw)
        c2 = b2 - (n + 2) / (aw * n) + a2 / aw**2
        e1 = c1 / aw
        e2 = c2 / (aw**2 + a2)
        snp_count = 0
        tot_count = 0
        num_sing = 0
        p = []
        Ehet = []
        afs = [0 for _ in range(AN + 1)]
        for line in loci:
            gt = [x for x in line[6:] if x != "-9"]
            if len(gt) >= sampind:
                sgt = np.random.choice(gt, size=sampind, replace=False)
                sac = sum([int(x) for x in sgt])
                tot_count += 1
                if sac != 0 and sac != AN:
                    if sac == 1:
                        num_sing += 1
                    p1 = float(sac) / float(AN)
                    p.append(p1)
                    Ehet.append(p1 * (1 - p1))
                    afs[sac] += 1
                    snp_count += 1
        Pi = 0.0
        h = 0.0
        L = 0.0
        S = float(sum(afs[1:-1]))
        W = S / aw if aw != 0 else 0.0
        W2 = S * (S - 1) / ((aw**2) + a2) if (aw**2 + a2) != 0 else 0.0
        for j in range(1, AN):
            Pi += afs[j] * j * (AN - j)
            h += afs[j] * (j**2)
            L += afs[j] * j
        Pi = 2 * Pi / (n * (n - 1)) if n > 1 else 0.0
        h = 2 * h / (n * (n - 1)) if n > 1 else 0.0
        L = L / (n - 1) if n > 1 else 0.0
        div = Pi / float(tot_count) if tot_count > 0 else 0.0
        varPi_W = e1 * S + e2 * S * (S - 1)
        varPi_L1 = (((n - 2) / (6 * n - 6)) * W) if n > 2 else 0.0
        varPi_L2 = (
            (
                ((18 * n**2 * (3 * n + 2) * bw) - (88 * n**3 + 9 * n**2 - 13 * n + 6))
                / (9 * n * (n - 1) ** 2)
                * W2
            )
            if n > 1
            else 0.0
        )
        varL_W = (
            (((n / (2 * n - 2)) - 1 / aw) * W)
            + (
                (
                    a2 / aw**2
                    + (2 * (n / (n - 1)) ** 2) * a2
                    - 2 * (n * a2 - n + 1) / ((n - 1) * aw)
                    - (3 * n + 1) / (n - 1)
                )
                * W2
            )
            if n > 1 and aw != 0
            else 0.0
        )
        varPi_L = varPi_L1 + varPi_L2
        try:
            D = (Pi - W) / (varPi_W**0.5) if varPi_W > 0 else 0.0
            H = (Pi - L) / (varPi_L**0.5) if varPi_L > 0 else 0.0
            E = (L - W) / (varL_W**0.5) if varL_W > 0 else 0.0
        except ZeroDivisionError:
            D = H = E = 0.0
        return {
            "pop": pop,
            "ploidy": ploidy,
            "sampind": sampind,
            "scaff": loci[0][2],
            "start": window_start,
            "end": window_end,
            "window_size": window_size,
            "num_snps": snp_count,
            "num_sites": tot_count,
            "num_singletons": num_sing,
            "avg_freq": float(np.mean(p)) if p else 0.0,
            "avg_Ehet": float(np.mean(Ehet)) if Ehet else 0.0,
            "Diversity": div,
            "ThetaW": W,
            "Pi": Pi,
            "ThetaH": h,
            "ThetaL": L,
            "D": D,
            "H": H,
            "E": E,
        }

    # Main loop
    for rec in records:
        pop, ploidy, scaff, pos, an, dp = rec[:6]
        pos = int(float(pos))
        if cur_scaff is None:
            cur_scaff = scaff
        if scaff != cur_scaff:
            if buffer:
                win = process_window(buffer)
                if win:
                    results.append(win)
            buffer = []
            cur_scaff = scaff
            window_start = 0.0
            window_end = window_size
        if pos > window_end:
            if buffer:
                win = process_window(buffer)
                if win:
                    results.append(win)
            buffer = []
            while pos > window_end:
                window_end += window_size
            window_start = window_end - window_size
        buffer.append(rec)
    if buffer:
        win = process_window(buffer)
        if win:
            results.append(win)
    return results
