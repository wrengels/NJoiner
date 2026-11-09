#!/usr/bin/env python3
"""
njoiner.py -- a line-by-line Python port of the alignment core of NJoiner
(W. R. Engels, University of Wisconsin, 2010; Objective-C/Cocoa, Mac OS X 10.6).

NJoiner reconstructs non-homologous end-joining (NHEJ) events.  Given the
original sequence with the cut marked by carets (^ for a blunt cut, ^XXXX^ for
a staggered cut with overhang XXXX) and the sequence of a repair product, it
makes two Smith-Waterman-Gotoh local alignments -- product against the
sequence left of the break, and product against the sequence right of it --
and combines them into a three-row alignment (left / product / right) from
which it reads off microhomology, filler (inserted) bases, bases missing on
each side, and whether a filler stretch has a template in the original.

This port follows NJAligner.m (10 Sep 2010) and NJoin.m (7 Sep 2010)
statement by statement, including their quirks, so that it reproduces the
original program's classifications rather than an idealized version of them.
Quirks kept deliberately are marked "QUIRK".

Usage:
    python3 njoiner.py INPUT.txt [-o OUT.csv] [--diagrams OUT.txt] [--html OUT.html]
                       [--params poster2011|april2010|output2010]

INPUT.txt is an NJoiner input file: a line "Original, <seq with ^ marks>"
(or LEFT/RIGHT/OVERHANG lines), then one line per product "name, sequence".
Lines whose name contains "//" in the first five characters are comments.
"""
import sys, re, csv, argparse, html

PARAMS = {
    # The settings saved in parameters.plist and userPrefs.flist (12-14 April 2010),
    # which are those in force when "No109Rr3EJ output.csv" was written.
    'april2010': dict(match=4, mismatch=-8, ambiguity=0, delOpen=-16, delCont=-2,
                      insOpen=-12, insCont=-2, minFillerLength=6, minFillerPercentMatch=70),
    # A set that reproduces every row of "No109Rr3EJ output.csv" (12 April 2010):
    # deletion-open -14 instead of -16 (other values as saved).
    'output2010': dict(match=4, mismatch=-8, ambiguity=0, delOpen=-14, delCont=-2,
                       insOpen=-12, insCont=-2, minFillerLength=6, minFillerPercentMatch=70),
    # The default set printed on the Keystone poster (January 2011); closest to
    # the summary tabulation behind Table 3 of the Position Effects manuscript.
    'poster2011': dict(match=4, mismatch=-8, ambiguity=0, delOpen=-16, delCont=-4,
                       insOpen=-12, insCont=-4, minFillerLength=6, minFillerPercentMatch=70),
}

SPACE, HYPHEN, MATCH, AMB = ' ', '-', '|', ':'


def swg(a, b, P):
    """Smith-Waterman-Gotoh as in -doSmithWatermanGotohFor:... (NJAligner.m).
    a, b are strings (0-based here; the original used 1-based arrays).
    Returns (vmax, ai, bi): the lists of 1-based indices along the optimal
    local path, in increasing order, one entry per alignment column."""
    na, nb = len(a), len(b)
    R = nb + 1
    N = (na + 1) * R
    v = [0] * N; e = [0] * N; f = [0] * N
    bt = [0] * N; btE = [0] * N; btF = [0] * N
    A = ' ' + a; B = ' ' + b
    M, MM, AMBg = P['match'], P['mismatch'], P['ambiguity']
    dO, dC, iO, iC = P['delOpen'], P['delCont'], P['insOpen'], P['insCont']
    for i in range(1, na + 1):
        for j in range(1, nb + 1):
            ij = i * R + j
            ev = v[ij - 1] + dO; ee = e[ij - 1] + dC
            if ev > ee: e[ij] = ev; btE[ij] = j - 1
            else: e[ij] = ee; btE[ij] = btE[ij - 1]
            fv = v[ij - R] + iO; ff = f[ij - R] + iC
            if fv > ff: f[ij] = fv; btF[ij] = i - 1
            else: f[ij] = ff; btF[ij] = btF[ij - R]
            mv = M if A[i] == B[j] else MM
            # QUIRK: the original tests (a[i]=='N' || b[i]=='N') -- b indexed by i,
            # not j.  b is the reference, which never contains N, so in effect only
            # an N in the product sequence counts as ambiguous.  (Beyond the end of
            # b the original read unallocated memory; treated here as "not N".)
            if A[i] == 'N' or (i <= nb and B[i] == 'N'):
                mv = AMBg
            vt = v[ij - R - 1] + mv; ft = f[ij]; et = e[ij]
            if vt >= max(et, ft): vb = vt; bt[ij] = 0
            elif ft >= max(vt, et): vb = ft; bt[ij] = -1
            else: vb = et; bt[ij] = 1
            v[ij] = vb if vb > 0 else 0
    imax = jmax = 1; vmax = 0
    for i in range(1, na + 1):
        for j in range(1, nb + 1):
            if v[i * R + j] > vmax: vmax = v[i * R + j]; imax, jmax = i, j
    # backtrace (QUIRK: 'pointer' is uninitialized in the original; the
    # starting cell of a local maximum always has bt==0, so 0 is used)
    i, j = imax, jmax
    path = []
    pointer = 0; iT = jT = 10000
    while v[i * R + j] > 0:
        path.append((i, j))
        newp = bt[i * R + j]
        if pointer == 0 and newp == -1: iT = btF[i * R + j]
        if pointer == 0 and newp == 1: jT = btE[i * R + j]
        if i <= iT and j <= jT: pointer = newp
        if i > iT: pointer = -1
        if j > jT: pointer = 1
        if pointer < 1: i -= 1
        if pointer > -1: j -= 1
    path.reverse()
    return vmax, [p[0] for p in path], [p[1] for p in path]


def complement(c):
    return {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}.get(c, 'N')


def matcher(a, b):
    if a == SPACE or b == SPACE: return SPACE
    if a == HYPHEN or b == HYPHEN: return SPACE
    if a == 'N' or b == 'N': return AMB
    if a == b: return MATCH
    return SPACE


class Trip:
    __slots__ = ('top', 'mid', 'bot')
    def __init__(self, t, m, b): self.top, self.mid, self.bot = t, m, b
    def merge(self, o):
        if self.top == 0: self.top = o.top
        if self.bot <= 0: self.bot = o.bot


def gaps_to_minus1(x):
    # "Convert gaps to -1": a repeated index means the other sequence advanced alone
    x = list(x)
    for k in range(len(x) - 1, 0, -1):
        if x[k] == x[k - 1]: x[k] = -1
    return x


def align(joinS, leftS, rightS, original, P):
    """Port of -[NJAligner main] and the bookkeeping in -[NJoin doAlign].
    Returns a dict with the CSV fields and the display rows."""
    join, left, right = joinS.upper(), leftS.upper(), rightS.upper()
    nleft = len(left)
    _, jl, lj = swg(join, left, P)
    _, jr, rj = swg(join, right, P)
    jl, lj = gaps_to_minus1(jl), gaps_to_minus1(lj)
    jr, rj = gaps_to_minus1(jr), gaps_to_minus1(rj)
    # 1-based views, padded with a trailing 0 as the C arrays were
    JL = [None] + jl + [0]; LJ = [None] + lj + [0]
    JR = [None] + jr + [0]; RJ = [None] + rj + [0]
    jllength = len(jl) + 1; jrlength = len(jr) + 1

    trip = []
    newt = None
    for i in range(1, jllength):
        newt = Trip(LJ[i], JL[i], 0); trip.append(newt)
    lastJ = newt.mid
    fill_loc, fill_len = min(newt.top, newt.mid), 0      # QUIRK: location in mixed coordinates
    while lastJ < JR[1] - 1:                               # filler between the two alignments
        lastJ += 1; trip.append(Trip(0, lastJ, 0)); fill_len += 1
    j = 1
    micro_loc, micro_len = 0, 0
    if lastJ >= JR[1]:                                     # the two alignments overlap: microhomology
        i = 0
        k = trip[i].mid
        while JR[1] < k:
            k -= 1; trip.insert(0, Trip(0, k, 0))
        while JR[1] != trip[i].mid and i < len(trip) - 1: i += 1
        micro_loc = i
        if i < jllength:
            while JR[j] != lastJ + 1 and j < jrlength:
                if i < len(trip): oldt = trip[i]
                else:
                    oldt = Trip(0, -1, 0); trip.insert(i, oldt)
                nt = Trip(0, JR[j], RJ[j])
                if oldt.mid < 0 and nt.mid > 0:
                    oldt.bot = -1; i += 1
                elif oldt.mid > 0 and nt.mid < 0:
                    nt.top = -1; trip.insert(i, nt); jllength += 1; i += 1; j += 1
                elif oldt.mid == nt.mid:
                    oldt.merge(nt); j += 1; i += 1
                else:
                    i += 1
                micro_len += 1
        i -= 1
        while trip[i].mid < 0:                             # trim microhomology ending in a gap
            micro_len -= 1; i -= 1
    while j < jrlength:                                    # rest of the right alignment
        trip.append(Trip(0, JR[j], RJ[j])); j += 1

    # missing left bases
    jj = jllength - 2
    lastJ = trip[jj].top
    ml_loc, ml_len = jllength - 1, 0
    while lastJ < nleft and jj < len(trip) - 1:
        jj += 1; lastJ += 1; trip[jj].top = lastJ; ml_len += 1
    # missing right bases
    jj = 0
    while not trip[jj].bot: jj += 1
    lastJ = trip[jj].bot
    mr_loc, mr_len = jj, 0
    while lastJ > 1 and jj > 0:
        jj -= 1; lastJ -= 1; trip[jj].bot = lastJ; mr_loc -= 1; mr_len += 1

    def ch(s, x):
        return SPACE if x == 0 else HYPHEN if x == -1 else s[x - 1]
    t = ''.join(ch(leftS, q.top) for q in trip)
    m = ''.join(ch(joinS, q.mid) for q in trip)
    b = ''.join(ch(rightS, q.bot) for q in trip)
    T, Mm, Bb = t.upper(), m.upper(), b.upper()
    mt = [matcher(T[i], Mm[i]) for i in range(len(T))]
    mb = [matcher(Mm[i], Bb[i]) for i in range(len(T))]
    for i in range(ml_loc, ml_loc + ml_len):
        if i < len(mt): mt[i] = SPACE
    for i in range(mr_loc, mr_loc + mr_len):
        if i < len(mb): mb[i] = SPACE
    perfect = True
    for i in range(micro_loc, micro_loc + micro_len):
        if i >= len(mt) or mt[i] != MATCH or mb[i] != MATCH: perfect = False
    if micro_len == 0: perfect = True
    imperfect = not perfect
    complexA = False
    rs = ls = rd = ld = False
    for i in range(len(T)):
        if HYPHEN in (t[i], m[i], b[i]): complexA = True
        if ld and mt[i] == MATCH: complexA = True
        if rd and mb[i] == MATCH: complexA = True
        if mt[i] != SPACE: ls = True
        if mb[i] != SPACE: rs = True
        if ls and mt[i] == SPACE: ld = True
        if rs and mb[i] == SPACE: rd = True

    # filler template search
    templ = 0
    fill_loc2, fill_len2 = fill_loc, fill_len
    if fill_len >= P['minFillerLength']:
        orig = original.upper().replace('^', '')
        lo = max(0, fill_loc); hi = min(len(Mm), fill_loc + fill_len)
        fD = Mm[lo:hi] if hi > lo else ''
        fG = ''.join(complement(c) for c in reversed(fD))
        if fD:
            fitD = swg(orig, fD, P)[0]; fitG = swg(orig, fG, P)[0]
            score = max(fitD, fitG) * 100 / (fill_len * P['match'])
            templ = 0 if score < P['minFillerPercentMatch'] else (1 if fitD >= fitG else 2)
    # -[NJoin doAlign]
    if fill_len == 0: fill_loc2 = 0
    if micro_len == 0: micro_loc = 0
    lo = max(0, fill_loc2); hi = min(len(m), fill_loc2 + fill_len)
    filler = m[lo:hi] if hi > lo else ''
    nfill = max(0, hi - lo)
    micro = m[micro_loc:micro_loc + micro_len]
    oh = ''
    parts = original.split('^')
    if len(parts) == 3: oh = parts[1].upper()
    loh = T.rfind(oh) if oh else -1
    roh = Bb.find(oh) if oh else -1
    return dict(complex=int(complexA), microhomology=micro, micro_length=micro_len,
                imperfect_micro=int(imperfect), filler=filler, filler_length=nfill,
                filler_template=templ, missing_left=ml_len, missing_right=mr_len,
                rows=(t, ''.join(mt), m, ''.join(mb), b),
                ranges=dict(micro=(micro_loc, micro_len), filler=(lo, nfill),
                            mleft=(ml_loc, ml_len), mright=(mr_loc, mr_len),
                            loh=(max(loh, 0), len(oh) if loh >= 0 else 0),
                            roh=(max(roh, 0), len(oh) if roh >= 0 else 0)))


def read_input(path):
    txt = open(path, encoding='utf-8', errors='replace').read().replace('\r\n', '\n').replace('\r', '\n')
    left = right = overhang = original = None
    joins = []
    for line in txt.split('\n'):
        items = [x for x in re.split(r'[\t,]+', line.strip()) if x.strip()]
        items = [x.strip() for x in items]
        if len(items) < 2: continue
        name, seq = items[0], items[1].replace(' ', '')
        U = name.upper()
        if U == 'LEFT': left = seq
        elif U == 'RIGHT': right = seq
        elif U == 'OVERHANG': overhang = seq
        elif U == 'ORIGINAL': original = seq
        elif '//' in name[:6]: continue
        else: joins.append((name, seq))
    if original is None:
        oh = overhang or ''
        original = left[:len(left) - len(oh)] + '^' + oh + '^' + right[len(oh):]
    parts = original.split('^')
    if len(parts) == 2: left, right = parts
    elif len(parts) == 3: left, right = parts[0] + parts[1], parts[1] + parts[2]
    return original, left, right, joins


def diagram(name, r, maxlen=200, per_side=30):
    """The five-line text display of -[NJoin doAlign] (abridged when long)."""
    t, mt, m, mb, b = r['rows']
    rows = [t, mt, m, mb, b]
    if len(t) > maxlen:
        r1 = t.find(' '); r1 = r1 if r1 >= 0 else len(t) - 2
        r2 = re.search('[ATCGatcg]', b); r2 = r2.start() if r2 else len(b) - 2
        mid1, mid2 = min(r1, r2), max(r1, r2)
        Lr = (3, max(mid1 - 3 - per_side, 0)); mid2 -= Lr[1]
        rl = min(mid2 + per_side, len(t) - Lr[1]); rlen = len(t) - Lr[1] - rl - 2
        def ab(s):
            if Lr[1] > 0: s = s[:Lr[0]] + '…' + s[Lr[0] + Lr[1]:]
            if 0 < rlen <= len(s) - rl: s = s[:rl] + '…' + s[rl + rlen:]
            return s
        rows = [ab(x) for x in rows]
    return '\n'.join([name] + rows)


CSS = """body{font-family:-apple-system,Helvetica,sans-serif;margin:16px;color:#222;background:#fff}
pre{font-family:Menlo,Consolas,monospace;font-size:13px;line-height:1.25;margin:4px 0 18px;overflow-x:auto}
.n{font-weight:600;font-family:-apple-system,Helvetica,sans-serif}
.miss{color:#999;text-decoration:line-through}.mh{color:#c00;font-weight:700}
.fill{background:#cfe3ff}.oh{text-decoration:underline}.key span{padding:0 4px}
@media (prefers-color-scheme: dark){body{background:#1b1b1b;color:#ddd}.fill{background:#284870}.miss{color:#777}.mh{color:#ff6b6b}}"""


def html_rows(r):
    t, mt, m, mb, b = r['rows']; R = r['ranges']
    def mark(s, specs):
        cls = [set() for _ in s]
        for c, (lo, n) in specs:
            for i in range(lo, min(lo + n, len(s))): cls[i].add(c)
        out = []; cur = None; buf = ''
        for ch_, c in zip(s, cls):
            k = ' '.join(sorted(c))
            if k != cur:
                if buf: out.append(('<span class="%s">%s</span>' % (cur, html.escape(buf))) if cur else html.escape(buf))
                cur, buf = k, ''
            buf += ch_
        if buf: out.append(('<span class="%s">%s</span>' % (cur, html.escape(buf))) if cur else html.escape(buf))
        return ''.join(out)
    return '\n'.join([mark(t, [('miss', R['mleft']), ('oh', R['loh'])]), html.escape(mt),
                      mark(m, [('mh', R['micro']), ('fill', R['filler'])]), html.escape(mb),
                      mark(b, [('miss', R['mright']), ('oh', R['roh'])])])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('input'); ap.add_argument('-o', '--out')
    ap.add_argument('--diagrams'); ap.add_argument('--html')
    ap.add_argument('--params', default='poster2011', choices=sorted(PARAMS))
    a = ap.parse_args()
    P = PARAMS[a.params]
    original, left, right, joins = read_input(a.input)
    rows = []; diags = []; hrows = []
    for name, seq in joins:
        if len(seq) <= 5: continue
        r = align(seq, left, right, original, P)
        rows.append([name, seq, r['complex'], r['microhomology'], r['micro_length'], r['imperfect_micro'],
                     r['filler'], r['filler_length'], r['filler_template'], r['missing_left'], r['missing_right']])
        diags.append(diagram(name, r))
        info = 'micro %s (%d)  filler %s (%d)  missing L %d / R %d%s' % (
            r['microhomology'] or '-', r['micro_length'], r['filler'] or '-', r['filler_length'],
            r['missing_left'], r['missing_right'], '  complex' if r['complex'] else '')
        hrows.append('<div class="n">%s <small>%s</small></div><pre>%s</pre>' % (html.escape(name), html.escape(info), html_rows(r)))
    hdr = ['name', 'sequence', 'complex', 'microhomology', 'micro length', 'imperfect micro', 'filler',
           'filler length', 'filler template', 'missing left', 'missing right']
    out = open(a.out, 'w', newline='') if a.out else sys.stdout
    w = csv.writer(out); w.writerow(hdr); w.writerows(rows)
    if a.html:
        key = ('<p class="key">Key: <span class="miss">missing bases</span> <span class="mh">microhomology</span> '
               '<span class="fill">filler</span> <span class="oh">overhang</span> &nbsp; rows: left end / match / product / match / right end</p>')
        open(a.html, 'w').write('<!doctype html><meta charset="utf-8"><title>NJoiner alignments</title><style>%s</style>'
                                '<h2>NJoiner alignments: %s</h2><p>Original: <code>%s</code> &nbsp; parameters: %s</p>%s%s'
                                % (CSS, html.escape(a.input.split('/')[-1]), html.escape(original), a.params, key, ''.join(hrows)))
    if a.diagrams:
        open(a.diagrams, 'w').write('Original: %s\n\n' % original + '\n\n'.join(diags) + '\n')


if __name__ == '__main__':
    main()
