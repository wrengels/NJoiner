#!/usr/bin/env python3
"""
validate.py -- checks the Python port of NJoiner against the surviving record.

1. Row-by-row comparison with "test/No109Rr3EJ output (NJoiner 2010).csv", the
   output NJoiner itself wrote on 12 April 2010 for test/No109Rr3EJ.txt.
2. Re-derivation of the junction classes behind Table 3 of the Position Effects
   manuscript from the plate files in data/NHEJ_junction_sequences/, compared
   with the authors' summary tabulation (transcribed below from
   "NHEJ plates 1-8 NJoiner.xlsx").

Run from this folder:  python3 python/validate.py
Writes test/No109Rr3EJ comparison.csv and table3/*.csv, *.html.
"""
import csv, os, re, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import njoiner as N

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(HERE)

# ---------- 1. No109 ----------
def load_orig(fn):
    txt = open(fn, encoding='utf-8', errors='replace').read().replace('\r', '\n')
    rows = collections.defaultdict(list)
    for r in csv.reader(txt.splitlines()):
        if not r or r[0] == 'name': continue
        r = [x.lstrip("'").strip() for x in r]
        rows[r[0]].append(r)
    return rows

orig_rows = load_orig('test/No109Rr3EJ output (NJoiner 2010).csv')
original, left, right, joins = N.read_input('test/No109Rr3EJ.txt')
P = N.PARAMS['output2010']
fields = ['complex', 'microhomology', 'micro length', 'imperfect micro', 'filler', 'filler length',
          'filler template', 'missing left', 'missing right']
seen = collections.Counter(); out = []; nrow = nsame = ninv = 0
for name, seq in joins:
    if len(seq) <= 5: continue
    o = orig_rows[name][seen[name]]; seen[name] += 1
    r = N.align(seq, left, right, original, P)
    port = [r['complex'], r['microhomology'], r['micro_length'], r['imperfect_micro'], r['filler'],
            r['filler_length'], r['filler_template'], r['missing_left'], r['missing_right']]
    old = o[2:11]
    agree = [str(a) == b for a, b in zip(port, old)]
    # the 2010 build wrote 1 in "imperfect micro" for a PERFECT microhomology
    perfect_flag = int(r['micro_length'] > 0 and not r['imperfect_micro'])
    inv_ok = str(perfect_flag) == old[3]
    nrow += 1
    if all(a for i, a in enumerate(agree) if i != 3): nsame += 1
    if inv_ok: ninv += 1
    out.append([name, seq] + [x for pair in zip(old, port) for x in pair] + [int(all(a for i, a in enumerate(agree) if i != 3)), int(inv_ok)])
hdr = ['name', 'sequence'] + [f'{f} ({w})' for f in fields for w in ('NJoiner 2010', 'port')] + \
      ['all other fields agree', '"imperfect micro" = perfect-microhomology flag']
with open('test/No109Rr3EJ comparison.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(hdr); w.writerows(out)
print(f'No109: {nrow} junctions; all fields except "imperfect micro" agree in {nsame}; '
      f'"imperfect micro" equals a perfect-microhomology flag in {ninv}.')

# ---------- 2. Table 3 ----------
PLATES = [(1, '21C, young'), (2, '21C, old'), (3, '82D, young'), (6, '48C, young'),
          (8, '48C, 3 weeks'), (7, '48C, 4-6 weeks'), (5, '48C, held to old age')]
TAB = {1: (38, 9, [5, 5, 3, 2, 11, 0, 0, 0, 3]), 2: (20, 6, [7, 0, 4, 1, 1, 0, 0, 0, 1]),
       3: (13, 4, [6, 1, 1, 0, 1, 0, 0, 0, 0]), 6: (91, 28, [11, 17, 6, 1, 22, 0, 2, 0, 4]),
       8: (89, 16, [17, 13, 10, 4, 22, 0, 1, 1, 5]), 7: (82, 23, [23, 10, 6, 1, 12, 1, 0, 0, 6]),
       5: (78, 23, [30, 9, 2, 1, 12, 0, 0, 1, 0])}
P = N.PARAMS['poster2011']
os.makedirs('table3/alignments', exist_ok=True)
per = []; summ = []
for p, lab in PLATES:
    fn = f'data/NHEJ_junction_sequences/Rr3 NHEJ CRP {p} sib names.txt'
    original, left, right, joins = N.read_input(fn)
    n = indel = 0; mh = [0] * 9; aggg = ataa = 0
    for name, seq in joins:
        m = re.fullmatch(r'LongDel(Odd)?MH(\d)', seq)
        if seq == 'LongDelPlusInsert':
            n += 1; indel += 1; per.append([lab, name, seq, 'in-del (hand-scored long deletion)', '', '', '', '', '']); continue
        if m:
            n += 1; mh[int(m.group(2))] += 1
            per.append([lab, name, seq, 'EJ (hand-scored long deletion)', '', m.group(2), '', '', '']); continue
        if not re.fullmatch('[ACGTNacgtn]+', seq) or len(seq) <= 5:
            per.append([lab, name, seq, 'not analysed (not a clean sequence)', '', '', '', '', '']); continue
        r = N.align(seq, left, right, original, P); n += 1
        cls = 'in-del' if r['filler_length'] > 0 else 'EJ'
        if cls == 'in-del': indel += 1
        else:
            mh[min(r['micro_length'], 8)] += 1
            if r['micro_length'] == 4 and r['microhomology'].upper() in ('AGGG', 'AGNG'): aggg += 1
            if r['micro_length'] == 4 and r['microhomology'].upper() == 'ATAA': ataa += 1
        per.append([lab, name, seq, cls, r['microhomology'], r['micro_length'], r['filler'],
                    r['missing_left'], r['missing_right']])
    T = TAB[p]; ej = n - indel; tej = T[0] - T[1]
    pct = lambda a, b: round(100 * a / b, 1) if b else ''
    summ.append([lab, T[0], n, T[1], indel, tej, ej, tej - T[2][0], ej - mh[0],
                 pct(tej - T[2][0], tej), pct(ej - mh[0], ej), T[2][4], mh[4], pct(T[2][4], tej), pct(mh[4], ej), aggg, ataa])
    os.system(f'python3 python/njoiner.py "{fn}" --params poster2011 -o /dev/null '
              f'--html "table3/alignments/{lab.replace(",", "")}.html"')
with open('table3/per-junction classes (port, poster parameters).csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['sample', 'name', 'sequence', 'class', 'microhomology', 'micro length', 'filler', 'missing left', 'missing right'])
    w.writerows(per)
with open('table3/Table 3 re-derived vs tabulation.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['sample', 'junctions (tab)', 'junctions (port)', 'in-del (tab)', 'in-del (port)', 'EJ (tab)', 'EJ (port)',
                'EJ with mh (tab)', 'EJ with mh (port)', '% of EJ (tab)', '% of EJ (port)', '4-bp mh (tab)', '4-bp mh (port)',
                '4-bp % of EJ (tab)', '4-bp % of EJ (port)', '4-bp = AGGG repeat (port)', '4-bp = ATAA overhang (port)'])
    w.writerows(summ)
for s in summ: print(s)
