#!/usr/bin/env python3
# Genere le bloc de donnees de la Figure 1 de l'article a partir de
# fullspgemm.csv (medianes des essais par beta). A executer par les auteurs ;
# coller la sortie a la place des trois \addplot de la figure.
# Usage : python3 make_figure1.py fullspgemm.csv
import csv, sys

def med(xs):
    xs = sorted(xs)
    return xs[len(xs) // 2]

rows = list(csv.DictReader(open(sys.argv[1])))
betas = []
for r in rows:
    if r['kind'] == 'beta' and r['param'] not in betas:
        betas.append(r['param'])

def series(col):
    out = []
    for b in betas:
        rs = [r for r in rows if r['kind'] == 'beta' and r['param'] == b]
        out.append((b, med([float(r[col]) for r in rs])))
    return out

sl  = series('t_stateless_ms')
gd  = series('t_guarded_ms')
slg = series('t_stateless_gallop_ms')
rr  = series('rewind_rate')

def fmt(pairs):
    return ' '.join(f"({b},{v:.3f})" for b, v in pairs)

print("\\addplot+[mark=o] coordinates")
print("  {" + fmt(rr) + "};")
print("\\addlegendentry{rewind rate}")
print("\\addplot+[mark=square*] coordinates")
print("  {" + fmt([(b, s / g) for (b, s), (_, g) in zip(sl, gd)]) + "};")
print("\\addlegendentry{stateless / persistent}")
print("\\addplot+[mark=triangle*] coordinates")
print("  {" + fmt([(b, s / g) for (b, s), (_, g) in zip(sl, slg)]) + "};")
print("\\addlegendentry{plain / galloping (stateless)}")
