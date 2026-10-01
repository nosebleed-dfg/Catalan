r"""flip_families.py (2026-10-01, real date): "square each face, the carry is the triangle left over", made exact.

Every torsion-free genus-0 subgroup of index 12 in PSL2(Z) is an ideal triangulation of the 4-punctured sphere:
12 darts, u = rotation inside a triangle (4 triangles), s = the other side of an edge (6 edges), and the cusps are the
cycles of T = S U (first s, then u); the cusp width is the vertex degree, so the widths add up to 12 and each ideal
triangle has area pi.  A FLIP replaces one edge PQ by the other diagonal RS of its quadrilateral: P and Q lose one
triangle corner each, R and S gain one.  "Squaring a pentagon and carrying the triangle" is a flip.
This script
  1. lists all such triangulations up to relabelling (= subgroups up to conjugacy), with their cusp widths;
  2. names the congruence ones by building the coset actions of Gamma(3), Gamma_0(4) n Gamma^0(2), Gamma_1(5),
     Gamma_0(6), Gamma_0(8), Gamma_0(9) (Beauville's six families);
  3. builds the flip graph, and prints for each triangulation the edge products width(x)*width(y) (the invariant of the
     two-parabolic group of an adjacent cusp pair).
Usage: python flip_families.py"""
import itertools
from collections import deque

N = 12
u0 = tuple((i//3)*3 + (i + 1) % 3 for i in range(N))

def cycles(p):
    seen = [False]*len(p); out = []
    for i in range(len(p)):
        if not seen[i]:
            c = []; j = i
            while not seen[j]: seen[j] = True; c.append(j); j = p[j]
            out.append(c)
    return out
def comp(p, q):                     # first p, then q
    return tuple(q[p[i]] for i in range(len(p)))
def canon(s, u):
    best = None
    for start in range(N):
        lab = {start: 0}; order = [start]; dq = deque([start])
        while dq:
            x = dq.popleft()
            for g in (s, u):
                y = g[x]
                if y not in lab: lab[y] = len(order); order.append(y); dq.append(y)
        if len(order) < N: return None
        key = (tuple(lab[s[order[i]]] for i in range(N)), tuple(lab[u[order[i]]] for i in range(N)))
        if best is None or key < best: best = key
    return best
def widths(s, u):
    return tuple(sorted((len(c) for c in cycles(comp(s, u))), reverse=True))
def involutions(pts):
    if not pts: yield {}; return
    x = pts[0]
    for y in pts[1:]:
        rest = [z for z in pts[1:] if z != y]
        for d in involutions(rest): yield {**d, x: y, y: x}

classes = {}
for d in involutions(list(range(N))):
    s = tuple(d[i] for i in range(N))
    if len(cycles(comp(s, u0))) != 4: continue
    key = canon(s, u0)
    if key is None: continue
    classes.setdefault(key, widths(*key))
print("1. torsion-free genus-0 subgroups of index 12, up to conjugacy: %d" % len(classes))

# ---- named congruence groups ----
def mat_mul(A, B, m):
    return ((A[0]*B[0] + A[1]*B[2]) % m, (A[0]*B[1] + A[1]*B[3]) % m, (A[2]*B[0] + A[3]*B[2]) % m, (A[2]*B[1] + A[3]*B[3]) % m)
S = (0, -1, 1, 0); U = (0, -1, 1, 1)
def perm_from_action(points, act):
    idx = {p: i for i, p in enumerate(points)}
    return tuple(idx[act(p, S)] for p in points), tuple(idx[act(p, U)] for p in points)
def p1(m):
    pts = []
    for c in range(m):
        for d in range(m):
            from math import gcd
            if gcd(gcd(c, d), m) != 1: continue
            # normalise by units
            rep = min(((c*k) % m, (d*k) % m) for k in range(1, m) if gcd(k, m) == 1)
            if rep not in pts: pts.append(rep)
    def act(p, X):
        from math import gcd
        c, d = (p[0]*X[0] + p[1]*X[2]) % m, (p[0]*X[1] + p[1]*X[3]) % m
        return min(((c*k) % m, (d*k) % m) for k in range(1, m) if gcd(k, m) == 1)
    return pts, act
def rows_pm(m):
    pts = []
    for c in range(m):
        for d in range(m):
            if (c, d) == (0, 0): continue
            rep = min((c, d), ((-c) % m, (-d) % m))
            if rep not in pts: pts.append(rep)
    def act(p, X):
        c, d = (p[0]*X[0] + p[1]*X[2]) % m, (p[0]*X[1] + p[1]*X[3]) % m
        return min((c, d), ((-c) % m, (-d) % m))
    return pts, act
def psl_cosets(m, inH):
    # all of SL2(Z/m), cosets H g with H = {h : inH(h)}, right action
    els = [(a, b, c, d) for a in range(m) for b in range(m) for c in range(m) for d in range(m) if (a*d - b*c) % m == 1]
    H = [h for h in els if inH(h)]
    def coset(g): return min(mat_mul(h, g, m) for h in H)
    pts = sorted(set(coset(g) for g in els))
    def act(p, X): return coset(mat_mul(p, tuple(x % m for x in X), m))
    return pts, act
named = {}
for name, (pts, act) in (("Gamma_0(6)", p1(6)), ("Gamma_0(8)", p1(8)), ("Gamma_0(9)", p1(9)), ("Gamma_1(5)", rows_pm(5)),
                         ("Gamma(3)", psl_cosets(3, lambda h: h in ((1, 0, 0, 1), (2, 0, 0, 2)))),
                         ("Gamma_0(4) n Gamma^0(2)", psl_cosets(4, lambda h: h[2] == 0 and h[1] % 2 == 0))):
    assert len(pts) == N, (name, len(pts))
    s, u = perm_from_action(pts, act)
    key = canon(s, u)
    assert key in classes, name
    named[key] = name
    print("   %-24s widths %s" % (name, widths(s, u)))

# ---- flips ----
def flips(s, u):
    out = []
    for a in range(N):
        b = s[a]
        if a > b: continue
        a1, a2 = u[a], u[u[a]]; b1, b2 = u[b], u[u[b]]
        if b in (a, a1, a2): continue                      # folded triangle: no flip
        v = list(u)
        v[b2] = a1; v[a1] = a; v[a] = b2; v[a2] = b1; v[b1] = b; v[b] = a2
        out.append(canon(s, tuple(v)))
    return out
keys = sorted(classes, key=lambda k: (classes[k], k))
label = {}
count = {}
for k in keys:
    w = classes[k]; count[w] = count.get(w, 0) + 1
    label[k] = "%s%s%s" % ("".join(str(x) for x in w) if max(w) < 10 else "-".join(str(x) for x in w),
                           "" if count[w] == 1 else chr(ord('a') + count[w] - 1), (" [" + named[k] + "]") if k in named else "")
adj = {k: sorted(set(label[f] for f in flips(*k) if f is not None and f != k)) for k in keys}
print("\n2. classes, edge products, and flips")
for k in keys:
    s, u = k
    t = comp(s, u); cyc = cycles(t); deg = {x: len(c) for c in cyc for x in c}
    prods = sorted(set(deg[a]*deg[s[a]] for a in range(N) if deg[a] != 0 and not any(a in c and s[a] in c for c in cyc)))
    print("   %-36s edge products %-22s flips to: %s" % (label[k], prods, ", ".join(adj[k])))

# ---- distances from Apery's family ----
start = [k for k in keys if named.get(k) == "Gamma_1(5)"][0]
dist = {label[start]: 0}; dq = deque([start])
lab2key = {label[k]: k for k in keys}
while dq:
    k = dq.popleft()
    for nb in adj[k]:
        if nb not in dist: dist[nb] = dist[label[k]] + 1; dq.append(lab2key[nb])
print("\n3. flip distance from Apery's family (5,5,1,1):")
for lab, dd in sorted(dist.items(), key=lambda kv: (kv[1], kv[0])): print("   %d  %s" % (dd, lab))
