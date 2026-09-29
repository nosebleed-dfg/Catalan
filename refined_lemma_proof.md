# A refined arithmetic lemma for Zudilin's linear forms in even β-values, with proof

Companion to the results log section "The multi-β route" (dedekind_carry_results-1.md) and to multibeta.py. Written 2026-09-28.

**Priority note (added later on 2026-09-28).** The statement of Theorem 3 below — at least one of β(2), β(4), β(6), β(8), β(10) is irrational — was proved by Li Lai and Li Zhou, arXiv 2103.00904 (2021), Theorem 6.1, with Zudilin's construction, s = 11, η = (94; 32,32,32,32,33,34,35,36,37,38,39) and a different normalization of the numerator (η₀ blocks of length n), margin ≈ 0.21 per n. Normalizations only rescale the forms and cannot change the primitive integer form; the refined bound below is covariant under them and gives their direction a margin of 6.21 per n. The contribution of this note is the arithmetic (Theorems 1, 4, 5, 6), not the five-value statement, which it recovers with a simpler direction and a larger margin.

**Summary.** For the very-well-poised linear forms r_n = Σ_{i even} a_i β(i) + a_0 of W. Zudilin ("Arithmetic of Catalan's constant and its relatives", arXiv 1804.09922; Lemmas 4–5 there) we prove a p-adic lower bound for the coefficients that is sharper than his Φ_n^{-1}d_M^{s-i} a_i ∈ Z at the primes p in (M/2, M]. The two new ingredients are elementary: (i) a partial sum Σ(−1)^ℓ(ℓ+½)^{−i} carries p in its denominator only at poles at distance ≥ (p+1)/2 from the centre, and those poles have low multiplicity; (ii) a Laurent coefficient at a pole can lose a factor p only through a factor of the rational function divisible by p, and when no denominator factor at that pole is divisible by p the loss is capped by the number of p-divisible numerator factors. With the exact Legendre count of the top coefficient this gives Theorem 1; Theorem 2 turns it into a denominator sequence D_n with (log D_n)/n → s·M/n − κ_ref; Theorem 3 applies it to η = (33; 9,10,10,11,11,12,12,12,13,13,14), s = 11, where the decay 76.8118 per n beats the height 76.1255 per n: **at least one of β(2), β(4), β(6), β(8), β(10) is irrational** (Zudilin: one of β(2), …, β(12)). What is cited from Zudilin: Lemma 2 (integral representation, positivity, asymptotics), Lemma 3 (the saddle), and the first part of Lemma 4 (d_M^{s−i}a_{i,k} ∈ Z) for the primes p ≤ √(2h₀) only, whose total contribution is O(√n). His cancellation factor Φ_n is not used; it is recovered as a special case (Remark 4).

Notation: v_p is the p-adic valuation; ⌊·⌋ the floor; [P] = 1 if P holds, 0 otherwise. Throughout p is an **odd prime with p² > 2h₀** unless said otherwise.

---

## 1. Setting (Zudilin's construction, his §3)

Let s ≥ 3 be odd and η = (η₀; η₁, …, η_s) integers with 0 < η_j < η₀/2 and η₁ + … + η_s ≤ (s−1)η₀/2. For n ≥ 1 with η₀n even put

- h₀ = η₀n + 1 (odd), m = (h₀−1)/2 = η₀n/2 (the centre), h_j* = η_jn (so h_j = h_j* + ½), L_j = h₀ − 2h_j* = (η₀−2η_j)n + 1, N = min_j h_j* = η_min n,
- R_n(t) = γ_n · (2t+h₀) · ∏_{i=1}^{h₀−1}(t+i) / ∏_{j=1}^{s} ∏_{i=0}^{L_j−1}(t + h_j* + ½ + i),  γ_n = 4^{h₀−1} ∏_{j=2}^{s}(L_j−1)! / (h₁*!)²,
- r_n = Σ_{ν≥0} (−1)^ν R_n(ν).

Degree count: deg(denominator) − deg(numerator) = (s−1)h₀ − 2nΣη_j ≥ s−1 ≥ 2, so R_n(t) = O(t^{−2}) and everything below converges absolutely except the i = 1 tails, which are alternating with monotone terms.

**Poles.** R_n has poles only at t = −(k+½), k ∈ K := {N, N+1, …, h₀−1−N}. The factor j covers k iff h_j* ≤ k ≤ h₀−1−h_j*; write J(k) for the set of such j. The pole order is s_k = |J(k)| − [k = m] (the factor 2t+h₀ = 2(t + h₀/2) vanishes at the centre). Always s_k ≥ 1, s_m = s−1.

**Palindrome.** R_n(−t−h₀) = R_n(t) (each of the s+1 factor blocks is reversed by t ↦ −t−h₀ up to a sign, and the signs multiply to (−1)^{1+ΣL_j} = +1 since ΣL_j ≡ sh₀ ≡ 1 mod 2).

**Partial fractions.** R_n(t) = Σ_{k∈K} Σ_{i=1}^{s_k} a_{i,k} (t+k+½)^{−i} with a_{i,k} ∈ Q (no polynomial part). The palindrome gives a_{i,k} = (−1)^i a_{i,h₀−1−k}.

**The linear form.** R_n vanishes at ν = −1, …, −(h₀−1), so r_n = Σ_{ν≥−m}(−1)^νR_n(ν). Substituting the partial fractions and ℓ = ν + k,

  r_n = Σ_{i} 2^i β(i) · Σ_k (−1)^k a_{i,k} + a₀,   with  a_i := 2^i Σ_k (−1)^k a_{i,k},  a₀ := Σ_{k∈K} Σ_{i=1}^{s_k} (−1)^k a_{i,k} S_{i,k},

where β(i) = Σ_{ℓ≥0}(−1)^ℓ(2ℓ+1)^{−i} and the partial sums are

  S_{i,k} = Σ_{ℓ=k−m}^{−1} (−1)^ℓ (ℓ+½)^{−i}  (k < m),   S_{i,k} = −Σ_{ℓ=0}^{k−m−1} (−1)^ℓ (ℓ+½)^{−i}  (k > m),   S_{i,m} = 0.

(Here Σ_{ℓ≥k−m} = Σ_{ℓ≥0} + Σ_{k−m≤ℓ≤−1} for k < m and Σ_{ℓ≥0} − Σ_{0≤ℓ≤k−m−1} for k > m; the rearrangement is a finite sum of convergent series.) Since h₀−1 is even, k and h₀−1−k have the same parity, so a_i = (−1)^i a_i and **a_i = 0 for odd i**: r_n ∈ Q + Qβ(2) + Qβ(4) + … + Qβ(s−1). This is Zudilin's Lemma 5; the rest of this note is about denominators.

---

## 2. The local structure at a pole

Fix k ∈ K and put t = −(k+½) + u. Each linear factor of R_n is (value at u = 0) × (1 + u/value) when the value is nonzero. The nonzero values are:

- A_k (numerator): α_i = i − k − ½ = (2i−2k−1)/2 for 1 ≤ i ≤ h₀−1 (half-integers, never 0), and, if k ≠ m, α = m − k from the factor 2t+h₀ = 2(m−k)(1 + u/(m−k)); at the centre this factor is 2u.
- B_k (denominator): β = h_j* + i − k for 1 ≤ j ≤ s, 0 ≤ i ≤ L_j−1, β ≠ 0 (integers; the zero values are exactly the |J(k)| poles).

Then

  R_n(t) = u^{−s_k} B_k(u),  B_k(u) = C_k ∏_{α∈A_k}(1 + u/α) / ∏_{β∈B_k}(1 + u/β),  C_k = 2γ_n ∏_{α∈A_k} α / ∏_{β∈B_k} β,

and  **a_{i,k} = [u^{s_k−i}] B_k(u)** for 1 ≤ i ≤ s_k; in particular a_{s_k,k} = C_k (the top coefficient).

Sizes: |2i−2k−1| ≤ 2h₀−3, |m−k| < h₀, |β| ≤ h₀−1−2N < h₀. Hence for an odd prime with p² > 2h₀ every element of A_k ∪ B_k has p-adic valuation 0 or 1 (for α = o/2 with o odd, v_p(α) = v_p(o)).

### Lemma A (derivative saturation)

Let p be an odd prime with p² > 2h₀. Put A_k^p = {α ∈ A_k : v_p(α) = 1}, B_k^p = {β ∈ B_k : v_p(β) = 1}, and

  D_k(p) = |A_k^p| if B_k^p = ∅,  D_k(p) = +∞ otherwise.

Then for 0 ≤ r ≤ s_k − 1:  v_p([u^r] B_k(u)) ≥ v_p(C_k) − min(r, D_k(p)).  Equivalently, **v_p(a_{i,k}) ≥ v_p(C_k) − min(s_k − i, D_k(p))** for 1 ≤ i ≤ s_k.

*Proof.* Expanding the product, [u^r] B_k/C_k = Σ ∏_{α∈S} α^{−1} ∏_{β∈B_k} (−β^{−1})^{r_β}, the sum over subsets S ⊆ A_k and exponents r_β ≥ 0 with |S| + Σ_β r_β = r. A term has valuation −|S ∩ A_k^p| − Σ_{β∈B_k^p} r_β. This is ≥ −(|S| + Σr_β) = −r always, and if B_k^p = ∅ it equals −|S ∩ A_k^p| ≥ −|A_k^p|. The valuation of a finite sum is at least the minimum over its terms. ∎

Remark. The bound −r is Zudilin's "each derivative costs d_M"; the cap |A_k^p| is new. When B_k^p = ∅ the p-singular part of B_k/C_k is the polynomial ∏_{α∈A_k^p}(1+u/α) of degree |A_k^p|, so no Taylor coefficient can carry more than |A_k^p| inverse powers of p.

### Lemma B (the top coefficient is an exact carry count)

Let p be an odd prime with p² > 2h₀, x = n/p, y = k/p (k ∈ K). Then

  v_p(C_k) = Σ_{j=2}^{s} ⌊(η₀−2η_j)x⌋ − 2⌊η₁x⌋ + ⌊η₀x − y + ½⌋ + ⌊y + ½⌋ − Σ_{j=1}^{s} ( ⌊(η₀−η_j)x − y⌋ − ⌊η_jx − y⌋ + [p | h_j*−k] − [j ∈ J(k)] ) + [k≠m]·[p | m−k].

Consequently, with Zudilin's carry function

  φ(x,y) = ⌊2(η₀x−y)⌋ + ⌊2y⌋ − ⌊η₀x−y⌋ − ⌊y⌋ − 2⌊η₁x⌋ − ⌊(η₀−2η₁)x⌋ + Σ_{j=1}^{s} ( ⌊(η₀−2η_j)x⌋ − ⌊y−η_jx⌋ − ⌊(η₀−η_j)x−y⌋ ),

one has the exact identity

  **Φ_p(k) := v_p(C_k) + (s − s_k) = φ(n/p, k/p) + ε_k,  ε_k := [k = m] + [k ≠ m]·[p | m−k] ∈ {0, 1}.**

In particular v_p(C_k) ≥ φ(n/p, k/p) − (s − s_k), with equality unless k ≡ m (mod p).

*Proof.* Since p² exceeds every factor, v_p of a product of integers from an interval is the number of multiples of p in it, and v_p(N!) = ⌊N/p⌋ for N < p².

(a) γ_n: L_j − 1 = (η₀−2η_j)n and h₁* = η₁n, so v_p(γ_n) = Σ_{j≥2}⌊(η₀−2η_j)x⌋ − 2⌊η₁x⌋; 4^{h₀−1} and the constant 2 are units.

(b) The factor 2t+h₀ contributes v_p(2(m−k)) = [p | m−k] when k ≠ m, and 0 at the centre.

(c) Numerator: the values 2α_i = 2i−2k−1 run through the odd integers from 1−2k to 2h₀−3−2k. An odd integer is a multiple of p iff it is pq with q odd, so the count is N_k = #{q odd : (1−2k)/p ≤ q ≤ (2h₀−3−2k)/p}. For reals A ≤ B, #{q odd : A ≤ q ≤ B} = ⌊(B+1)/2⌋ − ⌊(A+1)/2⌋ + [A ∈ 2Z+1]. Here (B+1)/2 = (2h₀−3−2k+p)/(2p) and (A+1)/2 = (p+1−2k)/(2p). Three simplifications:
 (i) ⌊(2h₀−3−2k+p)/(2p)⌋ = ⌊(2h₀−2−2k+p)/(2p)⌋: an integer M with 2pM = 2h₀−2−2k+p would give p(2M−1) = 2(h₀−1−k), odd = even; and (2h₀−2−2k+p)/(2p) = η₀x − y + ½ because 2h₀−2 = 2η₀n.
 (ii) −⌊(p+1−2k)/(2p)⌋ + [p | 2k−1] = −⌊(p−2k)/(2p)⌋: if p ∤ 2k−1 the two floors agree (an integer M with 2pM = p+1−2k forces p | 2k−1); if 2k−1 = pq with q odd, then (p+1−2k)/(2p) = (1−q)/2 ∈ Z while (p−2k)/(2p) = (1−q)/2 − 1/(2p), so the floors differ by one and the indicator compensates.
 (iii) −⌊(p−2k)/(2p)⌋ = −⌊½ − y⌋ = ⌊y + ½⌋, because y + ½ = (2k+p)/(2p) is never an integer (2k+p is odd).
 Hence N_k = ⌊η₀x − y + ½⌋ + ⌊y + ½⌋.

(d) Denominator, factor j: the values β = h_j* + i − k, 0 ≤ i ≤ L_j−1, run through the integers of [h_j*−k, h₀−1−h_j*−k], and the nonzero multiples of p among them number
 M_{j,k} = ⌊(h₀−1−h_j*−k)/p⌋ − ⌊(h_j*−k−1)/p⌋ − [j ∈ J(k)] = ⌊(η₀−η_j)x − y⌋ − ⌊η_jx − y⌋ + [p | h_j*−k] − [j ∈ J(k)],
 using h₀−1−h_j* = (η₀−η_j)n and ⌊z − 1/p⌋ = ⌊z⌋ − [z ∈ Z] for z = (h_j*−k)/p (its fractional part is a multiple of 1/p).

Adding (a)+(b)+(c)−(d) gives the formula for v_p(C_k). For the identity: s − s_k = s − |J(k)| + [k=m], and −⌊y − η_jx⌋ = ⌊η_jx − y⌋ + 1 − [η_jx − y ∈ Z] = ⌊η_jx − y⌋ + 1 − [p | h_j*−k], so Σ_j(1 − [p | h_j*−k] + ⌊η_jx−y⌋) = −Σ_j⌊y−η_jx⌋; and ⌊2z⌋ − ⌊z⌋ = ⌊z+½⌋ turns ⌊2(η₀x−y)⌋ + ⌊2y⌋ − ⌊η₀x−y⌋ − ⌊y⌋ into ⌊η₀x−y+½⌋ + ⌊y+½⌋. Collecting terms, v_p(C_k) + s − s_k = φ(x,y) + [k=m] + [k≠m][p | m−k]. ∎

*Numerical confirmation of the identity.* The exact top coefficients a_{s_k,k} were computed for (η, n) = ((31;10⁵,11⁴,12⁴), 8), ((33;9,…,14), 8), ((33;9,…,14), 12), ((3;1⁷), 80), ((5;2⁵), 20), ((36;10,…,15), 5) and the identity Φ_p(k) = φ(n/p,k/p) + ε_k was checked at every pole k and every odd prime √(2h₀) < p ≤ M: 12062 (prime, pole) pairs, 0 violations (script lemmaB_check.py, session tmp directory).

Remark (the tower). s − s_k is the number of floors missing from the tower of multiplicities at pole k. Zudilin's estimate v_p(a_{i,k}) ≥ −(s−i) + φ is his Lemma 4; Lemma B says his φ is exactly the top coefficient's carry count plus one unit per missing floor, and the missing floors are precisely the derivatives that are never taken at that pole — the slack that Lemma A cashes in.

### Lemma C (the partial-sum window)

Let p be an odd prime with p² > 2h₀. Then v_p(S_{i,k}) ≥ −i if |k − m| ≥ (p+1)/2, and v_p(S_{i,k}) ≥ 0 otherwise.

*Proof.* (ℓ+½)^{−i} = 2^i (2ℓ+1)^{−i}, and the odd integers |2ℓ+1| occurring in S_{i,k} are 1, 3, …, 2|k−m|−1 ≤ h₀−2−2N < p². The prime p divides one of them iff p ≤ 2|k−m|−1, i.e. iff |k−m| ≥ (p+1)/2, and then at most to the first power. ∎

---

## 3. Theorem 1 (the refined bound at one prime)

For x > 0 put P(x) = [η_min x, (η₀−η_min)x], W(x) = {y ∈ P(x) : |y − η₀x/2| ≥ ½}, and for y ∈ P(x):

- ŝ(x,y) = #{j : η_jx ≤ y ≤ (η₀−η_j)x}  (the tower height),
- Ĉ(x,y) = 1 if there are j and an integer M ≠ 0 with η_jx − y ≤ M ≤ (η₀−η_j)x − y, else 0,
- D̂(x,y) = +∞ if Ĉ(x,y) = 1, and ⌊η₀x − y + ½⌋ + ⌊y + ½⌋ otherwise,
- **κ̂(x) = min( min_{y∈W(x)} φ(x,y),  min_{y∈P(x)} [ φ(x,y) + max(1, ŝ(x,y) − D̂(x,y)) ] )**, the first minimum omitted when W(x) = ∅.

(The minima are over closed sets of step functions with finitely many values, hence attained.)

**Theorem 1.** Let p be an odd prime with p² > 2h₀ and x = n/p. Then

  v_p(a₀) ≥ κ̂(x) − s  and  v_p(a_i) ≥ κ̂(x) − s  for every even i ≥ 2.

Hence the exponent of p in the common denominator of a₀, a₂, …, a_{s−1} is at most max(0, s − κ̂(n/p)).

*Proof.* Write y_k = k/p. For k ∈ K we have y_k ∈ P(x) (K = [η_min n, (η₀−η_min)n]); ŝ(x,y_k) = |J(k)| (the covering conditions divided by p); Ĉ(x,y_k) = [B_k^p ≠ ∅] (a nonzero multiple pM of p lies in [h_j*−k, h₀−1−h_j*−k] iff η_jx − y_k ≤ M ≤ (η₀−η_j)x − y_k); and when B_k^p = ∅, D_k(p) = |A_k^p| = N_k + [k≠m][p | m−k] = D̂(x,y_k) + ε'_k with ε'_k := [k≠m][p | m−k] (Lemma B(c),(b)). Put ε_k = Φ_p(k) − φ(x,y_k) ∈ {0,1} (Lemma B), so v_p(C_k) = φ(x,y_k) + ε_k − (s − s_k), and note ε_k − ε'_k = [k = m].

Each term of a₀ is (−1)^k a_{i,k} S_{i,k}; by Lemmas A and C its valuation is at least v_p(C_k) − min(s_k − i, D_k(p)) + v_p(S_{i,k}).

*Window poles* (|k−m| ≥ (p+1)/2, so k ≠ m and y_k ∈ W(x) because |y_k − η₀x/2| ≥ ½ + 1/(2p)): using min(s_k−i, D_k) ≤ s_k − i and v_p(S) ≥ −i, the valuation is ≥ v_p(C_k) − s_k = Φ_p(k) − s ≥ φ(x,y_k) − s ≥ min_{W(x)} φ − s.

*Other poles* (S_{i,k} is p-integral, i ≥ 1): the valuation is ≥ v_p(C_k) − min(s_k − 1, D_k(p)) = φ(x,y_k) − s + ε_k + [s_k − min(s_k−1, D_k)] = φ(x,y_k) − s + ε_k + max(1, s_k − D_k(p)).
 If B_k^p ≠ ∅: D_k = ∞, so this is ≥ φ − s + 1 = φ − s + max(1, ŝ − D̂) with D̂ = ∞.
 If B_k^p = ∅: D_k = D̂ + ε'_k, and max(1, a − ε') ≥ max(1, a) − ε' for ε' ∈ {0,1}, so the valuation is ≥ φ − s + (ε_k − ε'_k) + max(1, s_k − D̂) = φ − s + [k=m] + max(1, s_k − D̂). For k ≠ m, s_k = ŝ(x,y_k). For k = m, s_m = ŝ − 1 and 1 + max(1, ŝ−1−D̂) ≥ max(1, ŝ − D̂). In both cases the valuation is ≥ φ(x,y_k) − s + max(1, ŝ(x,y_k) − D̂(x,y_k)) ≥ min_{P(x)}[φ + max(1, ŝ − D̂)] − s.

So every term of a₀ has valuation ≥ κ̂(x) − s, hence so does a₀. For a_i with i ≥ 2: v_p(a_i) ≥ min_k v_p(a_{i,k}) and the same computation with i in place of 1 gives v_p(a_{i,k}) ≥ φ(x,y_k) − s + max(i, ŝ − D̂) ≥ κ̂(x) − s (at the centre use 1 + max(i, ŝ−1−D̂) ≥ max(i, ŝ−D̂)). ∎

**Remark 1 (evaluating κ̂).** φ(x,·), ŝ, Ĉ, D̂ are step functions of y whose breakpoints lie on the lines y = η₀x − M/2, y = M/2, y = η_jx + M, y = (η₀−η_j)x − M (M ∈ Z). Each term of φ is either right-continuous (⌊2y⌋−⌊y⌋, −⌊y−η_jx⌋) or left-continuous (⌊2(η₀x−y)⌋−⌊η₀x−y⌋, −⌊(η₀−η_j)x−y⌋) in y, so at a breakpoint of one family the value equals a one-sided limit; the value can be strictly below both one-sided limits only at a common breakpoint of a right- and a left-continuous term, which forces (η₀−η_j−η_l)x ∈ Z, (η₀−η_j)x ∈ ½Z or η₀x ∈ ½Z, i.e. p | n when p > 2η₀. The same holds for ŝ, Ĉ, D̂ (their breakpoints are consistent with φ's continuity directions). Hence for x ∉ Z the minima in κ̂ may be taken over the open intervals between consecutive breakpoints (midpoint evaluation, as multibeta.py does); the primes dividing n contribute at most s·log n to log D_n below and are negligible.

**Remark 2 (two special cases).** (a) If x ≥ 2/μ′ with μ′ = η₀ − 2η_min, the two components of W(x) have length (μ′x−1)/2 ≥ ½ each and their union modulo 1 is an interval of length μ′x − 1 ≥ 1, so min_{W(x)} φ = min_{y∈R} φ(x,y) =: φ₀(x) by the 1-periodicity of φ in y; the second minimum is ≥ φ₀(x) + 1. Hence **κ̂(x) = φ₀(x) for x ≥ 2/μ′** and the refinement lives entirely on the primes p > μ′n/2 = M′/2 (M′ = h₀−1−2N). (b) Dropping the cap (D̂ = ∞ everywhere) and the window (W = P) gives κ̂ = φ₀, i.e. Zudilin's bound v_p(a₀) ≥ −s + φ₀(n/p), which is his Lemma 5 with Φ_n = ∏ p^{φ₀(n/p)}.

---

## 4. Theorem 2 (denominators and their growth)

Let M = max(h₀−1−2N, h₁*) = μn, μ = max(η₀−2η_min, η₁), and define

  D_n = ∏_{p ≤ √(2h₀)} p^{s·⌊log_p M⌋} · ∏_{√(2h₀) < p ≤ M} p^{max(0, s − κ̂(n/p))}.

**Theorem 2.** (i) D_n a_i ∈ Z for i = 0, 2, 4, …, s−1. (ii) As n → ∞ along the admissible n,

  (log D_n)/n → H_ref := ∫_{1/μ}^{∞} max(0, s − κ̂(x)) dx/x².

If κ̂(x) ≤ s for all x ≥ 1/μ (which is checked for every η used below), then H_ref = s·μ − κ_ref with κ_ref = ∫_{1/μ}^{∞} κ̂(x) dx/x², and by Remark 2(a)

  κ_ref = κ_Z + ∫_{1/μ}^{2/μ′} (κ̂(x) − φ₀(x)) dx/x²,  κ_Z = ∫_{1/μ}^{∞} φ₀(x) dx/x²  (Zudilin's constant, log Φ_n/n → κ_Z).

*Proof.* (i) For p ≤ √(2h₀) (this includes p = 2) we use Zudilin's Lemma 4, first part: d_M^{s−i} a_{i,k} ∈ Z, together with d^i_{h₀−2N−2} S_{i,k} ∈ Z (the odd numbers in S_{i,k} are ≤ h₀−2N−2 ≤ M); hence d_M^{s} a₀ ∈ Z and d_M^{s−i} a_i ∈ Z, so the exponent of p in the denominators is at most s·v_p(d_M) = s⌊log_p M⌋. For √(2h₀) < p ≤ M use Theorem 1. For p > M: p exceeds h₀−1−2N ≥ |β| for all β ∈ B_k, so B_k^p = ∅; p exceeds L_j−1 and h₁*, so v_p(γ_n) = 0 and v_p(C_k) ≥ |A_k^p|; Lemma A gives v_p(a_{i,k}) ≥ |A_k^p| − min(s_k−i, |A_k^p|) ≥ 0; and S_{i,k} is p-integral (its odd numbers are < p). So a₀ and a_i are p-integral.
(ii) The small primes contribute at most s·log M·π(√(2h₀)) = O(√n). For the main range: κ̂ is a step function of x on [1/μ, ∞) — its breakpoints are the x-coordinates of intersections of the breakpoint lines of Remark 1 with each other and with the boundaries y = η_min x, y = (η₀−η_min)x, y = η₀x/2 ± ½, together with the jumps of the x-only floors; all are rationals with denominator ≤ 2η₀ — with finitely many pieces on every bounded interval, and equal to the 1-periodic φ₀ for x ≥ 2/μ′. On a piece [a,b) the prime number theorem gives Σ_{n/b<p≤n/a} log p = n(1/a − 1/b)(1+o(1)); the primes p ≤ n/X (x ≥ X) contribute at most s·Σ_{p≤n/X} log p ~ sn/X, which is o(n) after X → ∞; and the primes dividing n (where κ̂ is to be taken over closed sets, Remark 1) contribute O(log n). Summing the pieces gives the integral, and the formula for κ_ref follows from Remark 2(a). ∎

*Computation of κ_ref.* multibeta.py (`asym`) evaluates κ̂ on the Farey grid of order 2η₀ over [1/μ, 2/μ′] (midpoints of consecutive fractions, exact rational endpoints), integrates the resulting step function against dx/x² exactly, and adds the periodic tail Σ_{pieces} v·(ψ(X+b) − ψ(X+a)) with the digamma function (Σ_{m≥X}(1/(m+a) − 1/(m+b)) = ψ(X+b) − ψ(X+a)). The routine reproduces Zudilin's published constants: for his η = (31; 10⁵,11⁴,12⁴), κ_Z = 42.7664565011 (height 100.2335434989) and for his §2 construction (3; 1¹⁷), κ_Z = 0.9411124762; his φ₀ table for η = (31; …) is reproduced interval by interval.

---

## 5. Theorem 3 (five even β-values)

Take s = 11 and η = (33; 9, 10, 10, 11, 11, 12, 12, 12, 13, 13, 14); n even. The constraints hold (η_j < 16.5, Ση_j = 127 ≤ 165). Here η_min = η₁ = 9, μ′ = μ = 15, M = 15n, and the forms are r_n = a₀ + a₂β(2) + a₄β(4) + a₆β(6) + a₈β(8) + a₁₀β(10).

**Analytic side** (Zudilin, Lemma 2 and Lemma 3, cited): r_n > 0 and lim r_n^{1/n} = e^{−δ} with

  δ = −log[(4η₀)^{η₀}/(η₁^{2η₁}(η₀−2η₁)^{η₀−2η₁})] − log max_{t∈[0,1]^s} ∏_j t_j^{η_j}(1−t_j)^{η₀−2η_j}/(1+∏t_j)^{η₀} = 76.8117957305…,

the maximum being given by the unique zero x₀ = 0.001013931742… of x∏_j((η₀−η_j) − η_jx) − ∏_j(η_j − (η₀−η_j)x) in (0,1) (uniqueness checked numerically: one sign change on a 20000-point grid) and x_j = (η_j − (η₀−η_j)x₀)/((η₀−η_j) − η_jx₀) = 0.374128, 0.43396, 0.43396, 0.499239, 0.499239, 0.570745 (×3), 0.649414 (×2), 0.736378. (Exact measured closeness log|r_n|/n = −84.44, −81.57, −78.90, −78.48 at n = 2, 4, 12, 16, where the exact forms and the direct alternating sum agree to 242–1510 digits.)

**Arithmetic side** (Theorem 2). κ̂ ≤ 11 = s throughout (the maximum, 11, is attained on [1/15, 1/14)), so H_ref = 11·15 − κ_ref. The steps of φ₀ and of κ̂ on the window [1/15, 2/15):

| x-piece | φ₀ (Zudilin) | κ̂ (refined) | gain | ∫ dx/x² | rebate |
|---|---|---|---|---|---|
| [1/15, 1/14) | 2 | 11 | 9 | 15 − 14 = 1 | 9 |
| [1/14, 1/13) | 2 | 10 | 8 | 1 | 8 |
| [1/13, 1/12) | 4 | 10 | 6 | 1 | 6 |
| [1/12, 1/11) | 4 | 7 | 3 | 1 | 3 |
| [1/11, 1/10) | 7 | 7 | 0 | | 0 |
| [1/10, 1/9) | 5 | 5 | 0 | | 0 |
| [1/9, 4/33) | 3 | 3 | 0 | | 0 |
| [4/33, 1/8) | 4 | 4 | 0 | | 0 |
| [1/8, 2/15) | 3 | 3 | 0 | | 0 |

So κ_ref − κ_Z = 26 exactly, κ_Z = 62.8744847732…, κ_ref = 88.8744847732…, and

  H_ref = 165 − 88.8744847732… = 76.1255152268…  (Zudilin's bound would give 165 − 62.8745 = 102.1255).

Independent check: brute-force integration of κ̂ and φ₀ on 400000 points of the variable u = 1/x gives 88.8739 and 62.8739 (difference 26.0000).

**Conclusion.** D_n r_n ∈ Z + Zβ(2) + Zβ(4) + Zβ(6) + Zβ(8) + Zβ(10) is positive and (log D_n)/n − log(1/r_n)/n → 76.1255 − 76.8118 = −0.6863 < 0, so D_n r_n → 0. If β(2), β(4), β(6), β(8), β(10) were all rational with common denominator q, then q·D_n r_n would be a positive integer tending to 0. Hence **at least one of β(2), β(4), β(6), β(8), β(10) is irrational.**

Corroborating directions (same computation): (36; 10,10,11,11,12,12,13,13,14,14,15): δ = 86.6348, H_ref = 176 − 92.2457 = 83.7543, net −2.88 per n; (40; 11,12,12,13,13,14,14,15,15,16,17): δ = 94.7080, H_ref = 198 − 108.9555 = 89.0445, net −5.66; (54; 16,17,17,18,19,19,20,21,21,22,23): net −8.56. For Zudilin's own η = (31; 10⁵,11⁴,12⁴), s = 13: κ_ref − κ_Z = 11/3, H_ref = 96.5669 against δ = 100.7397 (his margin 0.506 becomes 4.173). Under Zudilin's bound alone none of the s = 11 directions closes (their Zudilin-heights exceed δ by 17 to 26 per n): the gain comes entirely from spreading the tower so that the poles feeding the top primes have one or two floors.

---

## 6. Exact verification of Theorem 1 (multibeta.py `forms`)

For each η and n the script computes every a_{i,k} exactly (log-series of B_k in fmpq), checks R(t) = Σ a_{i,k}(t+k+½)^{−i} at random rational t, checks the odd a_i vanish and that Σ a_iβ(i) + a₀ equals the direct alternating sum to hundreds of digits, factors the denominators, and prints, for every odd prime √(2h₀) < p ≤ M, the exact exponent of p in den(a₀) and den(a₂) next to Zudilin's bound s − φ₀(n/p), the finite-n refined bound (using the exact v_p(C_k), D_k(p), s_k and window) and the asymptotic s − κ̂(n/p).

| η | n | primes tested | exponent > refined bound | exponent = refined bound | exponent < refined bound (mechanism II) |
|---|---|---|---|---|---|
| (31; 10⁵,11⁴,12⁴), s=13 | 6 | 10 | 0 | 10 | 0 |
| same | 8 | 15 | 0 | 14 | 1 (p = 23) |
| (5; 2⁵), s=5 | 20 | 6 | 0 | 6 | 0 |
| (33; 9,10,10,11,11,12,12,12,13,13,14), s=11 | 4 | 8 | 0 | 8 | 0 |
| same | 8 | 20 | 0 | 19 | 1 (p = 41) |
| same | 12 | 31 | 0 | 29 | 2 (p = 37, 43) |
| same | 16 | 40 | 0 | 40 | 0 |
| (36; 10,10,11,11,12,12,13,13,14,14,15), s=11 | 3, 5 | 8, 14 | 0 | all | 0 |
| (3; 1⁵), s=5 | 60, 100 | 9, 16 | 0 | 7, 14 | 2, 2 |
| (3; 1⁷), s=7 | 80 | 14 | 0 | 8 | 6 |

The bound was never exceeded, and it is attained at most primes; the occasional extra unit ("mechanism II", an alternating window sum of the a_{i,k} vanishing mod p) is not used here and could only lower denominators further. At the finite-n level the exact exponent at, e.g., η = (31; …), n = 8, p = 83 is 3 against Zudilin's 7: the poles feeding 1/83 are k = 80, 81, 82 of multiplicity 5 with v_p(C_k) = 2 (giving 5 − 2 = 3), while the multiplicity-13 poles have D_k = 3 (numerator factors −83, 83, 249 and no denominator factor divisible by 83), so their a_{i,k} carry at most 83^{−3}.

---

## 7. Remarks

3. **What the refinement is.** Zudilin charges every pole s − i derivatives and lets every pole feed the partial sums; the refinement charges each pole only its own floors (Lemma B), caps the derivative cost by the count of p-divisible numerator factors when no denominator factor is p-divisible (Lemma A), and lets only the poles beyond distance p/2 from the centre feed p into the partial sums (Lemma C). All three are exact bookkeeping of the same factors his proof uses.
4. **Recovering his bound.** With W(x) = P(x) and D̂ = ∞, κ̂ = φ₀ and Theorem 1 is v_p(a₀) ≥ −s + φ₀(n/p), i.e. Φ_n^{-1}d_M^s a₀ ∈ Z at the primes p² > 2h₀.
5. **Where it acts.** Only for x = n/p < 2/μ′, i.e. p > M′/2: at those primes the window W(x) misses part of a period and the poles inside it are the edge poles. For a tower with steps (η_j distinct), the edge poles have one floor and the top primes cost nothing: for η = (33; …) the primes in (14n, 15n] do not divide the denominators at all (κ̂ = s), which is what the exact runs show (n = 8: p = 97…109 have exponent 1 at that finite n, asymptotically 0).
6. **Lenses.** Hyperoperations: the tower s_k, floor by floor, is the whole content of Lemmas A and B. Hybrid base: φ, D̂ and N_k are base-p digit and carry counts of the pair (n/p, k/p); κ_ref is an integral of a step function with rational breakpoints, so all rebates are rational (11/3, 26, 20, 86/3). Palindrome: R(−t−h₀) = R(t) removes the odd β's; its local shadow mod p is mechanism II.
7. **What is not proved here.** Mechanism II; optimality of the η's (found by coordinate descent); anything about β(2) alone (s = 3, 5, 7 remain far from closing in this family: best per unit η₀ found +0.67 for s = 7, +0.23 for s = 9).
8. **Prior art checked (2026-09-28).** Rivoal–Zudilin 2003 (Math. Ann. 326), §5 and §7: uniform derivative cost (−λ per derivative), no window (every pole feeds the partial sums), and the minimum of their carry count ϖ₀ over all y; their Lemma 6 ties the partial-sum lcm ranges to the pole class, which is the only tower-aware element there. Their per-factor Lemmas 10–11 are the ingredients that Lemma B above re-derives as an exact identity. Zudilin 2018 Lemma 4/5 is the same structure. Fischler 2019 is a Padé/Siegel method. None contains Lemmas A–C or Theorem 1.

---

## 8. Theorem 4: the per-order mechanism II is the palindrome reduced mod p

Empirically (results log, "The per-order rule"), a₀ loses one more power of p than Theorem 1 predicts whenever the leading poles form blocks centred on poles c ≡ m (mod p). This section proves that statement. Throughout p is an odd prime with p² > 2h₀ and p ∤ n, so p > 2η₀ as well for all n considered.

**Definitions.** For a pole k write D(k) = m − k. The *factor families* at k are: the numerator family {α_i = i − k − ½ : 1 ≤ i ≤ h₀−1}, the centre factor 2D(k) (absent at k = m), and for each j the denominator family {β_{j,i} = h_j* + i − k : 0 ≤ i ≤ L_j − 1}. The *digit multiset* of a family at k is the multiset of the nonzero p-divisible members divided by p, reduced mod p (for the numerator, the odd integers q with pq = 2α; for a denominator family, the integers b with pb = β ≠ 0). A *digit cell* is a maximal interval of poles on which every family's digit multiset is constant. A *fixed-point pole* is a pole c ≡ m (mod p) with c ≠ m; write c = m − jp.

**Lemma D (reflection of the families).** Let c = m − jp and k* = 2c − k. Then, as multisets, {β_{j,·}(k*)} = {2jp − β : β ∈ {β_{j,·}(k)}} for every j, {α(k*)} = {2jp − α : α ∈ {α(k)}}, and 2D(k*) = 4jp − 2D(k).

*Proof.* β_{j,i}(k*) = h_j* + i − k*; substituting i ↦ L_j − 1 − i gives h_j* + L_j − 1 − i − k* = (h₀ − 1 − h_j*) − i − k* = −β_{j,i}(k) + (h₀ − 1 − k − k*) = −β_{j,i}(k) + 2(m − c) = 2jp − β_{j,i}(k). Likewise α_{h₀−i}(k*) = h₀ − i − k* − ½ = −α_i(k) + (h₀ − 1 − k − k*) = 2jp − α_i(k). The centre factor is linear in k. ∎

**Lemma E (digit cells around a fixed-point pole are symmetric, and their parity is even).** Let Z be a digit cell containing a fixed-point pole c = m − jp, and let k ∈ Z with k* = 2c − k ∈ Z, k ≠ c. Then:
(i) for each denominator family, the multiset of *all* multiples of p among its members at k (including the pole's own 0 if the family covers k) is a set of consecutive integers times p, symmetric about jp; hence its nonzero digit multiset D_j has even cardinality if the family covers k and odd cardinality otherwise;
(ii) the numerator's odd digits q are symmetric about 2j; hence their number |A^p| is even;
(iii) the number of p-unit members of each family is the same at k and at k*;
(iv) v_p(a₀-leading data): v_p(C_k) = v_p(C_{k*}) and the leading Taylor profile (Lemma A) is the same at k and k*.

*Proof.* By Lemma D the multiset at k* is the image of the multiset at k under β ↦ 2jp − β; by constancy on Z the two multisets coincide, so the multiset at k is symmetric under β ↦ 2jp − β. The multiples of p in an interval form consecutive multiples; a symmetric set of consecutive multiples of p about jp is {(j−w)p, …, (j+w)p}, which contains jp and has odd size 2w+1. If the family covers k, one of these is the pole's own zero (0 = (j−w)p requires w = j; in any case 0 lies in the range, and 0 is one of the symmetric multiples), so the nonzero digits number 2w, even; if not, all 2w+1 are nonzero, odd. For the numerator the members are half-integers α = o/2 with o odd; 2jp − α = (4jp − o)/2, so the odd multiples q of p (with pq = o) are symmetric about 2j, an even number: odd integers symmetric about an even integer come in pairs, so |A^p| is even. (iii) is clear from constancy. (iv): C_k and C_{k*} have the same p-divisible multisets, hence the same valuation by Lemma B's count (the centre factor is a unit at both, as p ∤ D(k) for k ≢ m), and by Lemma A the leading coefficient of [u^r]B_k at valuation v_p(C_k) − r is unit(C_k)·[w^r]F(w) with F(w) = ∏_{α∈A^p}(1 + w/a_α)/∏_{β∈B^p}(1 + w/b_β), where a_α = α/p, b_β = β/p are the digits; F depends on the digit multisets only. ∎

**Lemma F (the sign).** In the situation of Lemma E, unit(C_{k*}) ≡ −unit(C_k) (mod p).

*Proof.* C_k = 2γ_n · 2D(k) · ∏_i α_i(k) / ∏_{j,i} β_{j,i}(k) (nonzero members). Compare family by family using Lemma D. Centre factor: 2D(k*) = 4jp − 2D(k) ≡ −2D(k), ratio −1. For any other family, its members at k* are 2jp − β over the members β at k. A p-unit member satisfies 2jp − β ≡ −β, contributing a factor −1 to the ratio of units; a p-divisible member β = pb has 2jp − β = p(2j − b), whose unit part is 2j − b, and by Lemma E the multiset {2j − b} equals the multiset {b}, so the p-divisible members contribute a total ratio of 1. Hence the ratio of units of a family is (−1)^{N} with N its number of p-unit members. Summing over families, the total ratio is −(−1)^{ΣN}, and ΣN ≡ (h₀ − 1 − |A^p|) + Σ_j (L_j − [j covers k] − |D_j|) (mod 2). Here h₀ − 1 = η₀n is even, each L_j = (η₀ − 2η_j)n + 1 is odd, Σ_j [j covers k] = s_k, and by Lemma E(i),(ii) |A^p| is even while |D_j| is even for the s_k covering families and odd for the s − s_k others. So ΣN ≡ 0 + s + s_k + (s − s_k) ≡ 2s ≡ 0 (mod 2), and the total ratio is −1. ∎

(Numerical confirmation: at (3;1⁵), n = 80, p = 31, block centre c = 89, pair (85, 93): the centre factor has ratio −1 and all six other families ratio +1; at the winner η, n = 20, p = 61, c = 208, pair (205, 211): centre −1, twelve families +1. In the non-cancelling pair (150, 156) around c = 153 at n = 16, p = 37, the three η_j = 12 families have ratio 32: the pair straddles a digit boundary and Lemma E does not apply; the inner pair (152, 154) of the same centre is antisymmetric.)

**Lemma G (the partial-sum multipliers agree).** For k and k* = 2c − k in the lower window with |D(k) − jp| < p/2, the partial sums S_{i,k} and S_{i,k*} contain exactly the odd multiples p, 3p, …, (2j−1)p of p, so their p-singular parts are p^{−i}·c_i with the same c_i = Σ_{q odd ≤ 2j−1} (−1)^{(qp+1)/2}(−2/q)^i.

*Proof.* The odd numbers in S_{i,k} are those < 2D(k); the number of odd multiples of p below 2D is ⌊(2D − 1 + p)/(2p)⌋, which equals j for jp − p/2 < D < jp + p/2, and D(k*) = 2jp − D(k) lies in the same range. ∎

**Theorem 4.** Let p be an odd prime with p² > 2h₀ and p ∤ n. Let Z be a digit cell containing a fixed-point pole c = m − jp, lying in the lower window (Z ⊆ [N, (h₀−2−p)/2]), and symmetric about c. Then for every k ∈ Z ∖ {c} and every i, the terms (−1)^k a_{i,k} S_{i,k} and (−1)^{2c−k} a_{i,2c−k} S_{i,2c−k} of a₀ have leading parts (at valuation v_p(C_k) − s_k) that are exact negatives, so their sum has valuation ≥ v_p(C_k) − s_k + 1. Consequently, if the set of window poles k with Φ_p(k) = min_{window} Φ_p is a disjoint union of such cells Z (minus their centres, which are never leading since Φ_p(c) = φ + 1), then

  v_p(a₀) ≥ κ̂_W(x) − s + 1,  κ̂_W(x) := min_{window} Φ_p,

one unit better than Theorem 1 at that prime.

*Proof.* k and k* = 2c − k have the same parity. By Lemma E(iv) their coefficients a_{i,k}, a_{i,k*} have equal valuation and equal normalized leading coefficients h_{s_k − i}, and by Lemma F unit(C_{k*}) = −unit(C_k), so a_{i,k*} ≡ −a_{i,k} at leading order; by Lemma G the p^{−i} multipliers coincide; hence the leading parts cancel. All other terms of a₀ have valuation ≥ min_window Φ_p − s + 1 by the proof of Theorem 1 (non-leading window poles by one unit; centre poles by Lemma B's ε = 1; the band poles by the non-PS estimate). ∎

**Remark 9 (relation to the computable rule).** For p ∤ n the breakpoints of different families in y do not coincide (a coincidence forces (η_j ± η_l)x ∈ Z or η₀x ∈ ½Z, i.e. p | n for p > 2η₀), so the digit cells are exactly the cells on which φ(n/p, ·) is constant, and the hypothesis of Theorem 4 is the statement that the argmin set of φ over the window, reduced mod 1, is a single arc containing the fixed point y₀ = η₀x/2 and not cut by the window edge. This is the rule scored in the results log (precision 1.00 and recall 1.00 on the 36 events of the winner direction). Its asymptotic worth is ∫[rule holds] dx/x²: 0.717 per n for η = (33; 9,…,14), so the margin of Theorem 3 becomes 0.686 + 0.717 = 1.403 per n; 0.918 for Zudilin's η; 0.2545 for the simple family, where the rule reads "⌊n/p⌋ even and {n/p} ∉ [1/3, 1/2)".

**Remark 10 (what is not covered).** The exact centre m (j = 0) is different: there the digit multisets are symmetric under β ↦ −β, F(−w) = F(w), the odd Taylor coefficients vanish, and the pairing sign is (−1)^i, the global palindrome. The cross-order events (poles whose blocks vanish individually at the midpoint to the pole p away, and the cell-sum events) are not covered by Theorem 4.

---

## 9. Theorem 5: the leading residue of a₀, and the three shapes of mechanism II

The lemmas above give more than a bound: they give the leading p-adic term of a₀ explicitly. Throughout, p is odd with p² > 2h₀, p ∤ n, and we assume the *window-dominated regime*: the minimum φ_min := min_{k ∈ window} Φ_p(k) satisfies φ_min < min_{all k} [Φ_p(k) + max(1, s_k − D_k(p))] (so the leading order of a₀ in Theorem 1 comes from the partial sums; for x ≥ 2/μ′ this is automatic, and it holds at every prime of the data sets below).

**Lemma A′ (exact leading coefficient).** With the notation of Lemma A, write α = p·a for α ∈ A_k^p and β = p·b for β ∈ B_k^p (a, b are p-adic units, a half-integers). Then for 0 ≤ r ≤ s_k − 1,
  [u^r] B_k(u)/C_k = p^{−r} h_r^{(k)} + (terms of valuation ≥ −r + 1),  where  Σ_r h_r^{(k)} w^r = F_k(w) := ∏_{α∈A_k^p}(1 + w/a) / ∏_{β∈B_k^p}(1 + w/b).
*Proof.* In the expansion of Lemma A a term has valuation exactly −r iff every one of its r factors is p-divisible; those terms are exactly the expansion of ∏_{A^p}(1 + u/(pa))/∏_{B^p}(1 + u/(pb)) = F_k(u/p). ∎

**Lemma C′ (exact singular part of the partial sums).** For k in the lower window let J_k = ⌊(2(m−k) − 1 + p)/(2p)⌋ be the number of odd multiples of p below 2(m−k). Then the p-singular part of S_{i,k} is Σ_{q odd ≤ 2J_k − 1} (−1)^{(qp+1)/2} (−2/(qp))^i (the terms ℓ = −(qp+1)/2), all other terms being p-integral. ∎

**Theorem 5.** Let L be the set of lower-window poles with Φ_p(k) = φ_min. Then
  a₀ ≡ 2 p^{φ_min − s} · Σ_{k∈L} (−1)^k unit(C_k) · T_k   (mod p^{φ_min − s + 1}),
  T_k := Σ_{q odd ≤ 2J_k − 1} (−1)^{(qp+1)/2} (−2/q)^{s_k} Σ_{r=0}^{s_k−1} h_r^{(k)} (−q/2)^r.
Consequently the exponent of p in den(a₀) is s − φ_min unless Σ_{k∈L} (−1)^k unit(C_k) T_k ≡ 0 (mod p), in which case it is at most s − φ_min − 1.
*Proof.* a₀ = 2Σ_{k<m}Σ_i (−1)^k a_{i,k} S_{i,k} (the two halves of a₀ are equal by the palindrome). For k ∈ L the terms of valuation φ_min − s come from Lemma A′ and Lemma C′: Σ_i C_k h^{(k)}_{s_k−i} p^{−(s_k−i)} (−2/(qp))^i = C_k p^{−s_k} (−2/q)^{s_k} Σ_r h_r^{(k)} (−q/2)^r, and v_p(C_k) − s_k = Φ_p(k) − s = φ_min − s. All other terms have valuation ≥ φ_min − s + 1: non-leading window poles by one unit, the pole c ≡ m by Lemma B (ε_c = 1), band poles by the window-dominated hypothesis. ∎

T_k depends only on the digit multisets at k, on s_k and on J_k — it is constant on a digit cell — and it is a rational number with p-unit denominator (its denominators are products of digits). Three ways for the sum to vanish, which are the three observed shapes of mechanism II:
- **Per-order (Theorem 4):** L is a union of cells symmetric about fixed-point poles; T is constant on each cell and unit(C) is antisymmetric, so the sum vanishes pairwise, at every order i separately.
- **Per-pole:** p | numerator(T_k) for every k ∈ L; every pole vanishes on its own. For a fixed digit configuration T is a fixed rational, so the carrier primes are its large prime factors. In the simple family at x ∈ (1, 2) (J = 1, the cell at the window edge N): s = 5: A = {−1,1,3}, B = {1}⁵, T = 145/24 = 5·29/24; A = {−3,−1,1,3,5}: 781/120 = 11·71/120. s = 7: A = {−3,−1,1,3,5}: 15211/720 = 7·41·53/720. s = 9: A = {−1,1,3}: 8723/128 = 11·13·61/128. s = 11: A = {−3,−1,1,3,5}: 66781/256 = 11·13·467/256; A = {−1,1,3,5}: 187187/960 = 7·11²·13·17/960. The carriers observed in the sweeps — 29, 71; 41, 53; 61; 13, 17 — are exactly these prime factors, and primes sharing the configuration without dividing (31, 43, 47, 73) never fire. (Whenever the digit q = 1 is present, F(−1/2) = 0 exactly, so T is a truncation tail.)
- **Cell-sum:** the weighted sum over several cells with different T vanishes; an accident of the weights.

**Verification (cross_theory.py).** For every (n, p) pair of the seven sweep files — 818 + 594 + 1673 + 1070 (simple families s = 5, 7, 9, 11), 883 (η = (33; 9,…,14), n ≤ 40), 421 (η = (36; …)), 204 (Zudilin's η) — the residue Σ_{k∈L}(−1)^k unit(C_k) T_k was computed from the digits alone and compared with the exact exponent: agreement in 5663 of 5663 pairs (100%). Of the cross-order events, all leading cells have T ≡ 0 in 24/24 (s=5), 23/23 (s=7), 14/24 (s=9), 7/11 (s=11), 1/3, 1/6, 2/3 (the three η's); the rest are cell-sum.

**Remark 11.** Theorem 5 makes the exponent of every prime in den(a₀) computable from base-p digits alone, without the exact forms: Theorem 1 gives s − φ_min, and the residue decides the extra unit. The asymptotic weight of the per-order shape is ∫[Theorem 4's hypothesis] dx/x² (Remark 9); the per-pole and cell-sum shapes are supported on the primes dividing specific integers and carry no density.

---

## 10. Theorem 6: the asymptotic weight of Theorem 4, and the digit-only exponent

**Digit-only computation.** Theorems 1 and 5 together determine the exponent of every odd prime p with p² > 2h₀, p ∤ n in den(a₀) from base-p digits alone: Φ_p(k) by Lemma B's count, the leading set L, the units unit(C_k) mod p by a boundary recurrence along the poles (moving k to k+1 changes each factor family by one member at each end, so unit(C_{k+1})/unit(C_k) is a product of s+2 small ratios; the starting value at k = N uses the p-adic factorial units (Anton–Wilson: unit(M!) ≡ Π_i (m_i!)(−1)^{⌊M/p⌋+…} over the base-p digits m_i of M)), the digit profiles F and their Taylor coefficients mod p, the multipliers J_k, and the residue Σ_{k∈L}(−1)^k unit(C_k) T_k. Implemented as `multibeta.py digitexp ETA n p`; validated against the exact exponents on all sweep pairs with p ∤ n: 802 + 581 + 1647 + 883 + 421 + 204 = 4538 agreements out of 4538.

**Theorem 6 (rigorous per-order rebate).** Let η be admissible, μ′ = η₀ − 2η_min, and for x ∈ [1/μ, 2/μ′) let G(x) be the indicator of: "the set where φ(x, ·) attains its minimum over the lower window arc [η_min x, η₀x/2 − ½], reduced mod 1, is a single arc containing y₀ = η₀x/2 (mod 1) in its interior, and not containing the window endpoint η₀x/2 − ½ (mod 1)" — a step function of x with rational breakpoints of denominator ≤ 2η₀. Then the denominators D′_n := D_n ∏_{√(2h₀)<p≤M, p∤n, G(n/p)=1} p^{−1} still satisfy D′_n a_i ∈ Z for all i, and

  lim (log D′_n)/n = H_ref − ∫_{1/μ}^{2/μ′} G(x) dx/x².

*Proof.* Fix such a prime p with G(n/p) = 1 (x = n/p not a breakpoint of G, which excludes finitely many n per prime and hence a set of primes of total weight o(n)). By Remark 1 (p ∤ n) the cells of φ(n/p, ·) are the digit cells; the argmin arc is a union of digit cells symmetric about the fixed point y₀ (φ is symmetric under y ↦ η₀x − y), i.e. of cells symmetric about the poles c ≡ m (mod p) inside the window, and it stays away from the window edge by a fixed positive distance δ(x) > 1/(2p) for p large, so the finite window contains each such cell entirely. The leading set of a₀ at p is therefore a union of cells satisfying the hypothesis of Theorem 4, and v_p(a₀) ≥ κ̂(x) − s + 1. For the a_i, i ≥ 2, Theorem 1 already gives v_p(a_i) ≥ κ̂ − s + 1 (the non-PS term with i ≥ 2 exceeds the PS term by at least one). Hence one factor p can be removed from D_n at every such prime, and the prime number theorem over the pieces of G gives the integral, exactly as in Theorem 2. ∎

For η = (33; 9,10,10,11,11,12,12,12,13,13,14): ∫G dx/x² = 0.7172 (numerically: G = 1 on (0.308, 0.333), (0.445, 0.455), (0.462, 0.467), (0.485, 0.5), (0.571, 0.576), (0.615, 0.625), (0.667, 0.697), (0.923, 0.939), (0.95, 1.0) and on the periodic continuation — note G is then given by φ₀'s argmin structure for x ≥ 2/μ′, where the window arc covers a period and the condition reads "the argmin set of φ₀ is a single arc around y₀"), so the rigorous height of Theorem 3 becomes 76.1255 − 0.7172 = 75.4083 per n against the decay 76.8118: **margin 1.4035 per n**. For Zudilin's own direction the rebate is 0.9178 (margin 5.09); for the simple family 0.2545.
