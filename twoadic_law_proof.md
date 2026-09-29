# The 2-adic content of the lattice Catalan Hankel machines: what is proved and what is verified
(2026-09-28, real date.) Companion to `twoadic_law.py` and the results-log sections "The exact 2-adic structure…", "Both G machines are continuous dual Hahn weights", "Proof route, first stage", "THE PARTNER LAW, PINNED DOWN" and "Unified law and Catalan".

## 0. Status
- RIGOUR PASS (2026-09-29, prompted by the user's review of Theorem 7). Everything is now over ℚ₂, with integral monic bases.
  - The review: the old Theorem 7 wrote the CDH Gram matrix as D^{1/2}·π(J)·D^{1/2} with the orthonormal Jacobi matrix and L = W·√(h_q/h_q′). That needs √h_q and √λ_n, which do not lie in ℚ₂, and it did not show that L lies in the required integral ring.
  - Section 12, Theorem 7:
    - uses the monic Jacobi matrix 𝒥 (φ₀(f·P_pP_q) = h_p·f(𝒥)_pq, Lemma 12.1);
    - rescales by the rational diagonal diag(4^{−n}) (Catalan) or diag(2^{−n}) (twisted) to an integral matrix divisible by 4 (resp. 2);
    - runs Gaussian elimination in ℤ₂;
    - conjugates back, which gives an explicitly integral L with v₂(L_pq) ≥ 3(p−q) + 2 (resp. + 1) and v₂(Λ_p) = v₂(h_p).
  - The same device appeared in section 9. Its criterion is now Lemma 9.1, proved by determinant expansion, and its smoothness step uses balanced valuations v₂(M_ac) + ½(ν(a) − ν(c)).
  - Theorem 4 now states the explicit integer recurrences of the Krawtchouk polynomials (and corrects a slip in the twisted norm line). Theorem 6(a) spells out its analytic steps. Theorem 8 proves the Pascal reduction and det[C(a+i, j)] = 1.
  - New checks: proof_checks/th_thm7_rigorous.py, th_lemmaA_balanced.py, th_kraw_recurrence.py.
  - The statements are unchanged; only the proofs are tightened.
- UPDATE (2026-09-28, third stage, section 12): the MAIN THEOREM IS PROVED for both machines and all K. Lemmas A, B and C are all proved.
  - A(X) = A(C) + (X − C)B is a pencil of two Hankel forms:
    - A(C) is the CDH moment form at the 2-adic Catalan constant, with Smith form {h_q};
    - B is the Krawtchouk form, with Smith form equal to the band.
  - Mod 2 the two bases differ by the Pascal matrix, so Cauchy–Binet gives the coefficient law, the content and det laws (tie rule included) and the Smith form.
  - The key new facts:
    - the X-coefficient of the machine is Krawtchouk;
    - at X = C, the pole values are the 2-adic Stieltjes transform of the moments (Theorem 6).
  - Theorem 6(b) (twisted) is now written out in full, continuity step included. No gaps remain. (The proof is unrefereed, like everything else here.)
- UPDATE (2026-09-28, second stage, sections 8–11): the X-coefficient matrix B is exactly a Krawtchouk (binomial) Hankel matrix, and there is an exact split A(0) = −C·B + A(C), where C is a 2-adic Catalan constant with v₂(C) = −1.
  - Lemma A is PROVED for Catalan (the whole pole band, every m < K) and for the twisted first half (m < ⌊K/2⌋). The method is the Krawtchouk mirror plus Mahler smoothness.
  - Lemma C (the content +1 at ties) is PROVED from a coefficient law, and that law is verified.
  - Lemma B is reduced to two explicit statements about the moment corner in the Krawtchouk basis (section 11); these are verified but not proved.
- PROVED (sections 1–5):
  - both weights are continuous dual Hahn;
  - the moment-Hankel Smith form in the v-scale;
  - Mahler smoothness of the pole data;
  - the exact algebraic reductions (content bound, Uvarov form, divided-difference factorisation, the K-independent matrix T, pivots = Newton–Padé residuals, the u↔v invariance of the tropical function);
  - the closed form of the divided differences at r = 0, and its reduction to a single congruence.
- VERIFIED, NOT PROVED (section 6): the Main Theorem, i.e. the Smith form is the union of the pole band and the moment pivots, split by a tropical minimum, with the tie rule.
  - Exact Smith-form equality at 26 values of K (13 per machine), including 4 ties.
  - The content law at 134 values of K per machine: all K = 3…128 and K = 160, 180, 192, 224, 256, 320, 360, 384.
  - No failures.
- OPEN: Lemmas A (pole band), B (merge) and C (ties) in section 7. A proof of these three gives the Main Theorem.

## 1. Setting
- Two machines. In both, A_K is the K×K Hankel matrix φ₀(x^{i+j}/D(x)) with the pole value G replaced by X = 0.
  - Twisted: variable u = t², weight w_T = πt·sinh(πt)/cosh²(πt), poles u = −a², a = 0…K−1.
  - Catalan: variable v = 4t², weight w_C = πt²·cosh(πt)/sinh²(πt), poles v = −b², b = 1, 3, …, 2K−1.
- Notation: s₂ is the binary digit sum, v₂(n!) = n − s₂(n) (Legendre), and D(x) = Σ_{j<x} s₂(j).
- The content exponent is e₂(K) = −min_j v₂(coefficient of X^j in det(A + XB)).

## 2. Theorem 1 (continuous dual Hahn). PROVED.
1. The weights are
   - w_T(t) = |Γ(½+it)|⁶ / (4π|Γ(2it)|²);
   - w_C(t) = |Γ(1+it)² Γ(it)|² / (4π|Γ(2it)|²).
   These are continuous dual Hahn weights with (a, b, c) = (½, ½, ½) and (1, 1, 0).
2. The monic orthogonal polynomials in u = t² satisfy u·p_n = p_{n+1} + b_n p_n + λ_n p_{n−1}, with
   - twisted: b_n = 2n² + 2n + ¾, λ_n = n⁴;
   - Catalan: b_n = (2n+1)(n+1), λ_n = n³(n+1).
3. The norms (h₀ = μ₀ = ½) are
   - twisted: h_n = (n!)⁴/2;
   - Catalan: h_n = (n!)³(n+1)!/2.

Proof.
- Weights. Use |Γ(½+it)|² = π/cosh πt, |Γ(1+it)|² = πt/sinh πt, |Γ(it)|² = π/(t·sinh πt) and |Γ(2it)|² = π/(2t·sinh 2πt).
  - |Γ(½+it)|⁶/|Γ(2it)|² = (π³/cosh³)·(2t·2sinh·cosh/π) = 4π²t·sinh/cosh².
  - |Γ(1+it)|⁴|Γ(it)|²/|Γ(2it)|² = 4π²t²·cosh/sinh².
- Recurrence. For S_n (Koekoek–Lesky–Swarttouw 9.3) the three-term relation reads −(a²+x²)S̃_n = A_n S̃_{n+1} − (A_n+C_n)S̃_n + C_n S̃_{n−1}, with S̃_n = S_n/((a+b)_n(a+c)_n), A_n = (n+a+b)(n+a+c), C_n = n(n+b+c−1).
  - Put p_n = (−1)ⁿS_n (monic in u = x²).
  - Using (a+b)_{n+1}(a+c)_{n+1}/((a+b)_n(a+c)_n) = A_n, this gives u·p_n = p_{n+1} + (A_n + C_n − a²)p_n + A_{n−1}C_n p_{n−1}.
  - Substituting the parameters gives the stated b_n and λ_n.
- Norms. h_n = h₀·∏_{k≤n} λ_k with h₀ = μ₀ = ½; the Catalan parameter c = 0 is reached by continuity of the moments.
- Independent check: `twoadic_law.py cdh` confirms orthogonality and the norms against the exact moments for n < 30 (0 violations). ∎

## 3. Theorem 2 (moment Hankel Smith form in the v-scale). PROVED.
In v = 4t² the n×n moment Hankel matrix [μ_v(i+j)] has 2-adic elementary divisors exactly {v₂(h^v_q) : q < n}:
- twisted: h^v_q = 16^q (q!)⁴/2, so v₂(h^v_q) = 4q + 4v₂(q!) − 1 = 8q − 4s₂(q) − 1;
- Catalan: h^v_q = 16^q (q!)³(q+1)!/2, so v₂(h^v_q) = 4q + 3v₂(q!) + v₂((q+1)!) − 1.

Proof.
- The v-scale recurrence coefficients are 4b_n and 16λ_n: (8n²+8n+3, 16n⁴) for twisted and (4(2n+1)(n+1), 16n³(n+1)) for Catalan. Both are integral, so the monic orthogonal polynomials p^v_q have integer coefficients.
- The coefficient matrix C (rows p^v_q in the monomial basis) is unit lower triangular and integral, hence in GL_n(ℤ₂).
- Hence [μ_v(i+j)] = C⁻¹ diag(h^v_q) C⁻ᵀ, and the Smith form is the diagonal. ∎
- (Checked for n = 8, 16, 24.)

## 4. Lemma 3 (Mahler smoothness of the pole data). PROVED.
F(a) = −4(−1)^a β_a with β_a = Σ_{k<a} (−1)^k/(2k+1)². For every r ≥ 0 and n ≥ 1, v₂(ΔⁿF(r)) = n + 1.

Proof.
- Put G(a) = (−1)^a β_a. Then G(a) + G(a+1) = −g(a) with g(a) = (2a+1)⁻², i.e. (Δ+2)G = −g.
- For G_r(a) = G(a+r), g_r(a) = g(a+r), the Mahler coefficients M_n = ΔⁿG_r(0) satisfy M_{n+1} = −2M_n − Δⁿg_r(0), with M₀ = G(r). Hence M_n = (−2)ⁿG(r) − Σ_{k<n} (−2)^{n−1−k} Δ^k g_r(0).
- Expanding g_r(a) = (2r+1)⁻²(1 + 2a/(2r+1))⁻² 2-adically and using Σ_i (−1)^{k−i}C(k,i) i^j = k!·S(j,k), one gets Δ^k g_r(0) = k! Σ_{j≥k} C(−2,j)(2/(2r+1))^j (2r+1)⁻² S(j,k). So v₂(Δ^k g_r(0)) ≥ v₂((k+1)!) + k ≥ k + 1 for k ≥ 1.
- In M_n, the k = 0 term −(−2)^{n−1}g_r(0) has v₂ exactly n − 1, while (−2)ⁿG(r) and every k ≥ 1 term have v₂ ≥ n.
- So v₂(M_n) = n − 1, and v₂(ΔⁿF(r)) = n + 1. ∎

## 5. Exact reductions. PROVED (algebra).
- (5a) e₂(K) ≥ −v₂(det A_K), since det A_K is the X⁰ coefficient. Equality holds at every K ≤ 128 except K = 11, 61, 127 (twisted), where the content is larger by 1.
- (5b) Uvarov form: det M(X) = ±det[p_i(z_b)F_b(X) + R_i(z_b)] / det V, with R_i(z) = φ((p_i(u) − p_i(z))/(u − z)).
  - Proof: M = N·diag(1/D′(z_b))·V with N_{ib} = φ(u^i/(u − z_b)), because u^j = Σ_b z_b^j ℓ_b(u) for j < K; ∏_b D′(z_b) = ±det V²; then change the rows to the p_i basis.
- (5c) Divided differences on every row plus the Leibniz rule give det M^v(0) = ±det T_{[K]}, with T_{rm} = φ₀(N_r(v)/N_{m+1}(v)) and N_k = ∏_{j<k}(v − w_j).
  - The entries of T do not depend on K, so the K-pole machine is the leading K×K block of one infinite matrix.
  - For m ≥ r, T_{rm} = [w_r, …, w_m]F̃ (pole data). For m < r, T_{rm} = φ₀(∏_{j=m+1}^{r−1}(v − w_j)) (moments).
- (5d) The pivots of T are the Newton–Padé residuals.
  - Let P_k be the degree-k polynomial with φ₀(P_k/(v − w_j)) = 0 for j < k; it is the k-th orthogonal polynomial of φ₀/N_k.
  - Then τ_k = ρ_k/∏_{j<k}(w_k − w_j), with ρ_k = φ₀(P_k(v)/(v − w_k)).
  - Hence v₂(τ_k) = v₂(ρ_k) − 4k + 1 + s₂(k) (twisted nodes).
- (5e) u↔v invariance. For the twisted machine, T^v(K,M) − T^u(K,M) = −2K for every M, where
  - T^u(K,M) = Σ_{n<M}[−(2K−3) + 12n − 4s₂(n)] + Σ_{m<K−M}[1 + s₂(K−1−m) − s₂(m)];
  - T^v(K,M) = Σ_{q<M}[8q − 4s₂(q) − 1] + Σ_{m<K−M}[1 + s₂(K−1−m) − s₂(m) + 4m − 2K].
  - Proof: expand; all M-dependent terms cancel.
  - So the u-frame "moment species / ones / partners" description and the v-frame union below have the same minimiser and the same ties.
- (5f) Closed form at r = 0: [z₀, …, z_m]F̃ = −(8/(2m)!)·S_m, with S_m = Σ_{j=1}^{m} C(2m, m−j)β_j.
  - Proof: 1/∏_{i≠j}(i² − j²) = (−1)^j·2/((m−j)!(m+j)!).
- (5g) Mahler form of S_m.
  - S_m = Σ_{n=1}^{m} M_n(−1)ⁿ C(2m−n−1, m−n), using the identity Σ_j (−1)^j C(2m, m−j)C(j,n) = (−1)ⁿC(2m−n−1, m−n) (generating functions).
  - Substituting Lemma 3's recursion gives S_m = Σ_{k<m} (−1)^k Δ^k g(0)·Q_{k,m}, with Q_{k,m} = Σ_{p=0}^{m−k−1} 2^{m−k−1−p}C(m−1+p, p).
  - Q_{0,m} = 4^{m−1}, by Σ_{p≤n} C(n+p, p)2^{−p} = 2ⁿ.
  - So the divided-difference law at r = 0, v₂([z₀..z_m]F̃) = 1 + s₂(m), is equivalent to the congruence Σ_{j=1}^{m} C(2m, m−j)β_j ≡ 4^{m−1} (mod 2^{2m−1}). This is verified (e.g. m = 2: 44/9 − 4 = 8/9; m = 3: 4784/225 − 16 = 1184/225, v₂ = 5) but not proved: the k ≥ 1 terms cancel among themselves.

## 6. Main Theorem (unified law). VERIFIED, NOT PROVED.
In the v-scale, the 2-adic Smith form of A^v_K is the sorted union of two lists:
- the pole band κ(K, m) + 4m − 2K + c for m < K − M*, where κ(K, m) = 1 + s₂(K−1−m) − s₂(m) and c = 0 for integer poles (twisted), c = −1 for half-integer poles (Catalan);
- the moment pivots h_q (Theorem 2) for q < M*.

Here M* minimises Σ_{m<K−M} pole(m) + Σ_{q<M} h_q. If the minimum is attained twice, the top divisor is raised: by exactly 1 in the content, and by 1 or 2 in det A.

Consequently:
- e₂^T(K) = −min_M T^u(K, M) − [tie];
- e₂^C(K) = −min_M [Σ_{k<K−M} (κ(K,k) + 4k − 2K − 1) + Σ_{q<M} h^C_q] − [tie].

Evidence:
- The full Smith-form equality holds at K = 8, 11, 12, 16, 20, 24, 28, 32, 34, 36, 40, 44, 48 for both machines (`twoadic_law.py merge`).
  - This includes ties at twisted K = 11, 24, 34 and Catalan K = 12, where the top divisor is +1 (+2 at twisted K = 11).
- The content law holds at all 134 values of K for each machine (`twoadic_law.py law all`, `claw`).
- The u-frame elementary-divisor law was checked by Smith form at K = 4…56, 60, 64.
- The first-half Kummer pivots of the twisted machine: 0 violations in 400.
- The divided-difference law D_{r,m} = 2·C(2m,m)·unit: 0 violations for r, m ≤ 30.
- Catalan's pole-band determinant equals the sum of its pole-band formula at K = 8…48.

Consequences, given the theorem:
- the constant ⅙ (M* ≈ K/6 moment modes);
- the (K/3)·log₂K carry term (4Σ s₂ over the moment index: Delange);
- the exact law on K = 3·2^j (where the Kummer band values make the split clean);
- the powers-of-two drift;
- Catalan − twisted = (8/3)K + O(1).

## 7. What remains to prove
- Lemma A (pole band). For m < K − M*, the Newton-frame pole pivots have v₂ = κ(K, m) + 4m − 2K + c.
  - Route: the divided-difference law D_{r,m} = 2C(2m,m)·unit (reduced at r = 0 to the congruence in 5g), then the Kummer structure of Hankel matrices with entries 2^{1+s₂(K−1−i−j)}·unit.
- Lemma B (merge). The orthogonal-polynomial directions (moment species) and the Newton directions (pole band) are 2-adically orthogonal enough that the Smith form is the union, with the split at the tropical minimum.
  - Evidence: eliminating the pole half and the moment corner together decouples exactly (v₂(det A_EE) = K + Σε at every K tested), and the middle band then follows the pole-band law.
- Lemma C (tie). At a tie the two competing boundary terms have equal valuation, and their units cancel to exactly one extra power of 2 in the content.
  - The +2 cases in det A (K = 11, 61, 127) are exactly where content ≠ det A.
- (Superseded by sections 8–11 below: A and C are proved there; B is reformulated.)

## 8. The Krawtchouk structure (2026-09-28, second stage). PROVED.
Notation:
- N = 2K − 2 (twisted) or N = 2K − 1 (Catalan).
- x ∈ {0, …, N}, with centred variable y = x − N/2 (integers for twisted, half-integers for Catalan).
- b(x) = C(N, x).
- Q_n = monic symmetric Krawtchouk polynomials in y: y·Q_n = Q_{n+1} + λ_n·Q_{n−1}, λ_n = n(N+1−n)/4.
- ‖Q_n‖² = Σ_x b(x)Q_n(y)² = 2^N·n!·N!/((N−n)!·4ⁿ).
- ν(n) := v₂‖Q_n‖² = N − 2n + v₂(n!) + v₂(N!) − v₂((N−n)!).

**Theorem 4 (the X-coefficient is Krawtchouk).** B_K is the Hankel matrix of a discrete measure on the poles.
- Twisted, u-scale: nodes −a², weights β_a = 4ε_a/((K−1−a)!(K−1+a)!), with ε₀ = 1 and ε_a = 2 for a ≥ 1.
  - So N!·β is 4 × the symmetric binomial on y ∈ [−(K−1), K−1], pushed to u = −y².
- Catalan, v-scale: nodes −(2j+1)², weights β_j = (2j+1)²/(2·4^{K−1}(K−1−j)!(K+j)!).
  - So this is the binomial on the half-integers, times (2y)², pushed to v = −4y².
- The v-scale monic orthogonal polynomials p_n of B are the even (twisted) or odd (Catalan) Krawtchouk polynomials:
  - twisted: p_n(−4y²) = (−4)ⁿQ_{2n}(y);
  - Catalan: y·p_n(−4y²) = (−4)ⁿQ_{2n+1}(y).
- They satisfy v·p_n = p_{n+1} + b^v_n·p_n + λ^v_n·p_{n−1} with INTEGER coefficients:
  - twisted: b^v_n = −2(4n(K−1−n) + K − 1), λ^v_n = 4n(2n−1)(K−n)(2K−1−2n);
  - Catalan: b^v_n = −[(2n+1)(2K−1−2n) + 4(n+1)(K−1−n)], λ^v_n = 4n(K−n)(2n+1)(2K−1−2n).
  - So each p_n ∈ ℤ[v] is monic. The matrix taking 1, v, …, v^{K−1} to p_0, …, p_{K−1} is unit lower triangular in GL_K(ℤ).
- In the basis p_n, B is diagonal: B = diag(β̂_0, …, β̂_{K−1}), where β̂_n is the n-th norm.
  - So the leading pivots of B^v_K are the β̂_n, and its elementary divisors are the v₂(β̂_n).
  - Only a unimodular change of basis over ℤ₂ is used; no square roots or field extensions.
- v₂(β̂_n) = κ(K,n) + 4n − 2K + 1 + c, with c = 0 (twisted) or c = −1 (Catalan).
- Proof.
  - Weights = pole value of B divided by D′(z), where:
    - D′(−a²) = (−1)^a(K−1−a)!(K−1+a)!/ε_a;
    - D′(−(2j+1)²) = 4^{K−1}(−1)^j(K−1−j)!(K+j)!/(2j+1);
    - the pole values of B are 4(−1)^a (twisted, u-scale) and (−1)^j(2j+1)/2 (Catalan, v-scale).
  - Pushforward from the lattice x ∈ {0, …, N}; the pair ±y lands on one node:
    - twisted: Σ_a β_a f(−a²) = (4/N!)·Σ_x b(x) f(−y²);
    - Catalan: Σ_j β_j f(−(2j+1)²) = (1/(4^{K−1}N!))·Σ_x b(x) y² f(−4y²).
  - Recurrences.
    - For a symmetric measure, Q_{2n} = R_n(y²) and Q_{2n+1} = y·R′_n(y²).
    - From y·Q_n = Q_{n+1} + λ_nQ_{n−1}, with s = y²:
      - s·R_n = R_{n+1} + (λ_{2n} + λ_{2n+1})R_n + λ_{2n−1}λ_{2n}R_{n−1};
      - s·R′_n = R′_{n+1} + (λ_{2n+1} + λ_{2n+2})R′_n + λ_{2n}λ_{2n+1}R′_{n−1}.
    - Substituting v = −4s and p_n = (−4)ⁿR_n (resp. R′_n) gives b^v_n = −4·(the sum) and λ^v_n = 16·(the product).
    - With λ_k = k(N+1−k)/4 these are the integers above.
  - Norms.
    - Twisted: β̂^u_n = (4/N!)·‖Q_{2n}‖², so v₂ = 2 − v₂(N!) + ν(2n) = κ(K,n) + 1 in the u-scale. The v-scale adds 4n − 2K, since B^v = 4^{−K}·diag(4^i)·B^u·diag(4^j).
    - Catalan: β̂_n = 16ⁿ·‖Q_{2n+1}‖²/(4^{K−1}N!), so v₂ = 4n + ν(2n+1) − 2K + 2 − v₂(N!) = κ(K,n) + 4n − 2K.
    - Both simplifications use v₂(n!) = n − s₂(n) and s₂(2m+1) = 1 + s₂(m). ∎

**Theorem 5 (the 2-adic Catalan constant; the split).**
- Let g(y) = (2y+1)⁻² and G(y) = (−1)^y β_{|y|}. This is the even extension to ℤ; it satisfies G(y) + G(y+1) = −g(y) on all of ℤ, since g(−1−y) = g(y).
- Define G₁ := −Σ_{n≥0} (−1)ⁿΔⁿg/2^{n+1}.
  - It converges on ℤ₂ (Lemma 3 gives v₂(Δ^k g) ≥ v₂((k+1)!) + k).
  - It is even, and it satisfies the same equation.
  - Its smoothness: v₂(Δ^m G₁(x)) ≥ m − 1 + v₂((m+1)!) for all x ∈ ℤ.
- Hence G = G₁ + C·(−1)^y with C := −G₁(0) = Σ_n (−1)ⁿΔⁿg(0)/2^{n+1} and v₂(C) = −1.
  - C is the 2-adic (Boole-regularised) value of Σ(−1)^k/(2k+1)², the 2-adic Catalan constant.
  - β_{2^N} → 0 is the statement G₁(2^N) → G₁(0) = −C.
- Split: A(0) = −C·B + A(C) exactly. At X = C the pole values are smooth:
  - twisted: −4G₁(a);
  - Catalan: (2j+1)ΔG₁(j)/4 (because G₁ + g/2 = −ΔG₁/2).
- Twisted only: −4G₁(a) = 4C − 4(G₁(a) − G₁(0)).
  - The constant part adds the Hankel matrix 4C·h_{s−K+1}(z) (complete homogeneous symmetric functions of the nodes).
  - That is, it extends the moment functional by μ_u(−1) = 4C, the 2-adic analogue of the real μ(−1) = φ₀(1/u) = 4G.
- Catalan has no constant part: its smooth part is a Δ.

## 9. Lemma A. PROVED (Catalan: the full band; twisted: the first half).
- **(A-Cat)** The Catalan X = 0 pole Hankel H_ω (v-scale, K×K) has leading pivots of valuation κ(K,n) + 4n − 2K − 1 for every n < K, and SNF(H_ω) is the sorted band.
- **(A-Tw)** The twisted X = 0 pole Hankel (u-scale) has leading pivots of valuation κ(K,n) for n < ⌊K/2⌋.
- For 2k − 2 ≤ K − 1 the leading k×k block of A equals that of the pole Hankel. So the first-half Kummer pivots of A itself (previously checked 400/400) are proved for both machines.

Proof.
1. **Base plus perturbation.** The pole measure is ω = ω⁰(1 + ε(x)), with ω⁰ ∝ the Krawtchouk measure of Theorem 4 and ε(x) = ±(−1)^x·S(x). The smooth S is:
   - Catalan: S(x) = ΔG₁(x−K)/(2C). Here v₂(2C) = 0 and v₂(Δ^kS) ≥ k + v₂((k+2)!).
   - Twisted: S(x) = G₁(x−K+1)/C. Here S(0) = −1 and v₂(Δ^kS) ≥ k + v₂((k+1)!) for k ≥ 1.
   - (The ± signs and the evenness under y → −y follow from G₁ even and ΔG₁(−1−j) = −ΔG₁(j).)
2. **Criterion (diagonal dominance, over ℚ₂).** In the orthogonal polynomial basis p_n of ω⁰ (Theorem 4: a unimodular change of basis in the v-scale), Gram(ω) = D + E.
   - D = diag(h⁰_0, h⁰_1, …) is the base norms, and E_mn = Σ_nodes ω⁰·ε·p_m·p_n is the perturbation.

   **Lemma 9.1.** Let d_0, …, d_{n−1} ∈ ℚ₂^×, δ_i := v₂(d_i), and let E ∈ M_n(ℚ₂) satisfy v₂(E_ij) > (δ_i + δ_j)/2 for all i, j (diagonal included; E need not be symmetric). Put M = diag(d) + E. Then:
   - (i) for all I, J with |I| = |J|: v₂ det M[I, J] ≥ ½(Σ_{i∈I} δ_i + Σ_{j∈J} δ_j);
   - (ii) for every I: det M[I, I] = (Π_{i∈I} d_i)·(1 + η_I) with v₂(η_I) ≥ 1;
   - (iii) the k-th leading pivot of M has valuation δ_k, and the elementary divisors of M are the δ_i sorted.

   Proof.
   - det M[I, J] = Σ_σ ±Π_{i∈I} M_{i,σ(i)}, over bijections σ: I → J.
   - A factor with σ(i) = i is d_i + E_ii, of valuation exactly δ_i.
   - A factor with σ(i) ≠ i is E_{i,σ(i)}, of valuation > (δ_i + δ_{σ(i)})/2.
   - Since σ is a bijection, Σ_{i∈I} (δ_i + δ_{σ(i)})/2 = ½(Σ_I δ + Σ_J δ). This gives (i).
   - For I = J, the identity term is Π(d_i + E_ii) = Π d_i·Π(1 + E_ii/d_i), with every v₂(E_ii/d_i) ≥ 1.
   - Any σ ≠ id has a nonempty moved set S with σ(S) = S. The factors over S have total valuation strictly greater than the integer Σ_{i∈S} δ_i, hence at least that integer + 1. The fixed factors contribute exactly Σ δ_i. This gives (ii).
   - (iii): the leading minors come from (ii). The gcd of the k×k minors is ≥ the sum of the k smallest δ by (i), and it is attained by (ii) at the principal minor on those k indices. ∎

   - Lemma 9.1 uses no square roots and no field extension; the half-integers are only a bookkeeping device for comparing valuations.
   - Applied to D + E, it says the leading pivots are those of D. For Catalan (v-scale, unimodular basis) the Smith form is also that of D.
   - It remains to verify the hypothesis v₂(E_mn) > ½(v₂h⁰_m + v₂h⁰_n).
3. **Mirror.**
   - On the lattice (the pushforward of Theorem 4), E_mn = κσ_mσ_n·⟨Q_a, εQ_b⟩ and h⁰_m = κσ_m²·‖Q_a‖², where:
     - ⟨f, g⟩ := Σ_x b(x) f g;
     - κ is a common constant and σ_m is a power of ±4 (or ±1);
     - (a, b) = (2m+1, 2n+1) for Catalan and (2m, 2n) for twisted.
   - The constants cancel, so the hypothesis of Lemma 9.1 is exactly v₂⟨Q_a, εQ_b⟩ > ½(ν(a) + ν(b)).
   - From Σ_n C(N,n)K_n(x)tⁿ = (1−t)^x(1+t)^{N−x} for the p = ½ Krawtchouk polynomials: (−1)^x·Q_b = γ_b·Q_{N−b} on the nodes, with γ_b = (−1)^N·2^{N−2b}·b!/(N−b)!.
   - From the formula for ν, v₂(γ_b) = ½(ν(b) − ν(N−b)).
   - ε = (±1)·(−1)^x·S(x) (step 1). So the hypothesis for (a, b) is equivalent to v₂⟨Q_a, S·Q_c⟩ > ½(ν(a) + ν(c)), with c = N − b.
4. **Smoothness (balanced valuations, over ℚ₂).**
   - The functions on X = {0, …, N} have the basis Q_0, …, Q_N. The monic Q_{N+1} vanishes on X: it is the node polynomial of an (N+1)-point measure.
   - For a function f on X let M^f be its multiplication matrix in the monic basis: f·Q_c = Σ_a M^f_ac·Q_a on X.
     - Then ⟨Q_a, f·Q_c⟩ = M^f_ac·‖Q_a‖², and f ↦ M^f is a ring homomorphism.
     - M^x = 𝒥, the finite monic Jacobi matrix: 𝒥_{c+1,c} = 1, 𝒥_cc = N/2, 𝒥_{c−1,c} = λ_c.
   - **Balanced valuation.** For an (N+1)×(N+1) matrix M over ℚ₂ put β(M) := min_{a,c} [v₂(M_ac) + ½(ν(a) − ν(c))] ∈ ½ℤ ∪ {∞}.
     - β(MM′) ≥ β(M) + β(M′), because ½(ν(a) − ν(b)) + ½(ν(b) − ν(c)) telescopes.
     - β(M + M′) ≥ min(β(M), β(M′)), and β(tM) = v₂(t) + β(M).
   - The step 3 quantity is v₂⟨Q_a, SQ_c⟩ − ½(ν(a) + ν(c)) = v₂(M^S_ac) + ½(ν(a) − ν(c)), because v₂⟨Q_a, SQ_c⟩ = v₂(M^S_ac) + ν(a).
   - **Entries of 𝒥 − tI, t ∈ ℤ.** Since ‖Q_c‖² = λ_c‖Q_{c−1}‖², we have ν(c) − ν(c−1) = v₂(λ_c). So:
     - the (c+1, c) entry has balanced valuation ½v₂(λ_{c+1});
     - the (c−1, c) entry has v₂(λ_c) − ½v₂(λ_c) = ½v₂(λ_c);
     - the diagonal has v₂(N/2 − t).
     - Catalan (N odd): v₂(N/2 − t) = −1 and v₂(λ_c) = v₂(c(2K−c)) − 2 ≥ −2, so β(𝒥 − tI) ≥ −1.
     - Twisted (N even): v₂(N/2 − t) ≥ 0, and c(2K−1−c) is even, so v₂(λ_c) ≥ −1 and β(𝒥 − tI) ≥ −½.
   - **Binomial polynomials.** C(𝒥, k) = (1/k!)·Π_{t<k}(𝒥 − tI).
     - β(C(𝒥, k)) ≥ −k − v₂(k!) (Catalan) or ≥ −k/2 − v₂(k!) (twisted).
     - It is banded: C(𝒥, k)_ac = 0 when |a − c| > k.
   - **Finite Mahler expansion.** S = Σ_{k≤N} s_k·C(x, k) exactly on X, with s_k = Δ^kS(0). Hence M^S = Σ_k s_k·C(𝒥, k).
   - Catalan:
     - a is odd and c is even, so a ≠ c and the k = 0 term s₀·I does not contribute;
     - each k ≥ 1 term has balanced valuation ≥ v₂(s_k) − k − v₂(k!) ≥ v₂((k+1)(k+2)) ≥ 1, using v₂(s_k) ≥ k + v₂((k+2)!).
   - Twisted, leading ⌊K/2⌋ block:
     - a = 2m and c = N − 2n with m, n < ⌊K/2⌋ give c − a ≥ 2, so only k ≥ 2 contribute;
     - each has balanced valuation ≥ k + v₂((k+1)!) − k/2 − v₂(k!) = k/2 + v₂(k+1) ≥ 1.
   - So the hypothesis of Lemma 9.1 holds with margin ≥ 1. ∎
- **Verified link by link.**
  - proof_checks/th_lemmaA7.py, K = 4…20: the mirror identity (0 violations), the norm formula, the λ bound, the Mahler bound (tight), and the final margins.
  - proof_checks/th_lemmaA_balanced.py, K = 4, 5, 8, 11:
    - multiplication by x on X equals 𝒥;
    - the balanced bounds on C(𝒥, k) hold for all k ≤ N;
    - M^S = Σ s_k·C(𝒥, k) exactly;
    - the margins are ≥ 1.
  - The Catalan margin is exactly 1, attained at (m, n) = (0, K−1) through the k = 1 term. th_lemmaA5.py shows the same margins directly for K ≤ 24.
- **Why the twisted pole Hankel alone leaves the band past K/2.**
  - The k = 0 term s₀ = −1 lands on the anti-diagonal m + n = K − 1. This is the vanishing weight β₀ = 0, i.e. G₁(0) = −C.
  - The full matrix recovers the band through the μ(−1) = 4C extension (section 11).

## 10. Lemma C. PROVED (content form), given the coefficient law.
Setup:
- p(X) = det(A + XB) = q(X − C), where q(Y) = det(A(C) + YB).
- Put t_k := v₂(q_k) − k and m₀ := min_M T(K, M) (native scale).

**Coefficient law (Q).**
- t_k = m₀ for k ∈ {K − M : M minimises T}.
- t_k ≥ m₀ + 1 for every other k.

**Claim.** Given (Q):
- e₂ = −m₀ when the minimiser is unique, and e₂ = −m₀ − 1 at a tie;
- v₂(det A) = m₀ when the minimiser is unique, and ≥ m₀ + 1 at a tie.

Proof.
- Write C = u/2 with u a 2-adic unit. Then p_j = 2^j(−u)^{−j} Σ_k Q_k(−u)^k C(k, j), with Q_k = q_k2^{−k} and v₂(Q_k) = t_k.
- Unique minimiser: v₂(p₀) = m₀, and v₂(p_j) ≥ j + m₀ for j ≥ 1.
- Ties:
  - T(K,·) is strictly convex, because T(K,M+1) − T(K,M) = h_M − b_{K−1−M}.
    - h increases by at least 4 per step: 4 + 4v₂(q+1) (twisted), 4 + 3v₂(q+1) + v₂(q+2) (Catalan).
    - The band b increases by at least 2 per step.
  - So a tie is exactly two consecutive M, i.e. k₁ and k₂ = k₁ + 1 of opposite parity, and every other t is ≥ m₀ + 6.
  - Then v₂(p₀) ≥ m₀ + 1, since two unit multiples of 2^{m₀} add to an even multiple.
  - And v₂(p₁) = 1 + v₂(k₁Q̃₁ + k₂Q̃₂ + …) = m₀ + 1, since exactly one of k₁, k₂ is odd.
  - Also v₂(p_j) ≥ m₀ + 1 for all j. ∎

Remarks and checks:
- The +1 or +2 of det A at a tie is the next 2-adic digit of Q̃₁ + Q̃₂, which the tropical data does not fix. This is where K = 11, 61, 127 come from.
- (Q) is verified by th_lemmaC1.py: twisted and Catalan, K = 6, 8, 11, 12, 16, 20, 24, 28, 32, 34 (20 cases, 5 ties).
  - t = m₀ at the minimisers, all others ≥ m₀ + 1, and the content comes out exactly.
  - The strict tropical inequality t_k ≥ T(K, K−k) fails for M ≳ K/4. That is harmless: those t lie more than 1000 above m₀.

## 11. Lemma B, reduced. VERIFIED, NOT PROVED.
**The model.** A_model = −C·B + P_ext.
- P_ext is the true moment corner (Catalan), or the moment corner plus the μ(−1) = 4C term (twisted).
- SNF(A_model) = SNF(A) at K = 8, 11, 12, 16, 20, 24, 28, 32 for both machines, ties included (th_lemmaA4.py, th_lemmaA6.py).
  - Without the μ(−1) term the twisted model fails.
  - The Catalan smooth part is invisible.

**In the Krawtchouk basis.** A_model = D + P′.
- D = diag(−C·h^B_m), with v₂ = band b_m.
- P′_mn = Ψ(p_m p_n), where Ψ(f) = φ₀(⌊f/N_K⌋) (plus the residue term for twisted).
  - P′ vanishes for m + n < K (Catalan) or m + n < K − 1 (twisted).
  - Every nonzero entry has v₂ = −1 (K = 12 printed in full).
- Exact expansions:
  - det(D + P′) = Σ_U Π_{m∉U} d_m · det P′_U;
  - q_{K−M} = Σ_{|U|=M} Π_{m∉U} β̂_m · det(P′ + E′)_U.

**(B1) Top-block identity (verified).** The reversed top block Z_ij = Ψ(p_{K−1−i}p_{K−1−j}) has leading pivots exactly h₀, h₁, …, h_{r−1}.
- r = 2, 4, 4, 6, 7, 8 at K = 8, 12, 16, 20, 24, 32 (r ≈ K/4).
- M* < r in every case.
- So the top M*×M* block of the moment corner, in the Krawtchouk basis, is unitriangularly equivalent to the CDH moment Hankel.

**(B*) Weighted minor bound (conjectured).** For all I, J with |I| = |J|:
- v₂ det P′[I,J] − ½(Σ_I b + Σ_J b) ≥ m₀ − Σ_all b,
- with equality only at I = J = [K−M, K) for M a tropical minimiser.
- Evidence (principal minors, K = 12, all U with |U| ≤ 6): the minimum is attained only at the top blocks of the minimisers (Catalan K = 12: M = 2 and 3 tie at −29; twisted K = 12: M = 3 alone at −33).

**What (B1) and (B*) give.**
- The dominant term of each q_{K−M} (M near M*) is the top block, and every other term, including those carrying the Catalan perturbation E′ (margin ≥ 1 by section 9), is at least one power of 2 higher.
- That is the coefficient law (Q), hence the content law via section 10.
- The full Smith-form statement (the Main Theorem) also needs (B*) for non-principal minors.
- (Superseded by section 12, which proves the Main Theorem another way. The moment corner alone is the wrong object: the smooth pole part at X = C belongs with it.)

## 12. THE PENCIL THEOREM: the Main Theorem, proved (2026-09-28, third stage)
**Idea.** At X = C (the 2-adic Catalan constant) the machine becomes "the CDH moment functional divided by N_K", read 2-adically. So A(X) = A(C) + (X − C)·B is a pencil of two Hankel forms:
- A(C) has Smith form {h_q}, the moment pivots, in the CDH basis;
- B has Smith form {b_m + 1}, the Krawtchouk band, in the Krawtchouk basis;
- mod 2 the two bases are related by the Pascal matrix.

Cauchy–Binet then gives everything.

(Rigour pass 2026-09-29: sections 8, 9 and 12 now work entirely in ℚ₂ with integral monic bases. The earlier orthonormal-basis arguments needed √h_q and √λ_n, which do not lie in ℚ₂, and did not show that the resulting factor L is integral. Lemma 12.1 and Theorem 7 replace them; section 9 uses Lemma 9.1 and balanced valuations.)

**Lemma 12.1 (monic Jacobi matrices over ℤ₂).**
- Setup.
  - P_0, P_1, … are the monic CDH orthogonal polynomials in the v-scale: v·P_n = P_{n+1} + b_nP_n + λ_nP_{n−1}, with h_n := φ₀(P_n²) = ½·λ_1⋯λ_n.
    - Catalan: b_n = 4(2n+1)(n+1), λ_n = 16n³(n+1).
    - Twisted: b_n = 8n² + 8n + 3, λ_n = 16n⁴.
    - All are integers (Theorem 1), so P_n ∈ ℤ[v].
  - 𝒥 is the infinite tridiagonal matrix with 𝒥_{n+1,n} = 1, 𝒥_nn = b_n, 𝒥_{n−1,n} = λ_n.
    - Column convention: v·P_q = Σ_p 𝒥_pq P_p.
    - 𝒥 is not symmetric; no normalisation is taken.
- **(a)** For every polynomial f: f(v)·P_q = Σ_p f(𝒥)_pq P_p (a finite sum). Hence φ₀(f·P_p·P_q) = h_p·f(𝒥)_pq.
  - Proof: induction on the degree of f, then orthogonality φ₀(P_pP_r) = h_pδ_pr.
  - The symmetry h_p·f(𝒥)_pq = h_q·f(𝒥)_qp follows automatically.
- **(b) Rescaling by a rational diagonal matrix.** Put ρ = 4, v* = 0 (Catalan) or ρ = 2, v* = 1 (twisted), and S := diag(ρ^{−n})_{n≥0}. Then 𝒥̂ := S^{−1}(𝒥 − v*·I)S has
  - Catalan: 𝒥̂_{n+1,n} = 4, 𝒥̂_nn = 4(2n+1)(n+1), 𝒥̂_{n−1,n} = λ_n/4 = 4n³(n+1);
  - twisted: 𝒥̂_{n+1,n} = 2, 𝒥̂_nn = b_n − 1 = 2(2n+1)², 𝒥̂_{n−1,n} = λ_n/2 = 8n⁴.
  - So every entry of 𝒥̂ lies in ρℤ.
  - Each entry (𝒥̂^k)_pq is a finite sum, over lattice paths of length k, of products of k entries. Hence (𝒥̂^k)_pq ∈ ρ^kℤ.
  - Moreover ((𝒥 − v*I)^k)_pq = ρ^{q−p}·(𝒥̂^k)_pq.
- **(c) Moment decay.** Let m_k := φ₀((v − v*)^k) be the moments about v* (m_k = μ_k for Catalan). By (a) and (b), m_k = h_0·((𝒥 − v*I)^k)_00 = ½·(𝒥̂^k)_00. So:
  - v₂(m_k) ≥ 2k − 1 (Catalan);
  - v₂(m_k) ≥ k − 1 (twisted).
- **Definition (the 2-adic extension).**
  - Let 𝒜 be the power series f = Σ_k f_k(v − v*)^k with (f_k) bounded in ℚ₂.
  - Put φ₀^{(2)}(f) := Σ_k f_k m_k.
  - By (c) this converges, |φ₀^{(2)}(f)| ≤ 2·sup_k|f_k|, φ₀^{(2)} is linear, and it equals φ₀ on polynomials.
- Checked: proof_checks/th_thm7_rigorous.py (the rescaled matrices are integral and divisible by ρ).

**Theorem 6 (2-adic Stieltjes identity). PROVED (both machines).**
- **(a) Catalan, v-scale.** For every j ≥ 0: F_j(C) = −Σ_{k≥0} μ_k w_j^{−k−1} = φ₀^{(2)}(1/(v − w_j)), where
  - w_j = −(2j+1)²;
  - μ_k = 4^k(2^{2k+2}−1)|B_{2k+2}| = (k+1)T_{2k+1}/2 (T = tangent numbers), with v₂(μ_k) = 2k − 1 (von Staudt; only v₂ ≥ 2k − 1, from Lemma 12.1(c), is used);
  - the Taylor coefficients of 1/(v − w_j) at 0 are −w_j^{−k−1}, which are units, so 1/(v − w_j) ∈ 𝒜.
  - Proof.
    1. **The derivative form.** f_D(y) := −½Σ_{n≥0} E_n(0)·g^{(n)}(y)/n! = −½Σ_n a_n(2y+1)^{−n−2}, with a_n := E_n(0)(−2)ⁿ(n+1), using g^{(n)}(y)/n! = (−2)ⁿ(n+1)(2y+1)^{−n−2}.
       - E_n(0) = 0 for even n ≥ 2, and E_{2k+1}(0) = −(2^{2k+2}−1)B_{2k+2}/(k+1). So v₂(a_{2k+1}) = 2k + 1.
       - For |y| ≤ r in ℂ₂ with 1 < r < 2, |2y| < 1, so 2y + 1 is a unit and |(2y+1)^{−n−2}| = 1.
       - Hence the series converges uniformly on |y| ≤ r. f_D is analytic there, with Gauss norm ‖f_D‖_r ≤ max_n|a_n|.
    2. **Functional equation.**
       - f_D(y+1) = Σ_m f_D^{(m)}(y)/m!, because the step 1 is inside the disk of analyticity.
       - The double series Σ_{n,m} E_n(0)·g^{(n+m)}(y)/(n!m!) converges absolutely: its terms are E_n(0)·C(n+m, n)·g^{(N)}/N! with N = n + m, of size ≤ 2(n+1)·2^{−N}.
       - (1 + e^t)·Σ_n E_n(0)tⁿ/n! = 2 gives Σ_{n+m=N} E_n(0)/(n!m!) + E_N(0)/N! = 2δ_{N,0}.
       - Hence f_D(y) + f_D(y+1) = −g(y).
    3. **Uniqueness.** d := G₁ − f_D satisfies d(y+1) = −d(y) on ℤ. So d = c·(−1)^y, and Δ^m d(0) = c(−2)^m has v₂ = v₂(c) + m.
       - But v₂(Δ^m G₁(0)) ≥ m − 1 + v₂((m+1)!) (Theorem 5).
       - And v₂(Δ^m f_D(0)) ≥ v₂(m!) + m·log₂r − log₂‖f_D‖_r, since f_D = Σ a′_i y^i with |a′_i| ≤ ‖f_D‖_r·r^{−i} and Δ^m(y^i)(0) = m!·S(i,m).
       - Taking r close to 2, both bounds exceed m + v₂(c) for large m unless c = 0. So G₁ = f_D.
    4. **Evaluation.**
       - ΔG₁ = −2G₁ − g = Σ_k E_{2k+1}(0)·g^{(2k+1)}/(2k+1)!, the n = 0 term cancelling g.
       - Using g^{(2k+1)}(y) = −2^{2k+1}(2k+2)!(2y+1)^{−2k−3} and B_{2k+2} = (−1)^k|B_{2k+2}|, this gives F_j(C) = (2j+1)ΔG₁(j)/4 = Σ_k (−1)^k μ_k (2j+1)^{−2k−2}.
       - That series is −Σ_k μ_k w_j^{−k−1}. ∎
  - Corollary (a fast formula for the 2-adic Catalan constant): C = ½ + 2Σ_k (−1)^k μ_k = ½ + Σ_k (−1)^k (k+1)T_{2k+1}.
- **(b) Twisted, v-scale.** For every a ≥ 0: F_a(C)/4 = −G₁(a) = Σ_{k≥0} (−1)^k m_k (1+4a²)^{−k−1} = φ₀^{(2)}(1/(v + 4a²)), where
  - m_k = φ₀((v−1)^k) are the centred moments, with v₂(m_k) ≥ k − 1 (Lemma 12.1(c); in fact = k − 1);
  - here the moment functional sits 2-adically near v = 1 (the rescaled 𝒥 − I is ≡ 0 mod 2), while the nodes −4a² sit near 0.
  - Proof (complete).
    1. **Boole residue form on polynomials.**
       - For f ∈ ℚ[v] put F_f(m) = f(v_m) + 2v_m f′(v_m), with v_m = −(2m+1)². These are the double poles of sech²(πt) at t = i(m+½): F_f(m) = d/dζ[ζ·f(4ζ²)] at ζ = i(m+½), the local term of the weight at that pole.
       - Put ℬ(F) := ½Σ_{n≥0}(−½)ⁿΔⁿF(0), the Boole value (1+E)^{−1}F(0). For polynomial F this is the value at 0 of the unique polynomial solution of X(m) + X(m+1) = F(m).
       - Claim: ℬ(F_f) = φ₀(f).
         - For f = v^e, F_f(m) = (2e+1)(−1)^e(2m+1)^{2e}.
         - X(m) = 4^e E_{2e}(m+½)/2 solves X(m) + X(m+1) = (2m+1)^{2e}, because E_n(x+1) + E_n(x) = 2xⁿ.
         - So ℬ((2m+1)^{2e}) = 4^e E_{2e}(½)/2 = E_{2e}/2.
         - Hence ℬ(F_{v^e}) = (2e+1)(−1)^e E_{2e}/2 = (2e+1)|E_{2e}|/2 = μ_v(e).
    2. **Two continuous extensions.** Fix ½ < ρ < ρ′ < 1 and let 𝒜 be the functions analytic on the closed disk |v−1| ≤ ρ′.
       - E1(f) := Σ_k f̃_k m_k, where f = Σ f̃_k(v−1)^k (this is φ₀^{(2)} of Lemma 12.1).
         - |f̃_k| ≤ ‖f‖_{ρ′}·ρ′^{−k} and |m_k| ≤ 2^{1−k}.
         - So the series converges and |E1(f)| ≤ 2‖f‖_{ρ′}, because 2ρ′ > 1.
       - E2(f) := ℬ(F_f).
         - For |m| ≤ r := 2√ρ (and 1 < r < 2), |v_m − 1| = |4m² + 4m + 2| ≤ max(r²/4, r/4, ½) ≤ ρ. So F_f is analytic on |m| ≤ r.
         - ‖F_f‖_r ≤ ‖f‖_ρ + ‖f′‖_ρ ≤ 2‖f‖_{ρ′}/ρ (Cauchy estimate).
         - Its Mahler coefficients satisfy |ΔⁿF(0)| ≤ |n!|·r^{−n}·‖F‖_r.
         - The Boole series has terms of size ≤ 2·2ⁿ|n!|r^{−n}‖F‖_r = 2·2^{s₂(n)}r^{−n}‖F‖_r → 0.
         - So E2 is continuous on 𝒜.
       - E1 = E2 = φ₀ on polynomials (step 1).
       - f_a := 1/(v+4a²) lies in 𝒜: its pole is at distance |1 + 4a²|₂ = 1. Its Taylor truncations at v = 1 converge to it in ‖·‖_{ρ′}, with tail ≤ ρ′^{N+1}.
       - Hence E1(f_a) = E2(f_a).
    3. **Evaluation.**
       - F_{f_a}(m) = (4a² + s²)/(4a² − s²)² with s = 2m+1, which equals ½[(s−2a)^{−2} + (s+2a)^{−2}] = ½[g(m−a) + g(m+a)].
       - ℬ commutes with shifts, and ℬ(g(·+x)) = ½Σ(−½)ⁿΔⁿg(x) = −G₁(x) by the definition of G₁.
       - So E2(f_a) = −½[G₁(−a) + G₁(a)] = −G₁(a), since G₁ is even. ∎
  - Corollary: C = φ₀^{(2)}(1/v) = Σ_k (−1)^k m_k. So the twisted μ_v(−1) = C (μ_u(−1) = 4C) is literally the (−1)-st moment.
- Checked numerically to 2^{−179} (Catalan, j < 12) and 2^{−197} (twisted, a < 10): proof_checks/th_stieltjes.py and th_stieltjes_tw.py. Both C formulas agree with the Mahler-series C to 2^{−197}.

**Theorem 7 (the moment side). Everything is over ℚ₂; no square roots.**
- Let π := 1/N_K, expanded at v* (Lemma 12.1): π = Σ_k c_k(v − v*)^k.
  - c_k ∈ ℤ₂, because N_K ∈ ℤ[v] and N_K(v*) is odd: N_K(0) = Π(2j+1)² for Catalan, N_K(1) = Π(1+4a²) for twisted.
  - c_0 = 1/N_K(v*) is a unit.
- **(i)** A_K(C)_ij = φ₀^{(2)}(v^{i+j}·π) (v-scale).
- **(ii)** Let G_pq := φ₀^{(2)}(P_p·P_q·π) for p, q < K. This is the Gram matrix of A_K(C) in the CDH basis. Then G = LΛLᵀ, where:
  - L ∈ M_K(ℤ₂) is unit lower triangular, with v₂(L_pq) ≥ 3(p−q) + 2 (Catalan) or ≥ 3(p−q) + 1 (twisted) for p > q. In particular L ≡ I (mod 16).
  - Λ = diag(h_p·u_p) with u_p ∈ ℤ₂^× and u_p ≡ c_0 (mod ρ).
- **(iii)** Hence the leading pivots of A_K(C) are the h_p·u_p, and its Smith form is {v₂(h_p) : p < K}:
  - Catalan: v₂(h_p) = 4p + 3v₂(p!) + v₂((p+1)!) − 1;
  - twisted: v₂(h_p) = 8p − 4s₂(p) − 1.
  - With π = 1 this is Theorem 2 again.

Proof.
- **(i)**
  - Partial fractions: v^{i+j}/N_K = P + Σ_b r_b/(v − w_b), with P a polynomial and r_b = w_b^{i+j}/N_K′(w_b).
  - Each 1/(v − w_b) lies in 𝒜, since |w_b − v*|₂ = 1. For Catalan w_b = −(2b+1)², a unit; for twisted 1 − w_b = 1 + 4b², a unit.
  - φ₀^{(2)} is linear on 𝒜 and equals φ₀ on P.
  - Theorem 6 gives φ₀^{(2)}(1/(v − w_b)) = F_b(C), the v-scale pole value at X = C.
  - So φ₀^{(2)}(v^{i+j}π) = φ₀(P) + Σ_b r_b·F_b(C) = A_K(C)_ij.
- **(ii)**
  - P_p·P_q is a polynomial and the m_k tend to 0 (Lemma 12.1(c)), so the absolutely convergent double series may be rearranged. By Lemma 12.1(a), (b):
    - G_pq = Σ_k c_k·φ₀((v − v*)^k P_p P_q) = h_p·Σ_k c_k((𝒥 − v*I)^k)_pq = h_p·ρ^{q−p}·Π̃_pq, where Π̃ := Σ_k c_k·𝒥̂^k.
  - By Lemma 12.1(b) the series for Π̃ converges entrywise in ℤ₂, and Π̃ ≡ c_0·I (mod ρ).
  - **Gaussian elimination in ℤ₂** on the leading block Π̃_[K] = [Π̃_pq]_{p,q<K}:
    - invariant: the active submatrix stays ≡ c_0·I (mod ρ);
    - each pivot is ≡ c_0, a unit, and each multiplier is ≡ 0 (mod ρ);
    - subtracting a multiplier times a pivot-row entry changes the other entries by ≡ 0 (mod ρ).
    - Result: Π̃_[K] = 𝐋𝐔, with 𝐋 ∈ M_K(ℤ₂) unit lower triangular and ≡ I, and 𝐔 ∈ M_K(ℤ₂) upper triangular and ≡ c_0·I (mod ρ).
  - **Undo the rescaling.** With the diagonal matrices H := diag(h_p·ρ^{−p}) and Q := diag(ρ^q) over ℚ₂, G = H·Π̃_[K]·Q = (H𝐋H^{−1})·(H𝐔Q) =: L·U′.
    - L is unit lower triangular, with L_pq = (h_p/h_q)·ρ^{q−p}·𝐋_pq for p > q.
    - v₂(h_p/h_q) = Σ_{q<n≤p} v₂(λ_n), which is ≥ 5(p−q) for Catalan (v₂(16n³(n+1)) ≥ 5) and ≥ 4(p−q) for twisted.
    - With v₂(𝐋_pq) ≥ v₂(ρ), this gives the bounds on v₂(L_pq) stated in (ii).
    - U′ is upper triangular with diagonal h_p·𝐔_pp.
  - **Symmetrise.** G is symmetric, and its leading minors det G_[k] = Π_{p<k} h_p𝐔_pp are nonzero. A matrix with nonzero leading minors has a unique factorisation (unit lower)·(diagonal)·(unit upper). Hence U′ = ΛLᵀ with Λ = diag(h_p·𝐔_pp) and u_p := 𝐔_pp.
- **(iii)**
  - The CDH basis differs from the monomial basis by an integral unitriangular matrix (Theorem 2). That preserves leading minors and elementary divisors.
  - G = LΛLᵀ with L ∈ GL_K(ℤ₂), so G has the elementary divisors of Λ. ∎
- Checked: proof_checks/th_thm7_rigorous.py (Catalan and twisted, K = 8, 12):
  - 𝒥̂ is integral and divisible by ρ;
  - G = h_p·ρ^{q−p}·Π̃ (to truncation precision, 2^{−146} or better);
  - 𝐋 ≡ I and 𝐔 ≡ c_0·I (mod ρ);
  - the bounds on v₂(L_pq) hold;
  - G = LΛLᵀ and v₂(Λ_p) = v₂(h_p).
- End to end: th_catalan_struct.py (K = 6, 8, 12, 16, 20, 24) and th_twist_struct.py (K = 6, 8, 11, 12, 16) give SNF(A_K(C)) = {v₂(h_q)}.

**Theorem 8 (Pascal transition).**
- **Krawtchouk polynomials** (the orthogonal polynomials of B; the recurrences of Theorem 4):
  - Catalan: b^v_n = −[(2n+1)(2K−1−2n) + 4(n+1)(K−1−n)] is odd, and λ^v_n = 4n(K−n)(2n+1)(2K−1−2n) ≡ 0 (mod 4). So p_{n+1} = (v − b^v_n)p_n − λ^v_n p_{n−1} ≡ (v+1)p_n (mod 2), and p_m ≡ (v+1)^m (mod 2).
  - Twisted: b^v_n = −2(4n(K−1−n) + K − 1) is even, and λ^v_n ≡ 0 (mod 4). So p_m ≡ v^m (mod 2).
- **CDH polynomials** (Lemma 12.1):
  - Catalan: b_n ≡ 0 (mod 4) and λ_n ≡ 0 (mod 16), so P_q ≡ v^q (mod 4);
  - twisted: b_n ≡ 1 (mod 2) and λ_n ≡ 0 (mod 16), so P_q ≡ (v+1)^q (mod 2).
- **Transition.** p_m = Σ_q T_mq P_q, where T = (Krawtchouk coefficient matrix)·(CDH coefficient matrix)^{−1}. Both factors are integral and unit lower triangular, so T is too.
  - Mod 2 the P_q reduce to a basis of 𝔽₂[v], so T mod 2 is determined by the reductions:
    - Catalan: p_m ≡ (v+1)^m = Σ_q C(m,q)v^q ≡ Σ_q C(m,q)P_q;
    - twisted: p_m ≡ v^m = ((v+1) + 1)^m ≡ Σ_q C(m,q)(v+1)^q ≡ Σ_q C(m,q)P_q.
  - Either way T ≡ [C(m,q)] (mod 2).
- **Corner minors.** det[C(a+i, j)]_{i,j<M} = 1 for every a ≥ 0 and M ≥ 1.
  - Proof: for i = M−1 down to 1 replace row i by row i minus row i−1. Pascal's rule turns entry (i, j) into C(a+i−1, j−1), which vanishes at j = 0. Expanding along column 0 (whose only nonzero entry is C(a, 0) = 1 in row 0) leaves det[C(a+i′, j′)]_{i′,j′<M−1}. Induct.
  - Hence det T[[a, a+M), [0, M)] ≡ 1 (mod 2).

**MAIN THEOREM (both machines, all K). PROVED.** In the Krawtchouk basis (v-scale), A_K(X) = RΛRᵀ + (X − C)·diag(β̂_m), where:
- R := T·L, with T from Theorem 8 (p_m = Σ_q T_mq P_q) and L from Theorem 7(ii).
  - Indeed the Gram matrix of A_K(C) in the Krawtchouk basis is T·G·Tᵀ = (TL)Λ(TL)ᵀ.
  - R is integral and unit lower triangular.
  - R − T = T(L − I) ≡ 0 (mod 16), so R ≡ T ≡ Pascal (mod 2).
- Λ comes from Theorem 7(ii), with v₂(Λ_q) = v₂(h_q).
- v₂(β̂_m) = b_m + 1, with b_m = κ(K,m) + 4m − 2K + c (Theorem 4).
- Every object here is a matrix over ℚ₂ (R and T over ℤ₂). The consequences below use only Cauchy–Binet and valuations.

Consequences:
1. **Coefficient law.** q(Y) = det(A(C) + YB) has v₂(q_{K−M}) = Σ_{q<M} h_q + Σ_{m<K−M} (b_m + 1) for every 0 ≤ M ≤ K.
   - Proof: q(Y) = Σ_{Q,U} (det R[U^c, Q])²·Π_{q∈Q}Λ_q·Π_{m∈U}(Y·β̂_m) (Cauchy–Binet on [R | I]·diag(Λ, YB̂)·[R | I]ᵀ).
   - h and b are strictly increasing: h by at least 4 per step, b by at least 2. So the minimal term is unique, namely Q = [0, M), U = [0, K−M).
   - Its coefficient det R[[K−M, K), [0, M)]² is a unit.
   - (The earlier "violations" for M ≳ K/4 were truncation artefacts of C.)
2. **Content law with the tie rule** (section 10 is now unconditional): e₂(K) = −min_M T(K, M) − [minimum attained twice]. Lemma C is proved.
3. **det law.** v₂(det A_K) = min_M T(K, M) when the minimiser is unique, and ≥ min + 1 at a tie.
4. **Smith form (v-scale).** SNF(A^v_K(0)) = the K smallest elements of the multiset {b_m} ∪ {h_q}, i.e. the sorted union {b_m : m < K − M*} ∪ {h_q : q < M*}, with the top divisor raised at a tie.
   - The gcds of the k×k minors are ≥ the sum of the k smallest elements (Cauchy–Binet, R integral).
   - Equality is realised by the leading principal minor [0, k), whose unique minimal term has a Pascal-corner unit coefficient.
   - At an internal boundary tie between h_M and b_{k−M−1}, use I = J = [0, k−M−1) ∪ [k−M, k+1) instead. That excludes the tied partner, and the survivor's coefficient is det[C(k−M+i, j)]_{i,j≤M} = 1.
   - At k = K there is only one minor, hence the raised top divisor.
   - Lemma B is proved: the "2-adic orthogonality" is the Pascal matrix mod 2.

What this means:
- **The band is Krawtchouk:** the X-coefficient of the machine is a binomial measure, and 1 + s₂(K−1−m) − s₂(m) is v₂ of its (even or odd) Krawtchouk norms.
- **The moment species are CDH:** they are the moment pivots of the machine at X = C, the 2-adic Catalan constant, where the pole data are the moments' own Stieltjes transform.
- **The split at M* ≈ K/6** is just "take the K smallest of two increasing lists".
- **The "partner" and "pole band" pictures** in the u-frame and the Newton frame are shadows of this pencil.

**Numerical confirmation of the whole chain:**
- th_catalan_struct.py (Catalan K = 6, 8, 12, 16, 20, 24) and th_twist_struct.py (twisted K = 6, 8, 11, 12, 16), including the ties at Catalan K = 12 and twisted K = 11.
- Every check passes: A(C) = A(0) + CB to the precision of C, C via moments, SNF(A(C)) = h, the mod-2 structures, and the coefficient law for all M.
- Copies are kept in `proof_checks/`.
