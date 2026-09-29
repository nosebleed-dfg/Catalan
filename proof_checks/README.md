# proof_checks: numerical checks of each link in twoadic_law_proof.md, sections 8–12 (2026-09-28)
- Run with Python 3.12: C:\Users\PC\AppData\Local\Programs\Python\Python312\python.exe.
- Needs python-flint, sympy and the project's twisted_hankel.py; lemlib.py adds the parent folder to the path.

| script | what it checks | proof section |
|---|---|---|
| lemlib.py | shared helpers:<br>• 2-adic C (Mahler series)<br>• pole measures<br>• machine parts in the v-scale<br>• monic orthogonal polynomials, Smith form, band, h_q | — |
| th_lemmaA2.py `K1,..` | norms of the discrete pole measure vs κ(K, n): Catalan matches for all n; twisted for n < K/2 | 9 |
| th_lemmaA5.py `K1,..` | Catalan perturbation margins in the Krawtchouk basis (minimum exactly 1) | 9 |
| th_lemmaA7.py `K1,..` | each link of the Lemma A proof:<br>• mirror identity<br>• norm formula<br>• λ bound<br>• Mahler bound<br>• margins | 9 |
| th_lemmaA4.py / th_lemmaA6.py `K1,..` | Krawtchouk model + moment corner reproduces SNF(A) (twisted needs μ(−1) = 4C) | 11 |
| th_lemmaB1–B3.py | the moment corner in the Krawtchouk basis: entries have v₂ = −1, plus top-block pivots (superseded by 12) | 11 |
| th_lemmaC1.py `mach K1,..` | coefficient law near the minimum and the content derivation (tie rule) | 10 |
| th_stieltjes.py | Catalan: F_j(C) = −Σ μ_k w_j^{−k−1} (2-adic Stieltjes identity) | 12, Thm 6a |
| th_stieltjes_tw.py | twisted: −G₁(a) = Σ (−1)^k m_k (1+4a²)^{−k−1} | 12, Thm 6b |
| th_catalan_struct.py `K [KM]` | the Catalan chain end to end:<br>• A(C) = A(0) + CB<br>• C via moments<br>• SNF(A(C)) = h<br>• mod-2 structures, T ≡ Pascal<br>• coefficient law for every M | 12 |
| th_twist_struct.py `K [KM]` | the same chain for the twisted machine | 12 |
| th_thm7_rigorous.py `catalan\|twist K [KM]` | Theorem 7 without square roots (2026-09-29):<br>• rescaled monic Jacobi matrix integral, divisible by ρ<br>• G = h_p ρ^{q−p} Π̃<br>• elimination in ℤ₂<br>• integral L with v₂(L_pq) ≥ 3(p−q)+2 (resp. +1)<br>• G = LΛLᵀ, v₂(Λ_p) = v₂(h_p) | 12, Lemma 12.1 / Thm 7 |
| th_lemmaA_balanced.py `K1,..` | section 9 without square roots (2026-09-29):<br>• multiplication by x on the nodes = monic Jacobi matrix<br>• balanced-valuation bounds on C(𝒥, k)<br>• M^S = Σ s_k C(𝒥, k)<br>• margins ≥ 1 | 9, step 4 |
| th_paper_checks.py | the paper's supporting checks (2026-09-29):<br>• Lemma 3′ (the finite-difference identity for g, and the sharp bound v₂(Δ^k g) ≥ k + v₂((k+1)!))<br>• Appendix A integrals (40 digits)<br>• the continuous dual Hahn recurrence and norms against the exact moments (n < 30)<br>• the tangent-number formula for C<br>• the Theorem 3 growth table | paper §4, §11, App. A–B; proof §4 |
| th_kraw_recurrence.py `K1,..` | Theorem 4's explicit integer recurrences for the v-scale Krawtchouk polynomials, and the norm valuations | 8, Thm 4 |

Run log (2026-09-28): th_catalan_struct.py K = 6, 8, 12, 16, 20, 24 and th_twist_struct.py K = 6, 8, 11, 12, 16. All checks pass, including the coefficient law for every M and the ties at Catalan K = 12 and twisted K = 11. The default depth KM (6K+60 Catalan, 10K+60 twisted) is enough: only relative precision matters, i.e. KM > max h_q.
