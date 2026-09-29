# Catalan

Research notebook on Catalan's constant G = Σ (−1)^k/(2k+1)² and related L-values and odd zeta values. It studies Hankel-determinant constructions and reads them through carries, primes and base variance. Catalan's constant is the proof of concept; the open irrationality problems are the goal.

**Status:** working research notes and code, not refereed.

## Where to start
- `HANDOFF_catalan_hankel-1.md`: the entry point. It gives the state of every route, the commands to reproduce results, and what is closed or open.
- `dedekind_carry_results-1.md`: the dated results log, in order.
- `twoadic_law_proof.md`: the 2-adic content law of the lattice Catalan Hankel machines (the "pencil theorem"):
  - the machine is A(X) = A(C) + (X − C)·B, with C a 2-adic Catalan constant;
  - B is a Krawtchouk (binomial) form, and A(C) is the continuous dual Hahn moment form;
  - mod 2 their bases differ by the Pascal matrix, and Cauchy–Binet gives the Smith form and the content law.
  - `proof_checks/` has one numerical check script per proof step (see its README).
- `refined_lemma_proof.md`: the refined denominator lemma for the multi-β forms.

## Running the code
Python 3.12 with `python-flint` (0.9), `sympy`, `mpmath` and `numpy`. For example:

```
python twoadic_law.py law all
python proof_checks/th_catalan_struct.py 12
```

Scripts are referred to by name from the handoff and the results log. The `.json` and `.txt` files are cached outputs.
