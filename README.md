# Ropelength and crossing number: computer check

This program checks the rational inequalities in the paper
"A lower bound for ropelength in terms of crossing number" by Arseniy Varlamov.
The paper proves L ≥ 2.124 C^(3/4) for every knot and link, where L is the ropelength (thickness measured as a radius) and C is the crossing number.

## Run

```
python3 -B verify.py
```

Python 3 standard library only; every decision uses exact fractions. The run takes well under a second.
It prints `PASS` and writes the results to `certificate.json`.

## What is not checked

The geometric lemmas and the integration argument are proved in the paper, not here.
Section 6 of the paper lists exactly what the program checks.
