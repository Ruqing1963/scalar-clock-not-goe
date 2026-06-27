# scalar-clock-not-goe

**Narrow scalar clocks, not level repulsion: a fourth spurious-GOE hazard and a re-audit of metallogenic timing statistics**

Ruqing Chen — GUT Geoservice Inc., Montréal, Québec, Canada

A short audit note. It adds a **fourth** hazard to the three in the companion note
[spurious-level-repulsion](https://doi.org/10.5281/zenodo.20883107) (declustering,
seasonal gating, Wishart null), and uses it to re-audit a published claim of
**temporal GOE/GUE level repulsion in ore-forming systems**.

## The fourth hazard

A charge-and-release reservoir is a **relaxation oscillator = renewal process**:
its inter-event intervals are independent, so it has no spectral rigidity and
**cannot** produce genuine level repulsion. Yet a renewal clock with a narrow
interval distribution reproduces **any** value of the spacing ratio ⟨r⟩:

- ⟨r⟩ sweeps **0.39 → 0.53 (GOE) → 0.60 (GUE) → 0.88** purely as a function of the
  interval coefficient of variation (CV). **⟨r⟩ alone is uninformative about repulsion.**
- The genuine GOE fingerprints a scalar clock lacks — a negative shuffle deficit
  (⟨r⟩_obs < ⟨r⟩_shuf, rigidity) and logarithmic number variance — need **N ≳ 300**
  intervals to detect; geochronologic sequences have tens.

## Re-audit of orogenic-gold timing (metallogeny-rmt)

Reproducing the published ⟨r⟩ = 0.678 (n = 22) exactly, we find:

1. **Wrong quantity** — it measures within-deposit scatter among geochronometers
   (U-Pb / Re-Os / Ar-Ar on the same deposit), not province-scale pulse timing.
2. **Scalar clock** — CV = 0.31; a zero-repulsion renewal clock of CV = 0.38
   reproduces ⟨r⟩ = 0.68. No repulsion needed.
3. **Rigidity is an artifact** — the shuffle deficit −0.031 is reproduced by a
   zero-repulsion null through the same within-deposit-normalise-then-pool
   construction (it is a normalisation constraint, not spectral rigidity).
4. **Unresolvable** — jittering ages within ±4–6 Myr moves ⟨r⟩ across the entire
   Poisson-to-GUE range [0.39, 0.67]; the reported value sits at the top edge.

## Restated organising principle

| Configuration | Signature | Interpretation |
|---|---|---|
| Single relaxation oscillator | high ⟨r⟩, CV-set, shuffle-invariant | **scalar renewal clock (not GOE)** |
| Chaotic / strongly-mixed spectrum | ⟨r⟩≈0.53, ⟨r⟩_obs<⟨r⟩_shuf, rigid | GOE level repulsion |
| Superposition of many sources | ⟨r⟩→0.386 | Poisson |

Genuine Wigner–Dyson statistics require a chaotic, strongly-mixed spectrum (BGS);
a fluid reservoir does not possess one. The **superposition→Poisson** half of the
principle is unaffected.

## Implications (basis for two v2 errata)

- **metallogeny-rmt** (10.5281/zenodo.20768849): the "level repulsion" reading is
  not supported; reframe the single-source signal as a narrow scalar clock. The
  superposition→Poisson contrast and the data compilation survive.
- **sigma-spatial-rmt** (10.5281/zenodo.20778697): the spatial clustering result
  and the declustering-fabricates-GOE warning are **sound and unaffected**; only
  the "time–space asymmetry" framing is restated as "temporal scalar clock +
  spatial clustering".

## Reproduce

```bash
pip install -r requirements.txt
cd code && python scalar_clock_audit.py --outdir ..
```

Regenerates both figures and all reported numbers. Self-contained.

## Layout

```
code/scalar_clock_audit.py     mechanism demo + real-data re-audit + figures
data/orogenic_gold_ages.csv    ages re-analysed (from metallogeny-rmt; see data/README.md)
figures/fig_mechanism.{pdf,png}   Fig 1
        fig_realdata.{pdf,png}    Fig 2
paper/scalar_clock_audit.{pdf,tex}
```

## Citation
Cite the paper (`paper/scalar_clock_audit.pdf`) and this repository (Zenodo DOI on
release).

## License
CC BY 4.0 (see `LICENSE`).
