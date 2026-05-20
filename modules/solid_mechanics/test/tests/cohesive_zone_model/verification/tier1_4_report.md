# Tier 1.4 — Mixed-mode at β = 1: Verification Report

## Summary

`BiLinearMixedModeTraction` was verified under symmetric mixed-mode loading
(opening rate = shearing rate, hence mode-mix β = 1). After a single
remediation of the bulk material stiffness, the simulation reproduces the
closed-form Camanho–Dávila mixed-mode bilinear law to roughly **4×10⁻⁸**
in absolute traction (~10⁻⁹ relative to the peak), with the mode-mix
constraint β = 1 holding to **4×10⁻⁹**.

| Quantity | Tolerance | Achieved | Margin |
|---|---|---|---|
| `max \|tm − tm_ref\|` | 1e-6 | 4.05e-08 | 25× |
| `\|β − 1\| max` | 1e-6 | 4.33e-09 | 230× |
| `\|tx\| − \|ty\|` (component symmetry) | 1e-6 | 8.95e-08 | 11× |

Status: **PASS**.

## Test setup

- Mesh / bulk / solver: identical to Tier 1.1.
- Boundary conditions: top face pulled by the *same* function in both `disp_x`
  and `disp_y`, so on the cohesive interface `δ_x = δ_y` and therefore
  `β = δ_shear / δ_normal = 1`.
- Cohesive law: `BiLinearMixedModeTraction` with
  `K = 5000, N = 50, S = 30, GI_c = GII_c = 0.5, η = 1.45`.
- Loading function: piecewise linear with breakpoint at the per-axis
  mixed-mode onset displacement
  `u₀ = δ_m^0 / √2 ≈ 5.145·10⁻³`, end value
  `u_max = 1.5 · δ_m^f / √2 ≈ 2.916·10⁻²`.

The closed-form mixed-mode landmarks (B-K reduces to mode-I value because
`GI_c = GII_c` for this test) are:

```
δ_m^0 = δ_n^0 · δ_s^0 · √[(1 + β²) / (δ_s^0² + (β·δ_n^0)²)]
      = 0.01 · 0.006 · √[2 / (0.006² + 0.01²)]
      = 7.276·10⁻³

δ_m^f = (2 / (K · δ_m^0)) · [GI_c + (GII_c − GI_c)(β²/(1+β²))^η]
      = 2 · 0.5 / (5000 · 7.276·10⁻³)
      = 2.749·10⁻²

t_m^peak = K · δ_m^0 = 36.38
```

## What happened on the first run (failure)

With `E_bulk = 10⁹` (the value used successfully in Tier 1.1, 1.2, 1.3, and
1.5), the verifier reported:

```
beta observed range   : [0.999965, 1.000000]   max-err 3.48e-05
max |tx| - |ty|       : 8.95e-04
max |tm - tm_ref|     : 4.05e-04   at the elastic-softening kink
```

All three "errors" were ~4 orders of magnitude *above* the 10⁻⁶ tolerance,
yet the bilinear curve *visually* overlaid analytical perfectly. The
softening branch had a uniform ~10⁻⁴ offset that disappeared in both the
elastic regime and the post-decohesion regime.

This pattern — clean elastic, clean post-failure, biased in softening — is
the signature of a *kinematic* drift rather than a constitutive bug.

### Root cause

The pulled `1×2` rectangular mesh introduces an asymmetry: the bulk's
elastic response to imposed normal displacement differs from its response
to imposed shear displacement. Per unit interface area:

- Normal stiffness of one bulk block: `E/h = 10⁹ / 1 = 10⁹`
- Shear stiffness of one bulk block:  `G/h = E / (2(1+ν)) / h ≈ 3.85·10⁸`

The bulk is **2.6× stiffer in the normal direction than in shear**.

For the cohesive interface (penalty stiffness `K = 5·10³` in both channels),
the bulk-to-interface compliance ratios in the elastic regime are

- Normal:  `K / (E/h) = 5·10³ / 10⁹ = 5·10⁻⁶`
- Shear:   `K / (G/h) = 5·10³ / 3.85·10⁸ ≈ 1.3·10⁻⁵`

So when the top face is pulled identically in `x` and `y`, the bulk swallows
a slightly *larger* fraction of the imposed `disp_x` than of `disp_y`,
leaving slightly *less* tangential jump than normal jump on the cohesive
interface. The jumps differ by `~5·10⁻⁶` to `~10⁻⁵` relative — exactly the
~3·10⁻⁵ β-drift we measured.

That tiny β-drift then propagates into the law:

- `δ_m^0` (per the Camanho–Dávila formula) shifts slightly because it
  depends on β.
- Damage `d(δ_m_max)` is computed from the simulation's actual
  (slightly drifted) `δ_m_max`, while the verifier's closed-form
  reference uses the *exact* β = 1 landmarks.
- The two trajectories then differ by `O(β-drift)` — a small but
  systematic bias along the entire softening branch (~10⁻⁴ traction).

The error is *real* in the sense that the simulation's β was not exactly 1,
but it is *not* a CZM-law bug — the cohesive law is doing the right thing
for the slightly-drifted β it actually sees.

### The fix

Increased bulk Young's modulus from `E = 10⁹` to `E = 10¹³`. This
multiplies both bulk stiffnesses by `10⁴`, dropping the bulk-to-interface
compliance ratios to:

- Normal: `5·10⁻¹⁰`
- Shear:  `1.3·10⁻⁹`

Now the bulk is essentially a rigid backing. Both jumps equal the imposed
displacement to ~9 significant figures, and β drift collapses by 4 orders
of magnitude. Re-running:

```
beta observed range   : [1.000000, 1.000000]   max-err 4.33e-09   (was 3.48e-05)
max |tx| - |ty|       : 8.95e-08             (was 8.95e-04)
max |tm - tm_ref|     : 4.05e-08             (was 4.05e-04)
```

A clean 4-orders-of-magnitude improvement across all three metrics —
exactly what the bulk-stiffness scaling argument predicts.

The change is recorded in the input file with a comment noting why the
modulus was bumped, since the value is otherwise unphysical (most users
running 1.1 / 1.3 / 1.5 will not need it).

## Why Tier 1.1, 1.2, 1.3, 1.5 did not see this

In Tier 1.1, 1.3, 1.5 the loading is *pure* normal — only `disp_y` is driven,
`disp_x` is not. There is only one channel, so there is no opportunity for
asymmetric bulk swallowing.

In Tier 1.2 the loading is pure shear — `disp_x` is driven, `disp_y` is held
at zero on both faces. Again only one channel, no asymmetry.

The β = 1 condition is fragile in a way the single-channel tests never
exercise.

## Lessons

- For *single-mode* CZM verification, `E_bulk = 10⁹` is plenty stiff with
  `K = 5·10³` (compliance ratio ~5·10⁻⁶).
- For *mixed-mode* CZM verification at a fixed β, the bulk must be
  significantly stiffer in *all* relevant channels — or the analytical
  reference must be computed using the actual (drifted) β observed in the
  simulation.
- The first-failure trace (clean elastic + clean post-failure + biased
  softening) is diagnostic for kinematic/coupling errors versus
  constitutive errors. Constitutive bugs typically produce errors that
  scale with traction magnitude across all regimes; kinematic drift
  shows up only where the law's mode-mixity calculation is active.

## Files

```
verification/
  tier1_4_mixed_mode_beta1.i        — input (E_bulk = 1e13)
  verify_tier1_4.py                 — analytical comparator + plot
  tier1_4_mixed_mode_beta1_out.csv  — simulation CSV (201 rows)
  tier1_4_mixed_mode_beta1_compare.png — bilinear envelope vs sim, log-error
```
