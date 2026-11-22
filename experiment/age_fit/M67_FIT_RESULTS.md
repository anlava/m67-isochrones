# M67 Turnoff Fit Results

**Date:** 2026-11-22
**Data:** 257 M67 turnoff stars (11.5 ≤ G < 13.5) from Hunt & Reffert (2023)
**Colour cross-check:** GaiaDR3 IDs and BP−RP from TSV vs results/m67_members_selected.csv — 257/257 match
**Models:** PARSEC 1.2S + COLIBRI, Gaia EDR3 photometry, Cardelli Rv=3.1 extinction
**Metric:** median Euclidean distance to isochrone polyline segment in (G, BP−RP)
**Phases:** label ∈ {0,1,2} (PMS, MS, subgiants)

## Best-fit ages

| Grid | Parameters | Best logAge | Age (Gyr) | M (mag) |
|------|-----------|-------------|-----------|---------|
| Cavallo | [Fe/H]=+0.32, Av=0.18, dMod=9.64 | **9.35** | **2.24** | 0.02686 |
| Bossini | [Fe/H]=0.00, Av=0.115, dMod=9.726 | **9.55** | **3.55** | 0.02682 |

## Full metric table

| logAge | Cavallo M (mag) | Bossini M (mag) |
|--------|----------------|-----------------|
| 9.3 | 0.03602 | 0.12876 |
| 9.35 | 0.02686 | 0.10937 |
| 9.4 | 0.03518 | 0.08888 |
| 9.45 | 0.06193 | 0.05991 |
| 9.5 | 0.08891 | 0.03966 |
| 9.55 | 0.11151 | 0.02682 |
| 9.6 | 0.13642 | 0.04823 |
| 9.65 | 0.18285 | 0.07116 |
| 9.7 | 0.19252 | 0.08049 |

## Notes

- Cavallo minimum at logAge=9.35 (2.24 Gyr); Bossini at logAge=9.55 (3.55 Gyr).
- Both reach nearly the same M ≈ 0.027 mag — the median star-to-isochrone distance in CMD.
- The 1.3 Gyr age difference arises from different [Fe/H], Av, and distance modulus:
  Cavallo super-solar ([Fe/H]=+0.32) shifts the turnoff to bluer colours at fixed age,
  so a younger age is needed to match the observed colour of M67.
- Bossini (solar metallicity) gives 3.55 Gyr, consistent with literature estimates (3.5–4.5 Gyr).
- M ≈ 0.027 mag (~0.03 mag) corresponds to the typical Gaia DR3 photometric uncertainty
  in this brightness range (σ_G ≈ 0.005, σ_BP−RP ≈ 0.01 mag), confirming the fit is
  limited by model systematics, not by measurement noise.
