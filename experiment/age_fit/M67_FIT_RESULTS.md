# M67 Turnoff Fit Results

**Date:** 2026-11-22  
**Data:** 257 M67 turnoff stars (11.5 ≤ G < 13.5) from Hunt & Reffert (2023)  
**Models:** PARSEC 1.2S + COLIBRI, Gaia EDR3 photometry, Cardelli Rv=3.1 extinction  
**Metric:** median Euclidean distance in (G, BP−RP) space, projected to isochrone segments  
**Phases:** label ∈ {0,1,2} (PMS, MS, subgiants)

## Best-fit ages

| Grid | Parameters | Best logAge | Age (Gyr) | M (mag) | Boundary? |
|------|-----------|-------------|-----------|---------|-----------|
| Cavallo 2024 | [Fe/H]=+0.32, Av=0.18, dMod=9.64 | **9.70** | **5.01** | 11.2164 | ⚠ minimum at oldest grid node |
| Bossini 2019 | [Fe/H]=0.00, Av=0.115, dMod=9.726 | **9.70** | **5.01** | 11.3220 | ⚠ minimum at oldest grid node |

## Full metric table

| logAge | Cavallo M (mag) | Bossini M (mag) |
|--------|----------------|-----------------|
| 9.30   | 11.3273        | 11.4193         |
| 9.35   | 11.3099        | 11.4044         |
| 9.40   | 11.2881        | 11.3839         |
| 9.45   | 11.2736        | 11.3727         |
| 9.50   | 11.2564        | 11.3541         |
| 9.55   | 11.2456        | 11.3379         |
| 9.60   | 11.2346        | 11.3297         |
| 9.65   | 11.2234        | 11.3266         |
| 9.70   | **11.2164**    | **11.3220**     |

## Notes

- Both grids show a monotonic decrease in median distance from logAge=9.30 to 9.70.  
  The true minimum lies at logAge > 9.70 (age >5 Gyr), beyond the grid range.
- Cavallo fits systematically better than Bossini (ΔM ≈ 0.09–0.10 mag across all ages).
- The known literature age of M67 is ~3.5–4.5 Gyr (logAge ≈ 9.55–9.65).  
  Our fit is systematically older, possibly due to:
  - Different metallicity/extinction assumptions
  - Systematic offsets between PARSEC isochrones and Gaia photometry
  - The limited G-range (11.5–13.5) not including the full turnoff morphology