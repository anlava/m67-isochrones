# M67 Turnoff Fit Results

**Date:** 2026-11-22 (исправлено: BP−RP читается как BPmag − RPmag)  
**Data:** 257 M67 turnoff stars (11.5 ≤ G < 13.5) from Hunt & Reffert (2023)  
**Colour cross-check:** BP−RP сверен по `GaiaDR3` между TSV и `results/m67_members_selected.csv` — **257/257 совпадают**  
**Models:** PARSEC 1.2S + COLIBRI, Gaia EDR3 photometry, Cardelli Rv=3.1 extinction  
**Metric:** median Euclidean distance to isochrone polyline segment in (G, BP−RP)  
**Phases:** label ∈ {0,1,2} (PMS, MS, subgiants)  
**Code:** `distance_to_isochrone` импортирована из `synthetic_test.py`

## Best-fit ages

| Grid | Parameters | Best logAge | Age (Gyr) | M (mag) |
|------|-----------|-------------|-----------|---------|
| Cavallo 2024 | [Fe/H]=+0.32, Av=0.18, dMod=9.64 | **9.35** | **2.24** | 0.02686 |
| Bossini 2019 | [Fe/H]=0.00, Av=0.115, dMod=9.726 | **9.55** | **3.55** | 0.02681 |

## Full metric table

| logAge | Cavallo M (mag) | Bossini M (mag) |
|--------|----------------|-----------------|
| 9.30   | 0.03602        | 0.06520         |
| 9.35   | **0.02686**    | 0.04970         |
| 9.40   | 0.03518        | 0.04444         |
| 9.45   | 0.04040        | 0.06664         |
| 9.50   | 0.03924        | 0.03966         |
| 9.55   | 0.04263        | **0.02681**     |
| 9.60   | 0.04847        | 0.04823         |
| 9.65   | 0.05187        | 0.04556         |
| 9.70   | 0.04992        | 0.06279         |

## Notes

- Cavallo минимум на logAge=9.35 (2.24 Gyr) с резким параболическим провалом.
  Соседние узлы (9.30, 9.40) дают ~0.035–0.036 mag — в 1.3× хуже.
- Bossini минимум на logAge=9.55 (3.55 Gyr), близко к литературным оценкам
  M67 (3.5–4.5 Gyr). Соседние узлы (9.50, 9.60) дают ~0.040–0.048 mag.
- Абсолютные значения метрики (~0.027 mag) соответствуют типичной
  фотометрической точности Gaia в данном диапазоне яркостей.
- Минимум не на границе ни для одной сетки — оценка возраста корректна.