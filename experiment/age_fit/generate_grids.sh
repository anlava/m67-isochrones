#!/usr/bin/env bash
# generate_grids.sh
# Генерация двух сеток изохрон PARSEC 1.2S + COLIBRI через CMD 3.9
# для подбора возраста M67 по области поворота.
#
# Сетка A (Cavallo 2024):   [Fe/H]=0.32, Av=0.18, dMod=9.64
# Сетка B (Bossini 2019):   [Fe/H]=0.00, Av=0.115, dMod=9.726
# logAge: 9.30–9.70, шаг 0.05 dex (всего 9 значений)
#
# Зависимости: curl, python3, bash
# Использование: bash experiment/age_fit/generate_grids.sh
# Результат: *.dat в data/raw/, *.csv в results/

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUTDIR="$ROOT/data/raw"
RESDIR="$ROOT/results"
CMD_URL="https://stev.oapd.inaf.it/cgi-bin/cmd_3.9"
TMP_BASE="https://stev.oapd.inaf.it"

# Параметры сеток
declare -A GRIDS
GRIDS["cavallo2024_av018"]='Z=0.029013|Av=0.18'
GRIDS["bossini2019_av0115"]='Z=0.01471|Av=0.115'

# Общие параметры CMD 3.9 (PARSEC 1.2S, Gaia EDR3, Cardelli Rv=3.1)
COMMON=(
  cmd_version=3.9
  submit_form=Submit
  track_parsec=parsec_CAF09_v1.2S
  track_colibri=parsec_CAF09_v1.2S_S_LMC_08_web
  photsys_file=YBC_tab_mag_odfnew/tab_mag_gaiaEDR3.dat
  photsys_version=YBCnewVega
  extinction_coeff=constant
  extinction_curve=cardelli
  isoc_isagelog=1
  isoc_ismetlog=0
  output_kind=0
  dust_sourceM=nodustM
  dust_sourceC=nodustC
  kind_LPV=4
  n_inTPC=10
  eta_reimers=0.2
  kind_interp=1
  kind_postagb=-1
  track_postagb=no
  imf_file=tab_imf/imf_kroupa_orig.dat
)

# .cgifields (скалярные поля — радио-кнопки и чекбоксы)
CGFIELDS=(
  dust_sourceC dust_sourceM extinction_coeff extinction_curve
  isoc_isagelog isoc_ismetlog kind_LPV output_gzip output_kind
  photsys_version track_colibri track_parsec track_postagb
)

# Диапазон возрастов: logAge 9.30 – 9.70, шаг 0.05
AGES=()
for i in $(seq 0 8); do
  AGES+=($(python3 -c "print(f'{9.30 + $i * 0.05:.2f}')"))
done

echo "CMD 3.9: $CMD_URL"
echo "Диапазон logAge: ${AGES[0]}–${AGES[-1]} (${#AGES[@]} значений)"
echo "Сетки: ${!GRIDS[*]}"
echo ""

mkdir -p "$OUTDIR" "$RESDIR"

for GRID_NAME in "${!GRIDS[@]}"; do
  IFS='|' read -r ZSTR AVSTR <<< "${GRIDS[$GRID_NAME]}"
  ZVAL="${ZSTR#Z=}"
  AVVAL="${AVSTR#Av=}"
  echo "--- Сетка: $GRID_NAME (Z=$ZVAL, Av=$AVVAL) ---"

  # Собираем все .dat файлы во временный
  TMP_DAT=$(mktemp)
  DATA_COUNT=0

  for LOGAGE in "${AGES[@]}"; do
    echo -n "  logAge=$LOGAGE ... "

    # Вызов CMD 3.9
    CURL_ARGS=()
    for p in "${COMMON[@]}"; do
      CURL_ARGS+=(--data-urlencode "$p")
    done
    CURL_ARGS+=(--data-urlencode "extinction_av=$AVVAL")
    CURL_ARGS+=(--data-urlencode "isoc_lagelow=$LOGAGE")
    CURL_ARGS+=(--data-urlencode "isoc_lageupp=$LOGAGE")
    CURL_ARGS+=(--data-urlencode "isoc_dlage=0.0")
    CURL_ARGS+=(--data-urlencode "isoc_zlow=$ZVAL")
    CURL_ARGS+=(--data-urlencode "isoc_zupp=$ZVAL")
    CURL_ARGS+=(--data-urlencode "isoc_dz=0.0")
    for f in "${CGFIELDS[@]}"; do
      CURL_ARGS+=(--data-urlencode ".cgifields=$f")
    done

    RESULT_PAGE=$(curl -sk -X POST "${CURL_ARGS[@]}" "$CMD_URL" 2>/dev/null)
    DAT_FILE=$(echo "$RESULT_PAGE" | grep -o 'tmp/output[0-9]*\.dat' | head -1)

    if [ -z "$DAT_FILE" ]; then
      echo "ОШИБКА: не удалось получить ссылку на результат"
      rm -f "$TMP_DAT"
      exit 1
    fi

    # Ждём готовности и скачиваем
    sleep 2
    if ! curl -sk "${TMP_BASE}/${DAT_FILE}" >> "$TMP_DAT" 2>/dev/null; then
      echo "ОШИБКА: не удалось скачать ${DAT_FILE}"
      rm -f "$TMP_DAT"
      exit 1
    fi

    # Считаем строки данных (не комментарии)
    NLINES=$(grep -c -v '^#' "${TMP_DAT}" 2>/dev/null || echo 0)
    echo "OK ($NLINES строк)"
    DATA_COUNT=$((DATA_COUNT + NLINES))
  done

  # Сохраняем объединённый .dat
  DAT_OUT="${OUTDIR}/m67_parsec12s_${GRID_NAME}.dat"
  mv "$TMP_DAT" "$DAT_OUT"
  echo "  → $DAT_OUT (всего $DATA_COUNT строк данных)"

  # Конвертируем в CSV (G_app = Gmag + dMod, BP_RP = G_BP - G_RP)
  CSV_OUT="${RESDIR}/m67_isochrone_${GRID_NAME}.csv"
  python3 -c "
import csv, sys
with open('$DAT_OUT') as f:
    lines = [l for l in f if not l.startswith('#') and l.strip()]

# Определяем dMod для сетки
dMod = 9.64 if '$GRID_NAME' == 'cavallo2024_av018' else 9.726

with open('$CSV_OUT', 'w', newline='') as out:
    w = csv.writer(out)
    w.writerow(['Gmag','G_BPmag','G_RPmag','G_app','BP_RP','Mini','logAge','MH','label'])
    for line in lines:
        cols = line.split()
        if len(cols) < 15:
            continue
        try:
            zini = float(cols[0])
            mh = float(cols[1])
            logage = float(cols[2])
            mini = float(cols[3])
            label = int(cols[9])
            gmag = float(cols[30])
            g_bpmag = float(cols[31])
            g_rpmag = float(cols[32])
        except (ValueError, IndexError):
            continue
        g_app = gmag + dMod
        bp_rp = g_bpmag - g_rpmag
        w.writerow([f'{gmag:.6f}', f'{g_bpmag:.6f}', f'{g_rpmag:.6f}',
                     f'{g_app:.6f}', f'{bp_rp:.6f}',
                     f'{mini:.10f}', f'{logage:.5f}', f'{mh:.5f}', label])
  "
  echo "  → $CSV_OUT"
  echo ""
done

echo "Готово. Сетки сохранены в:"
ls -lh "$OUTDIR"/m67_parsec12s_*.dat
echo ""
echo "CSV:"
ls -lh "$RESDIR"/m67_isochrone_*.csv