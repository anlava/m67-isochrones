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
OUTDIR="$ROOT/experiment/age_fit"
RESDIR="$ROOT/experiment/age_fit"
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
  DAT_OUT="${OUTDIR}/grid_${GRID_NAME}.dat"
  mv "$TMP_DAT" "$DAT_OUT"
  echo "  → $DAT_OUT (всего $DATA_COUNT строк данных)"

  # Конвертируем в CSV (чтение колонок по заголовку)
  CSV_OUT="${RESDIR}/grid_${GRID_NAME}.csv"
  python3 -c "
import csv
with open('$DAT_OUT') as f:
    lines = f.readlines()

# Парсим заголовок колонок
header_line = None
for line in lines:
    if line.startswith('#') and 'Zini' in line and 'Gmag' in line:
        header_line = line.lstrip('#').strip()
        break
if header_line is None:
    print('ОШИБКА: заголовок колонок не найден', file=sys.stderr)
    exit(1)

col_names = header_line.split()
required = ['Gmag','G_BPmag','G_RPmag','Mini','logAge','MH','label']
indices = {name: col_names.index(name) for name in required}

# Определяем dMod для сетки
dMod = 9.64 if '$GRID_NAME' == 'cavallo2024_av018' else 9.726

data_lines = [l for l in lines if not l.startswith('#') and l.strip()]
with open('$CSV_OUT', 'w', newline='') as out:
    w = csv.DictWriter(out, fieldnames=['Gmag','G_BPmag','G_RPmag','G_app','BP_RP','Mini','logAge','MH','label'])
    w.writeheader()
    for line in data_lines:
        cols = line.split()
        if len(cols) != len(col_names):
            continue
        try:
            gmag = float(cols[indices['Gmag']])
            g_bp = float(cols[indices['G_BPmag']])
            g_rp = float(cols[indices['G_RPmag']])
        except (ValueError, IndexError):
            continue
        g_app = gmag + dMod
        bp_rp = g_bp - g_rp
        w.writerow({
            'Gmag': f'{gmag:.6f}', 'G_BPmag': f'{g_bp:.6f}', 'G_RPmag': f'{g_rp:.6f}',
            'G_app': f'{g_app:.6f}', 'BP_RP': f'{bp_rp:.6f}',
            'Mini': cols[indices['Mini']], 'logAge': cols[indices['logAge']],
            'MH': cols[indices['MH']], 'label': cols[indices['label']],
        })
print(f'→ $CSV_OUT')
" || exit 1
  echo ""
  echo ""
done

echo "Готово. Сетки сохранены в:"
ls -lh "$OUTDIR"/grid_*.dat
echo ""
echo "CSV:"
ls -lh "$RESDIR"/grid_*.csv