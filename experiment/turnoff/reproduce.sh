#!/usr/bin/env bash
#
# reproduce.sh — воспроизвести CMD области поворота и отбор субгигантов
#                из исходных данных (data/raw/m67_members_hr2023.tsv).
#
# Все команды из корня проекта.
# Использует python из активного окружения (первый python3 в PATH).

set -euo pipefail

# Автоопределение Python с matplotlib + numpy
for _py in python3 python3.9 python3.10 python3.11 python3.12; do
    _p=$(command -v $_py 2>/dev/null || true)
    if [ -n "$_p" ] && "$_p" -c "import matplotlib,numpy" 2>/dev/null; then
        PYTHON="$_p"
        break
    fi
done
if [ -z "${PYTHON:-}" ]; then
    echo "ОШИБКА: не найден Python с matplotlib и numpy" >&2
    exit 1
fi
cd "$(git rev-parse --show-toplevel 2>/dev/null || echo "$(dirname "$0")/../..")"

echo "Python: $($PYTHON --version)"
echo "Корень: $(pwd)"

echo ""
echo "=== 1. Убеждаемся, что results/m67_members_selected.csv существует ==="
if [ ! -f results/m67_members_selected.csv ]; then
    echo "Файл не найден. Запускаю src/plot_m67.py..."
    $PYTHON src/plot_m67.py
else
    echo "  OK: results/m67_members_selected.csv ($(wc -l < results/m67_members_selected.csv) строк)"
fi

echo ""
echo "=== 2. CMD + отбор кандидатов (правило S1) ==="
$PYTHON experiment/turnoff/plot_turnoff_cmd.py

echo ""
echo "=== 3. Проверка 1: восстановление ID из TSV ==="
$PYTHON experiment/turnoff/check_recover_ids.py

echo ""
echo "=== 4. Проверка 2: синтетический остаток (корректный знак) ==="
$PYTHON experiment/turnoff/check_synthetic_residual.py

echo ""
echo "=== 5. Проверка 2: синтетический остаток (обращённый знак, ожидаемый сбой) ==="
if $PYTHON experiment/turnoff/check_synthetic_residual.py --flip-sign; then
    echo "  НЕОЖИДАННО: скрипт завершился кодом 0 (ожидался 1)"
    exit 1
else
    echo "  ОЖИДАЕМЫЙ СБОЙ: код 1 (12/20 не совпали, знак обращён)"
fi

echo ""
echo "=== 6. CMD с изохронами и ΔG на label=2 ==="
$PYTHON experiment/turnoff/plot_turnoff_isochrones.py

echo ""
echo "=== Воспроизведение завершено. Сравнение моделей не производилось. ==="