#!/usr/bin/env bash
#
# reproduce.sh — воспроизвести CMD области поворота и отбор субгигантов
#                из исходных данных (data/raw/m67_members_hr2023.tsv).
#
# Все команды из корня проекта.

PYTHON="/usr/local/bin/python3.9"
set -euo pipefail
cd "$(git rev-parse --show-toplevel 2>/dev/null || echo "$(dirname "$0")/../..")"

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
echo "=== 5. Проверка 2: синтетический остаток (обращённый знак) ==="
$PYTHON experiment/turnoff/check_synthetic_residual.py --flip-sign

echo ""
echo "=== Воспроизведение завершено. Сравнение моделей не производилось. ==="