#!/usr/bin/env bash
set -e
set -x


DIR_A=$1
DIR_B=$2

echo $DIR_A $DIR_B

TABLES=("Table S1.xlsx" "Table S2.xlsx" "Table S3.xlsx" "Table S4.xlsx" "Table S5.xlsx")

for ((i = 0; i < ${#TABLES[@]}; i++))
do
    t="${TABLES[$i]}"
    echo "-- $t ---------------------------------------------------------------"
    excel_cmp "$DIR_A/tables/$t" "$DIR_B/tables/$t"
    echo "---------------------------------------------------------------------"
    echo ""
done
