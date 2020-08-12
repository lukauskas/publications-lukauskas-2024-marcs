#!/usr/bin/env bash

DIR_A=$1
DIR_B=$2

echo $DIR_A $DIR_B

excel_cmp $DIR_A/preprocessing/table-pulldowns.xlsx "$DIR_B/Table S1 - table-pulldowns.xlsx"
excel_cmp $DIR_A/preprocessing/table-heatmap.xlsx "$DIR_B/Table S2 - table-heatmap.xlsx"
excel_cmp $DIR_A/ptm-response/ptm-response.xlsx "$DIR_B/Table S3 - ptm-response.xlsx"
excel_cmp $DIR_A/networks/table-networks.xlsx "$DIR_B/Table S4 - table-networks.xlsx"
excel_cmp $DIR_A/ptm-response/ptm-response-complexes.xlsx "$DIR_B/Table S5 - ptm-response-complexes.xlsx"
