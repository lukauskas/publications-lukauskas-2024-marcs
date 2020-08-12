#!/usr/bin/env bash

DIR_A=$1
DIR_B=$2

echo $DIR_A $DIR_B

excel_cmp $DIR_A/networks/table-networks.xlsx $DIR_B/networks/table-networks.xlsx
excel_cmp $DIR_A/curated_complexes.xlsx $DIR_B/curated_complexes.xlsx
excel_cmp $DIR_A/preprocessing/table-pulldowns.xlsx $DIR_B/preprocessing/table-pulldowns.xlsx
excel_cmp $DIR_A/preprocessing/table-heatmap.xlsx $DIR_B/preprocessing/table-heatmap.xlsx
excel_cmp $DIR_A/ptm-response/ptm-response-complexes.xlsx $DIR_B/ptm-response/ptm-response-complexes.xlsx
excel_cmp $DIR_A/ptm-response/ptm-response.xlsx $DIR_B/ptm-response/ptm-response.xlsx
