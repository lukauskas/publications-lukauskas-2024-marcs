#!/usr/bin/env bash
set -e

DIR_A=$1
DIR_B=$2

if [ -z "$DIR_A" ]
  then
    DIR_A="output/sorted/tables"
fi

if [ -z "$DIR_B" ]
  then
    DIR_B="expected-outputs/first-submission"
fi

echo "Comparing outputs from $DIR_A"
echo "To the outputs in $DIR_B"

WORKDIR=`mktemp -d -t compare-excels.XXXXXXXXXX`
if [[ ! "$WORKDIR" || ! -d "$WORKDIR" ]]; then
  echo "Could not create temp dir"
  exit 1
fi

clean_up () {
    if [ -d "$WORKDIR" ]; then
        echo "Removing $WORKDIR"
        rm -rf "$WORKDIR"
    fi
}
trap clean_up EXIT

echo "Working directory: $WORKDIR"

compare_excels () {
    local file_actual=$1;
    local file_expected=$2;

    if [[ ! -f "${file_actual}" ]]; then
        echo "ERROR: File from DIR_A does not exist: ${file_actual}"
        return 1
    fi


    if [ ! -f "${file_expected}" ]; then
        echo "ERROR: File from DIR_B does not exist: ${file_actual}"
        return 1
    fi

    cp "$file_actual" "${WORKDIR}/actual.xlsx"
    cp "$file_expected" "${WORKDIR}/expected.xlsx"

    if docker run --rm -e JAVA_OPTIONS='-Xmx8G' -v "${WORKDIR}":/wd:ro lukauskas/docker-excelcompare:0.6.1 --diff_numeric_precision=0.00000001 actual.xlsx expected.xlsx; then
        echo "SUCCESS: EXCEL FILES ${file_actual} and ${file_expected} match"
    else
        local EXITCODE=$?
        echo "FAILURE: EXCEL FILES ${file_actual} and ${file_expected} differ, or the script crashed. See above for differences/errors"
        return $EXITCODE
    fi
}

for ((i = 1; i <= 6; i++))
do
    echo "-- Table S${i} ---------------------------------------------------------------"
    echo "(you can ignore warnings about illegal access)"
    compare_excels "$DIR_A/Table S${i}.xlsx" "$DIR_B/Lukauskas-et-al-2021-Suppl-Table-${i}.xlsx"
    echo "---------------------------------------------------------------------"
    echo ""
done
