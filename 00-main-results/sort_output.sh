#!/usr/bin/env bash
set -e

UNSORTED_OUTPUT_DIR=$1
SORTED_OUTPUT_DIR=$2
# To skip the confirmation
REPLY=$3

if [ -z $UNSORTED_OUTPUT_DIR ] | [ -z $SORTED_OUTPUT_DIR ]
then
    echo "Usage $0 [unsorted_output_dir] [sorted_output_dir] [optional: y to do not ask to override]"
    exit 1
fi
# If output directory exists
if [ -d ${SORTED_OUTPUT_DIR} ]
then
    # Add some warnings
    echo "warning: ${SORTED_OUTPUT_DIR} will be overwritten"

    # Check if the user specified REPLY
    if [ -z $REPLY ]
    then
        # Thanks https://stackoverflow.com/questions/1885525/how-do-i-prompt-a-user-for-confirmation-in-bash-script
        read -p "Are you sure? Type 'y' for yes" -n 1 -r
        echo
    else
        echo "Answer already specified: ${REPLY}. Not asking interactively"
    fi
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
        echo "Deleting directory";
        rm -rf "$SORTED_OUTPUT_DIR"
    else
        echo "No permission to delete unsorted directory. Quitting.";
        exit 1
    fi
fi

echo "Creating output directory"
mkdir -p "$SORTED_OUTPUT_DIR"

echo "Gathering data for Figure 1"
mkdir -p "$SORTED_OUTPUT_DIR/figure-01"

set -x
cp $UNSORTED_OUTPUT_DIR/scatterplots/scatter-H42-RdBu_r.pdf $SORTED_OUTPUT_DIR/figure-01/figure-01c.pdf
set +x

echo "Gathering data for Figure S1"
mkdir -p "$SORTED_OUTPUT_DIR/figure-S01"

set -x
cp $UNSORTED_OUTPUT_DIR/scatterplots/scatter-H03.pdf $SORTED_OUTPUT_DIR/figure-S01/figure-S01b.pdf
cp $UNSORTED_OUTPUT_DIR/scatterplots/scatter-H21.pdf $SORTED_OUTPUT_DIR/figure-S01/figure-S01c.pdf
cp $UNSORTED_OUTPUT_DIR/scatterplots/scatter-H40.pdf $SORTED_OUTPUT_DIR/figure-S01/figure-S01d.pdf
cp $UNSORTED_OUTPUT_DIR/scatterplots/grid-5-postprocessed.pdf $SORTED_OUTPUT_DIR/figure-S01/figure-S01e.pdf
set +x

echo "Gathering data for Figure 2"
mkdir -p "$SORTED_OUTPUT_DIR/figure-02"

set -x
cp $UNSORTED_OUTPUT_DIR/ptm-response/volcano-plots/H3K4me3-coloured.pdf $SORTED_OUTPUT_DIR/figure-02/figure-02b.pdf
cp $UNSORTED_OUTPUT_DIR/ptm-response/volcano-plots/H3K4me1-coloured.pdf $SORTED_OUTPUT_DIR/figure-02/figure-02c.pdf
cp $UNSORTED_OUTPUT_DIR/ptm-response/volcano-plots/H3K4me1-coloured.pdf $SORTED_OUTPUT_DIR/figure-02/figure-02c.pdf
cp $UNSORTED_OUTPUT_DIR/ptm-response/pairwise_scatterplot_H3K4me3vH3K4me1.pdf $SORTED_OUTPUT_DIR/figure-02/figure-02d.pdf
cp $UNSORTED_OUTPUT_DIR/ptm-response/significant-counts-barplots.pdf $SORTED_OUTPUT_DIR/figure-02/figure-02e.pdf
cp $UNSORTED_OUTPUT_DIR/ptm-response/significant-counts-histogram-truncated.pdf $SORTED_OUTPUT_DIR/figure-02/figure-02f.pdf
set +x

echo "Gathering data for Figure S2"
mkdir -p "$SORTED_OUTPUT_DIR/figure-S02"

set -x
cp $UNSORTED_OUTPUT_DIR/ptm-response/H3K4me3-ratio-travels.pdf $SORTED_OUTPUT_DIR/figure-S02/figure-S02a1.pdf
cp $UNSORTED_OUTPUT_DIR/ptm-response/H3K4me3-ratio-travels-huge.pdf $SORTED_OUTPUT_DIR/figure-S02/figure-S02a2.pdf
cp $UNSORTED_OUTPUT_DIR/ptm-response/network_df.pdf $SORTED_OUTPUT_DIR/figure-S02/figure-S02b.pdf
set +x


echo "Gathering data for Figure 3"
mkdir -p "$SORTED_OUTPUT_DIR/figure-03"

set -x
# No output for this generated in this pipeline
echo "No output from this pipeline in this figure" > $SORTED_OUTPUT_DIR/figure-03/intentionally_empty.txt
set +x

echo "Gathering data for Figure S3"
mkdir -p "$SORTED_OUTPUT_DIR/figure-S03"

set -x
cp "$UNSORTED_OUTPUT_DIR/ptm-response/volcano-plots/H3K27me2.pdf" "$SORTED_OUTPUT_DIR/figure-S03/figure-S03-r1c1-H3K27me2.pdf"
cp "$UNSORTED_OUTPUT_DIR/ptm-response/volcano-plots/H3K27me3.pdf" "$SORTED_OUTPUT_DIR/figure-S03/figure-S03-r1c2-H3K27me3.pdf"
cp "$UNSORTED_OUTPUT_DIR/ptm-response/volcano-plots/H3K9me2.pdf" "$SORTED_OUTPUT_DIR/figure-S03/figure-S03-r1c3-H3K9me2.pdf"
cp "$UNSORTED_OUTPUT_DIR/ptm-response/volcano-plots/H3K9me3.pdf" "$SORTED_OUTPUT_DIR/figure-S03/figure-S03-r1c4-H3K9me3.pdf"

cp "$UNSORTED_OUTPUT_DIR/ptm-response/volcano-plots/H3K9acK14ac.pdf" "$SORTED_OUTPUT_DIR/figure-S03/figure-S03-r2c1-H3K9acK14ac.pdf"
cp "$UNSORTED_OUTPUT_DIR/ptm-response/volcano-plots/H3K27ac.pdf" "$SORTED_OUTPUT_DIR/figure-S03/figure-S03-r2c2-H3K27ac.pdf"
cp "$UNSORTED_OUTPUT_DIR/ptm-response/volcano-plots/H3ac.pdf" "$SORTED_OUTPUT_DIR/figure-S03/figure-S03-r2c3-H3ac.pdf"
cp "$UNSORTED_OUTPUT_DIR/ptm-response/volcano-plots/DNA Methylation.pdf" "$SORTED_OUTPUT_DIR/figure-S03/figure-S03-r2c4-DNA Methylation.pdf"

cp "$UNSORTED_OUTPUT_DIR/ptm-response/volcano-plots/H4K16ac.pdf" "$SORTED_OUTPUT_DIR/figure-S03/figure-S03-r3c1-H4K16ac.pdf"
cp "$UNSORTED_OUTPUT_DIR/ptm-response/volcano-plots/H4K20me2.pdf" "$SORTED_OUTPUT_DIR/figure-S03/figure-S03-r3c2-H4K20me2.pdf"
cp "$UNSORTED_OUTPUT_DIR/ptm-response/volcano-plots/H4K20me3.pdf" "$SORTED_OUTPUT_DIR/figure-S03/figure-S03-r3c3-H4K20me3.pdf"
cp "$UNSORTED_OUTPUT_DIR/ptm-response/volcano-plots/H4ac.pdf" "$SORTED_OUTPUT_DIR/figure-S03/figure-S03-r3c4-H4ac.pdf"

set +x


echo "Gathering data for Figure 4"
mkdir -p "$SORTED_OUTPUT_DIR/figure-04"

set -x
cp $UNSORTED_OUTPUT_DIR/ptm-response/network-projections.pdf $SORTED_OUTPUT_DIR/figure-04/figure-04a.pdf
cp $UNSORTED_OUTPUT_DIR/ptm-response/pairwise_scatterplot_H3acvH4ac.pdf $SORTED_OUTPUT_DIR/figure-04/figure-04b.pdf
cp $UNSORTED_OUTPUT_DIR/ptm-response/complex-differences-plot.H3ac_vs_H4ac.pdf $SORTED_OUTPUT_DIR/figure-04/figure-04c.pdf
set +x

echo "Gathering data for Figure S4"
mkdir -p "$SORTED_OUTPUT_DIR/figure-S04"

set -x
cp "$UNSORTED_OUTPUT_DIR/ptm-response/volcano-plots/H2A.Z.pdf" "$SORTED_OUTPUT_DIR/figure-S04/figure-S04a.pdf"
cp "$UNSORTED_OUTPUT_DIR/ptm-response/pairwise_scatterplot_H4acvH2A.Z.pdf" "$SORTED_OUTPUT_DIR/figure-S04/figure-S04b.pdf"
cp "$UNSORTED_OUTPUT_DIR/ptm-response/pairwise_scatterplot_H3acvH2A.Z.pdf" "$SORTED_OUTPUT_DIR/figure-S04/figure-S04c.pdf"
set +x

echo "Gathering data for Figure S5"
mkdir -p "$SORTED_OUTPUT_DIR/figure-S05"

set -x
cp "$UNSORTED_OUTPUT_DIR/networks/graphs/graph.q.0.05.pdf" "$SORTED_OUTPUT_DIR/figure-S05/figure-S05ab-network-1-q0.05.pdf"
cp "$UNSORTED_OUTPUT_DIR/networks/graphs/graph.q.0.05.gephi.gexf" "$SORTED_OUTPUT_DIR/figure-S05/figure-S05ab-network-1-q0.05.gexf"

cp "$UNSORTED_OUTPUT_DIR/networks/graphs/graph.q.0.01.pdf" "$SORTED_OUTPUT_DIR/figure-S05/figure-S05ab-network-2-q0.01.pdf"
cp "$UNSORTED_OUTPUT_DIR/networks/graphs/graph.q.0.01.gephi.gexf" "$SORTED_OUTPUT_DIR/figure-S05/figure-S05ab-network-2-q0.01.gexf"

cp "$UNSORTED_OUTPUT_DIR/networks/graphs/graph.q.0.001.pdf" "$SORTED_OUTPUT_DIR/figure-S05/figure-S05ab-network-3-q0.001.pdf"
cp "$UNSORTED_OUTPUT_DIR/networks/graphs/graph.q.0.001.gephi.gexf" "$SORTED_OUTPUT_DIR/figure-S05/figure-S05ab-network-3-q0.001.gexf"

cp "$UNSORTED_OUTPUT_DIR/networks/graphs/graph.q.0.0001.pdf" "$SORTED_OUTPUT_DIR/figure-S05/figure-S05ab-network-4-q0.0001.pdf"
cp "$UNSORTED_OUTPUT_DIR/networks/graphs/graph.q.0.0001.gephi.gexf" "$SORTED_OUTPUT_DIR/figure-S05/figure-S05ab-network-4-q0.0001.gexf"

cp "$UNSORTED_OUTPUT_DIR/networks/graphs/graph.q.high-confidence.pdf" "$SORTED_OUTPUT_DIR/figure-S05/figure-S05ab-network-5-hc.pdf"
cp "$UNSORTED_OUTPUT_DIR/networks/graphs/graph.q.high-confidence.gephi.gexf" "$SORTED_OUTPUT_DIR/figure-S05/figure-S05ab-network-5-hc.gexf"


cp "$UNSORTED_OUTPUT_DIR/networks/training-score-summary.csv" "$SORTED_OUTPUT_DIR/figure-S05/figure-S05ab-network-stats.csv"

cp "$UNSORTED_OUTPUT_DIR/networks/training-prc-truncated.pdf" "$SORTED_OUTPUT_DIR/figure-S05/figure-S05b.pdf"
cp "$UNSORTED_OUTPUT_DIR/networks/training-score-summary.csv" "$SORTED_OUTPUT_DIR/figure-S05/figure-S05b.csv"

cp "$UNSORTED_OUTPUT_DIR/networks/training-score_breakdown_by_publication_count.pdf" "$SORTED_OUTPUT_DIR/figure-S05/figure-S05c.pdf"
cp "$UNSORTED_OUTPUT_DIR/networks/training-score_breakdown_by_experimental_system.pdf" "$SORTED_OUTPUT_DIR/figure-S05/figure-S05d.pdf"

cp "$UNSORTED_OUTPUT_DIR/networks/graphs/graph.q.0.001.pdf" "$SORTED_OUTPUT_DIR/figure-S05/figure-S05e.pdf"
cp "$UNSORTED_OUTPUT_DIR/networks/graphs/graph.q.0.0001.gephi.gexf" "$SORTED_OUTPUT_DIR/figure-S05/figure-S05e.gexf"

set +x

echo "Gathering data for Figure 5"
mkdir -p "$SORTED_OUTPUT_DIR/figure-05"

set -x
cp $UNSORTED_OUTPUT_DIR/networks/graphs/graph.q.high-confidence.pdf $SORTED_OUTPUT_DIR/figure-05/figure-05a.pdf
cp $UNSORTED_OUTPUT_DIR/ptm-response/barplots/barplot-horizontal-ino80_exclusive_subunits_.pdf $SORTED_OUTPUT_DIR/figure-05/figure-05d.pdf
set +x


echo "Gathering data for Figure S6"
mkdir -p "$SORTED_OUTPUT_DIR/figure-S06"

set -x
cp "$UNSORTED_OUTPUT_DIR/networks/graphs/graph.q.high-confidence.pdf" "$SORTED_OUTPUT_DIR/figure-S06/figure-S06.pdf"
cp "$UNSORTED_OUTPUT_DIR/networks/graphs/graph.q.high-confidence.gephi.gexf" "$SORTED_OUTPUT_DIR/figure-S06/figure-S06.gexf"
set +x

echo "Gathering data for Figure S7"
mkdir -p "$SORTED_OUTPUT_DIR/figure-S07"

set -x
# No output for this generated in this pipeline
echo "No output from this pipeline in this figure" > $SORTED_OUTPUT_DIR/figure-S07/intentionally_empty.txt
set +x


echo "Gathering tables"
mkdir -p "$SORTED_OUTPUT_DIR/tables"
mkdir -p "$SORTED_OUTPUT_DIR/tables/as_tsv"

set -x
cp $UNSORTED_OUTPUT_DIR/preprocessing/table-pulldowns.xlsx "$SORTED_OUTPUT_DIR/tables/Table S1.xlsx"
cp $UNSORTED_OUTPUT_DIR/preprocessing/table-pulldowns.sheet.*.tsv.gz "$SORTED_OUTPUT_DIR/tables/as_tsv/"
rename 's/table-pulldowns\./table-s1./' $SORTED_OUTPUT_DIR/tables/as_tsv/*.tsv.gz

cp $UNSORTED_OUTPUT_DIR/preprocessing/table-heatmap.xlsx "$SORTED_OUTPUT_DIR/tables/Table S2.xlsx"
cp $UNSORTED_OUTPUT_DIR/preprocessing/table-heatmap.sheet.*.tsv.gz "$SORTED_OUTPUT_DIR/tables/as_tsv/"
rename 's/table-heatmap\./table-s2./' $SORTED_OUTPUT_DIR/tables/as_tsv/*.tsv.gz

cp $UNSORTED_OUTPUT_DIR/ptm-response/ptm-response.xlsx "$SORTED_OUTPUT_DIR/tables/Table S3.xlsx"
cp $UNSORTED_OUTPUT_DIR/ptm-response/ptm-response.sheet.*.tsv.gz "$SORTED_OUTPUT_DIR/tables/as_tsv/"
rename 's/ptm-response\./table-s3./' $SORTED_OUTPUT_DIR/tables/as_tsv/*.tsv.gz

cp $UNSORTED_OUTPUT_DIR/ptm-response/ptm-response-heatmap.xlsx "$SORTED_OUTPUT_DIR/tables/Table S4.xlsx"
cp $UNSORTED_OUTPUT_DIR/ptm-response/ptm-response-heatmap.sheet.*.tsv.gz "$SORTED_OUTPUT_DIR/tables/as_tsv/"
rename 's/ptm-response-heatmap\./table-s4./' $SORTED_OUTPUT_DIR/tables/as_tsv/*.tsv.gz

cp $UNSORTED_OUTPUT_DIR/networks/table-networks.xlsx "$SORTED_OUTPUT_DIR/tables/Table S5.xlsx"
cp $UNSORTED_OUTPUT_DIR/networks/table-networks.sheet.*.tsv.gz "$SORTED_OUTPUT_DIR/tables/as_tsv/"
rename 's/table-networks\./table-s5./' $SORTED_OUTPUT_DIR/tables/as_tsv/*.tsv.gz

cp $UNSORTED_OUTPUT_DIR/ptm-response/ptm-response-complexes.xlsx "$SORTED_OUTPUT_DIR/tables/Table S6.xlsx"
cp $UNSORTED_OUTPUT_DIR/ptm-response/ptm-response-complexes.sheet.*.tsv.gz "$SORTED_OUTPUT_DIR/tables/as_tsv/"
rename 's/ptm-response-complexes\./table-s6./' $SORTED_OUTPUT_DIR/tables/as_tsv/*.tsv.gz

set +x
