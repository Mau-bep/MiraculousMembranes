#!/bin/bash
# Collects what is worth downloading from the chains of one batch (one-bead chains and tube chains alike).
#
#   bash Collect_chain_results.sh ../Results/Phase_one_flat                 # every chain folder of the batch (KA_0, KA_0.06667, ...)
#   bash Collect_chain_results.sh ../Results/Phase_one_flat KA_0.2 KA_1     # only these
#   DO_OBJ=0 bash Collect_chain_results.sh ../Results/Phase_one_vesicle     # only the csv's
#   for b in ../Results/Phase_one_* ; do bash Collect_chain_results.sh $b ; done
#
# For every chain folder <BATCH>/<TAG>/ (the one that has a chain_summary.csv) it makes
#   <BATCH>/Chain_data_<TAG>.csv                          a copy of the chain_summary.csv, so that the csv's of all the chains sit in one folder
#   <BATCH>/<TAG>/Final_states/<TAG>_Final_link_NNN.obj   the Final_state.obj of every POINT of the chain, in the order of the chain
# The files are COPIES: the chain folder stays as it is (a chain resumes from its Final_state.obj files). Run it again whenever more of the chain is
# finished; the files are overwritten.
#
# The points are the links of the list, not the transitional runs:
#   tube chains      chain_summary.csv has a column "kind": the rows with kind = point                      (waypoint and mid are left out)
#   one-bead chains  no such column: the rows whose run folder is not <...>_w1, _w2, ... (waypoints) nor <...>_mid (half step of a retry)
# and the folder of a run is the column "rundir" (tube) or, in the one-bead chains, the folder named in "run" with the run folder (1, 2, ...) inside.
# Change the awk lines in points() below to select other rows (e.g. only the way in: add   if ($col["leg"] != "in") next ).

FINAL_DIR=${FINAL_DIR:-Final_states}   # folder of the meshes, inside the chain folder
PAD=${PAD:-3}                          # digits of the link number in the file name (007), so that the files sort in the order of the chain
DO_CSV=${DO_CSV:-1}                    # 0: do not copy the csv's
DO_OBJ=${DO_OBJ:-1}                    # 0: do not copy the meshes

BATCH=${1:?usage: bash Collect_chain_results.sh BATCH_DIR [TAG ...]   (BATCH_DIR holds the chain folders, e.g. ../Results/Phase_one_flat)}
BATCH=${BATCH%/}
shift
if [ $# -gt 0 ]; then
    TAGS=("$@")
else
    TAGS=()
    for s in "$BATCH"/*/chain_summary.csv; do
        [ -f "$s" ] && TAGS+=("$(basename "$(dirname "$s")")")
    done
fi
if [ ${#TAGS[@]} -eq 0 ]; then
    echo "no <tag>/chain_summary.csv in $BATCH (BATCH_DIR is the folder that holds the chain folders)" >&2
    exit 1
fi

# prints "link rundir" for every point of the chain_summary.csv given
points() {
    awk -F, '
        { sub(/\r$/, "") }
        NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
        {
            run = $col["run"]
            if ("kind" in col) { if ($col["kind"] != "point") next }
            else if (run ~ /_w[0-9]+$/ || run ~ /_mid$/) next
            print $col["link"], (("rundir" in col) ? $col["rundir"] : run)
        }' "$1"
}

for TAG in "${TAGS[@]}"; do
    CHAIN="$BATCH/$TAG"
    if [ ! -f "$CHAIN/chain_summary.csv" ]; then
        echo "$TAG: no chain_summary.csv in $CHAIN, skipped"
        continue
    fi
    msg="$TAG:"
    if [ "$DO_CSV" = 1 ]; then
        cp "$CHAIN/chain_summary.csv" "$BATCH/Chain_data_$TAG.csv" && msg="$msg Chain_data_$TAG.csv"
    fi
    if [ "$DO_OBJ" = 1 ]; then
        mkdir -p "$CHAIN/$FINAL_DIR"
        copied=0
        missing=0
        while read -r LINK RUNDIR; do
            SRC="$CHAIN/$RUNDIR"
            # one-bead chains: the folder holds the numbered run folders, the last one is the run
            if [ ! -f "$SRC/Final_state.obj" ]; then
                SRC="$SRC/$(ls "$SRC" 2>/dev/null | grep -E '^[0-9]+$' | sort -n | tail -1)"
            fi
            if [ -f "$SRC/Final_state.obj" ]; then
                printf -v NAME "%s_Final_link_%0*d.obj" "$TAG" "$PAD" "$((10#$LINK))"
                cp "$SRC/Final_state.obj" "$CHAIN/$FINAL_DIR/$NAME"
                copied=$((copied + 1))
            else
                echo "$TAG: link $LINK: no Final_state.obj in $CHAIN/$RUNDIR" >&2
                missing=$((missing + 1))
            fi
        done < <(points "$CHAIN/chain_summary.csv")
        msg="$msg $copied final states -> $TAG/$FINAL_DIR/"
        [ $missing -gt 0 ] && msg="$msg ($missing missing)"
    fi
    echo "$msg"
done
