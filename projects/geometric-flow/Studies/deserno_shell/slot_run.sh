#!/usr/bin/env bash
# Global CPU semaphore shared by every agent/session of this study.
#   usage: slot_run.sh <command> [args...]
# Waits until one of NSLOTS lock files is free, runs the command with OMP_NUM_THREADS=1 at low priority
# while holding it. The machine has 8 cores and other people/sessions use it: never start a simulation
# without going through this script (and never raise NSLOTS above 6).
NSLOTS=${NSLOTS:-6}
LOCKDIR=${LOCKDIR:-/tmp/mm_deserno_slots}
mkdir -p "$LOCKDIR"
while true; do
  for i in $(seq 1 "$NSLOTS"); do
    exec 9>"$LOCKDIR/slot$i"
    if flock -n 9; then
      OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} nice -n 10 "$@"
      exit $?
    fi
    exec 9>&-
  done
  sleep 3
done
