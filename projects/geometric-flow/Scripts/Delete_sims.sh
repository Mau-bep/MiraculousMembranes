#!/bin/tcsh
#$ -cwd
#$ -j y

set start=123
set end=125

foreach Sim (`seq $start $end`)
    scancel $Sim
end





