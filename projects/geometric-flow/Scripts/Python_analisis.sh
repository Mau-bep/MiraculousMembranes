#!/bin/tcsh
#$ -cwd
#$ -j y



# set KB=0.005
foreach KB ( 0.05 0.01 0.005)
set Nsim=1

foreach v ( 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0 )
foreach Init_cond ( 1 )

python3 Read_E.py ${v} ${Init_cond} ${KB}

end
end



# foreach v ( 0.45 0.5 0.55 0.6 0.65 0.7 0.75 )
# # foreach v ( 0.475 0.525 0.575 0.625 0.675 0.725 )

# foreach Init_cond ( 2 )

# python3 Read_E.py ${v} ${Init_cond}

# end
# end


# foreach v ( 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 )
# foreach Init_cond ( 3 )

# python3 Read_E.py ${v} ${Init_cond}

# end
# end







end
