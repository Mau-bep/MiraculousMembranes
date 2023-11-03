import numpy as np
import matplotlib.pyplot as plt 



sigma=1.0
rc=2.0

r=np.linspace(0.5,2.0,100)
# rc=2
rc2=rc*rc

strength=0.1
alpha=2*(rc2/(sigma*sigma))*np.power(3/(2*( (rc2/(sigma*sigma)) -1)),3.0 ) 
# alpha=2*(rc2/(sigma*sigma))*pow( 3/(2*( (rc2/(sigma*sigma)) -1))  ,3.0 );
# alpha=2*(rc2/(sigma*sigma))*pow( 3/(2*( (rc2/(sigma*sigma)) -1))  ,3.0 );
F=strength*alpha*(-2)*( (sigma*sigma/(r*r*r))*np.pow(( (rc2/(r*r))-1),2.0 )+  + 2*((sigma*sigma/(r*r))-1)*( (rc2/(r*r))-1 )*(rc2/(r*r*r))     )

# Force=strength*alpha*(-2)*(  (sigma*sigma/(r*r*r))*pow(( (rc2/(r*r))-1),2.0 )+ 2*((sigma*sigma/(r*r))-1)*( (rc2/(r*r))-1 )*(rc2/(r*r*r)) )*dual_area*unit_r;
