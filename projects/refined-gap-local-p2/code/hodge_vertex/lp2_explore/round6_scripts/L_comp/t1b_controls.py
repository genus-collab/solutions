from t1_bilinear import *
Z=build(8,nshift=(0,1,2,-1))
for k in range(2,5):
    for d in (0,1,2):
        print(k,d,"(Z0,Z1)",nz(Pk_value(Z[0],Z[1],k,d)),"(Z0,Z2)",nz(Pk_value(Z[0],Z[2],k,d)),"(Z0,Z-1)",nz(Pk_value(Z[0],Z[-1],k,d)),"(Z0,Z0)",nz(Pk_value(Z[0],Z[0],k,d)))
