* Test 4-way coupling
Vin in 0 10
L1 1 0 1m
L2 2 0 1m
L3 3 0 1m
L4 4 0 1m
K1 L1 L2 L3 L4 0.9
R1 1 0 10
R2 2 0 10
R3 3 0 10
R4 4 0 10
.tran 1u 10u
.end
