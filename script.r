open big.xyz
bond 1.55
%output output.xyz
model MM3
energy
minimize 0.0001 0.00001 20000 100 min.xyz
energy
%heat .01
prod 0.0001 0.00001 100 0.1 50000 100 prod.xyz
