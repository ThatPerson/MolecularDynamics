open propane.xyz
bond 1.55
%output output.xyz
model MM3
energy
minimize 0.01 0.01 20000 100 min.xyz
energy
heat .01
prod 0.01 0.001 0 0 50000 100 prod.xyz
