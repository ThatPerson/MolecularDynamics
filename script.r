open file.xyz
bond 1.55
%output output.xyz
model MM3
energy
minimize 0.01 0.01 20000 100 min.xyz
energy
heat 3
prod 0.01 0.01 3 .1 10000 100 prod.xyz