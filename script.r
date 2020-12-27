open big.xyz
bond 1.55
model MM3
thermostat ANDERSEN
energy
minimize 0.0001 0.00001 40000 1000 min.xyz
energy
heat 300
temperature
energy
prod 0.0001 0.0001 300 .3 50000 10 prod.xyz
