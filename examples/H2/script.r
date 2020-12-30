open H2.xyz
bond 1.55
model GAFF
thermostat ANDERSEN
energy
minimize 0.0001 0.00001 40000 1000 min.xyz
energy
heat 300 50 400 0.01 0.01 0.0001 0.00001 0.3 heat.xyz
temperature
energy
prod 0.0001 0.0001 300 .3 50000 10 prod.xyz
