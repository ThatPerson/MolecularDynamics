open me.xyz
bond 1.55
hybridize 1 SP2
model GAFF 
thermostat ANDERSEN
energy
minimize 0.0001 0.00001 40000 1000 min.xyz
energy
heat 500 50 400 0.01 0.01 0.0001 0.00001 0.3 heat.xyz
temperature
energy
prod 0.0001 0.0001 500 .3 50000 10 prod.xyz
