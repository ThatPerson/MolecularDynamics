open water.xyz
bond 1.00
model TIP3P
thermostat ANDERSEN
energy
minimize 0.0001 0.00001 40000 1000 min.xyz
energy
heat 300
temperature
energy
prod 0.0001 0.0001 300 1 50000 10 prod.xyz
