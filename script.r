%open big.xyz
%bond 1.55
%output output.xyz
%model MM3
%energy
%minimize 0.0001 0.00001 40000 100 min.xyz
%output big_min.xyz
%open new.xyz
%bond 1.55
%model MM3
%energy
%minimize 0.0001 0.00001 40000 1000 min.xyz
%output new_min.xyz
open big_min.xyz
bond 1.55
model MM3
thermostat ANDERSEN
energy
heat 300
temperature
energy
prod 0.0001 0.0001 300 .3 50000 10 prod_and.xyz
