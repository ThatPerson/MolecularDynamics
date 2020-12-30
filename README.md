MolecularDynamics
=================

![](boat.gif)

Very simple molecular dynamics model. Currently attempts to implement a simplified version of the GAFF force field (e.g. no different atom environments). Has Andersen and Langevin thermostat kind of implemented, but Langevin doesn't really work.

  

Todo
----

* Fix units and get thermal control working
* Implement TIP3P water model
* Add the ability to have different models (e.g. define water and the molecule separately so bonding can be done)
* Different bond lengths for different atoms

Compilation
-----------

`gcc main.c -lm -fopenmp -o md -O3`

May be compiled with or without OpenMP support (this parallelises the atom loop).

Operation
---------

The user interface is as in `ThatPerson/Conformation`. In scripting mode (eg `./md script.r`) the following commands are implemented

* `open <filename>` Opens XYZ file `<filename>`
* `bond <bl A>` Introduced bonds between all atoms less than `<bl A>` apart (in Angstrom)
* `output <filename>` Writes XYZ file output into `<filename>`
* `model <model>` Sets up model `<model>`. Currently only `<model> = MM3` is allowed.
* `energy` Prints out the energy of the system calculated using the forcefield. 
* `minimize <dn> <dt> <steps> <outputsteps> <filename>` Runs minimization (e.g. velocity deleted each cycle, so system bond lengths and angles are relaxed) for `<steps>`, writing out the system to `<filename>` every `<outputsteps>`. Differentiation of the forcefield is done using `<dn>`, where this is a distance in Angstroms (generally ~0.01). Each step is `<dt>` picoseconds (generally 0.01).
* `heat <temperature> <temp_step> <steps> <dn> <dt> <dnM> <dtN> <viscosity> <filename>` Increments the temperature by `<temp_step>`, jumps forward `<dn>` `<dt>`, then runs minimization for `<steps>` with `dn=<dnM>` and `dt=<dtN>`. Repeats this until the tempeature reaches `<temperature>`
* `iterate <dn> <dt> <temp> <viscosity> <filename>` Jumps forward `<dn>` `<dt>`, allowing for more coarse iteration.
* `hybridize <id> <sp>` Changes atom hybridization state. This changes the expected bond angles. Atom ID is counted from XYZ file, indexed at 0. `<sp>` may be SP, SP2, or SP3 depending on site.
* `prod <dn> <dt> <temperature> <viscosity> <step> <outputstep> <filename>` Runs production MD. All inputs are as in `minimize`, except for `<temperature>` (*see warning*) and `<viscosity>` (*see warning*).
* `temperature` Outputs current system temperature
* `thermostat <thermostat>` `<thermostat>` can either be ANDERSEN or LANGEVIN. LANGEVIN is currently broken, see warning. With ANDERSEN temperature control seems okay but can explode. In ANDERSEN mode, the `v` parameter is taken from `viscosity`.



