#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#define TRUE 1
#define FALSE 0

#define ATOM_H 0
#define ATOM_C 1
#define ATOM_O 2
#define ATOM_N 3
#define ATOM_UNKNOWN 4
#define ATOM_TYPES 5

#define MOD_MM3 1
#define MOD_TIP3P 2
#define MOD_GAFF 3

#define INITIAL_BONDS 5 // Number of bonds per atom possible.
						// If you've got some weird compound may need to increase

	/* 1 kcal/mol = 4184 J/mol
	 * = 4184 kg m^2 s^-2 mol^-1
	 * = 4184 / (1.66*10^-27) amu m^2 s^-2 mol^-1
	 * = (4184 * (10^10)^2) / (1.66*10^-27) amu A^2 s^-2 mol^-1
	 * = (4184 * (10^10)^2 * (10^12)^-2) / (1.66*10^-27) amu A^2 ps^-2 mol^-1
	 * = (4184 * (10^10)^2 * (10^12)^-2) / (1.66*10^-27 * 6.02*10^23) amu A^2 ps^-2 
	 * = 418.68 amu A^2 ps^-2
	 */
#define KCAL_ATOMIC 418.68
#define BOLTZMANN_K 0.001985875 * KCAL_ATOMIC // kb = 0.001985875 kcal mol^-1 K^-1, * 418.68 conversion.

// CC = 9.987551 * 10^9 N m^2 C^-2
// 1 J = 1 N m
// eg units are J m C^-2
// so CC / 4184 kcal m / C^2 mol
// CC * 10^-10 / 4184 kcal Ang / C^2 mol
// = 
#define COULOMB_CONSTANT 2.3871e-04

#define LANGEVIN 0
#define ANDERSEN 1

#define SP1  0
#define SP2 1
#define SP3 2


struct Vector {
	double x, y, z;
};

struct Model {
	double eq_bond_length[ATOM_TYPES][ATOM_TYPES];
	double eq_bond_K[ATOM_TYPES][ATOM_TYPES];
	double angle_C[ATOM_TYPES];
	double angle_Z[ATOM_TYPES];
	double charge[ATOM_TYPES];
	double vdwM[ATOM_TYPES];
	double polarizability[ATOM_TYPES];
	double dihedral_pref[ATOM_TYPES][ATOM_TYPES];
	double dihedral_n[ATOM_TYPES][ATOM_TYPES];
	double dihedral_gamma[ATOM_TYPES][ATOM_TYPES];
	double vdwR[ATOM_TYPES];
	double vdw_EDEP[ATOM_TYPES];
};

struct Atom {
	struct Vector v;
	struct Vector vel;
	unsigned int n_bonds;
	int check;
	unsigned int lim_bonds;
	struct Atom ** bonds;
	char name[20];
	int type;
	int hybridization;
	unsigned int i;
	int andersen_f;
};

/* Molecule container - n_atoms in as struct. */
struct Molecule {
	struct Atom *as;
	unsigned int n_atoms;
	struct Model *model;
	int thermostat;
};

/** 
 * Generates uniform random long double from 0 to 1 inclusive.
 * @return long double
 *  Long double containing uniform random number.
 */
long double uniform_rand(void) {
	return ((long double) rand() + 1.) / ((long double) RAND_MAX + 1.);
}

/**
 * Uses Box-Muller method to generate a normally distributed random number.
 * @param mean
 *  Mean of normal distribution to select random variable from
 * @param std
 *  Standard deviation of normal distribution from which random variable selected
 * @return float
 *  Returns normally distributed random number
 */
double norm_rand(double mean, double std) {
	// Box-Muller method
	double rnd1, rnd2;
	rnd1 = (double) uniform_rand();
	rnd2 = (double) uniform_rand();
	double unadj = sqrt(-2 * log(rnd1)) * cos(2 * M_PI * rnd2);

	//printf("%f, %f, %f\n", rnd1, rnd2, mean + std*unadj);
	return mean + std * unadj;
}

void setup_model(struct Model * m, int model) {
	int i, k;

	for (i = 0; i < ATOM_TYPES; i++) {
		for (k = 0; k < ATOM_TYPES; k++) {
			m->eq_bond_length[i][k] = 0;
			m->eq_bond_K[i][k] = 0;
			m->dihedral_pref[i][k] = 0;
			m->dihedral_n[i][k] = 0;
			m->dihedral_gamma[i][k] = 0;
		}
		m->charge[i] = 0;
		m->angle_C[i] = 0;
		m->angle_Z[i] = 0;
		m->vdwM[i] = 0;
		m->polarizability[i] = 0;
		m->vdwR[i] = 0;
		m->vdw_EDEP[i] = 0;
	}

	if (model == MOD_GAFF) {
		m->eq_bond_length[ATOM_H][ATOM_H] = 0.738;
		m->eq_bond_K[ATOM_H][ATOM_H] = 4.661;
		
		m->eq_bond_length[ATOM_H][ATOM_C] = 1.090;
		m->eq_bond_K[ATOM_H][ATOM_C] = 6.217;
		m->eq_bond_length[ATOM_C][ATOM_H] = 1.090;
		m->eq_bond_K[ATOM_C][ATOM_H] = 6.217;
		
		m->eq_bond_length[ATOM_O][ATOM_H] = 0.960;
		m->eq_bond_K[ATOM_O][ATOM_H] = 5.794;
		m->eq_bond_length[ATOM_H][ATOM_O] = 0.960;
		m->eq_bond_K[ATOM_H][ATOM_O] = 5.794;
		
		m->eq_bond_length[ATOM_N][ATOM_H] = 1.010;
		m->eq_bond_K[ATOM_N][ATOM_H] = 6.057;
		m->eq_bond_length[ATOM_H][ATOM_N] = 1.010;
		m->eq_bond_K[ATOM_H][ATOM_N] = 6.057;
		
		m->eq_bond_length[ATOM_C][ATOM_O] = 1.440;
		m->eq_bond_K[ATOM_C][ATOM_O] = 7.347;
		m->eq_bond_length[ATOM_O][ATOM_C] = 1.440;
		m->eq_bond_K[ATOM_O][ATOM_C] = 7.347;
		
		m->eq_bond_length[ATOM_C][ATOM_N] = 1.470;
		m->eq_bond_K[ATOM_C][ATOM_N] = 7.504;
		m->eq_bond_length[ATOM_N][ATOM_C] = 1.470;
		m->eq_bond_K[ATOM_N][ATOM_C] = 7.504;
		
		m->eq_bond_length[ATOM_C][ATOM_C] = 1.526;
		m->eq_bond_K[ATOM_C][ATOM_C] = 7.643;
		
		m->angle_C[ATOM_C] = 1.339;
		m->angle_C[ATOM_N] = 1.300;
		m->angle_C[ATOM_O] = 1.249;
		
		m->angle_Z[ATOM_C] = 1.183;
		m->angle_Z[ATOM_N] = 1.212;
		m->angle_Z[ATOM_O] = 1.219;
		m->angle_Z[ATOM_H] = 0.784;
		
		m->vdwM[ATOM_C] = 12.;
		m->vdwM[ATOM_H] = 1.008;
		m->vdwM[ATOM_O] = 16.;
		m->vdwM[ATOM_N] = 14.;
		
		m->polarizability[ATOM_C] = 0.360;
		m->polarizability[ATOM_N] = 0.530;
		m->polarizability[ATOM_H] = 0.135;
		m->polarizability[ATOM_O] = 0.465;
		
		m->dihedral_pref[ATOM_C][ATOM_C] = 1.2/4.;
		m->dihedral_pref[ATOM_C][ATOM_N] = 10./4.;
		m->dihedral_pref[ATOM_N][ATOM_C] = 10./4.;
		m->dihedral_pref[ATOM_C][ATOM_O] = 4.6/4.;
		m->dihedral_pref[ATOM_O][ATOM_C] = 4.6/4.;
		
		m->dihedral_gamma[ATOM_C][ATOM_C] = 180. * (M_PI / 180.);
		m->dihedral_gamma[ATOM_C][ATOM_N] = 180. * (M_PI / 180.);
		m->dihedral_gamma[ATOM_N][ATOM_C] = 180. * (M_PI / 180.);
		m->dihedral_gamma[ATOM_C][ATOM_O] = 180. * (M_PI / 180.);
		m->dihedral_gamma[ATOM_O][ATOM_C] = 180. * (M_PI / 180.);
		
		m->dihedral_n[ATOM_C][ATOM_C] = 2.;
		m->dihedral_n[ATOM_C][ATOM_N] = 2.;
		m->dihedral_n[ATOM_N][ATOM_C] = 2.;
		m->dihedral_n[ATOM_C][ATOM_O] = 2.;
		m->dihedral_n[ATOM_O][ATOM_C] = 2.;
		
		m->vdwR[ATOM_H] = 1.3870;
		m->vdwR[ATOM_C] = 1.9080;
		m->vdwR[ATOM_O] = 1.6612;
		m->vdwR[ATOM_N] = 1.8240;
		
		m->vdw_EDEP[ATOM_H] = 0.0157;
		m->vdw_EDEP[ATOM_C] = 0.0860;
		m->vdw_EDEP[ATOM_N] = 0.1700;
		m->vdw_EDEP[ATOM_O] = 0.2100;
	}
	return;
}


/* reset_check()
 *  Resets all .check values to be FALSE. Should be run before _any_
 *  graph or rotate commands
 
 * Error modes;
 *  - If not all atoms are allocated, will segfault.
 */
void reset_check(struct Molecule *m) {
	unsigned int i;
	for (i = 0; i < m->n_atoms; i++) {
		m->as[i].check = FALSE;
	}
	return;
}


/* add_atom()
 *  Initialises atom *a.
 * arg *a: Pointer to atom struct to put data in.
 * arg x, y, z: positional coordinated.
 * arg name: atom name
 * returns: 1 in all cases.
 
 * Error modes;
 *  - If atom at *a not allocated memory, will segfault.
 */
int add_atom(struct Atom *a, float x, float y, float z, char name[20]) {
	a->v.x = x;
	a->v.y = y;
	a->v.z = z;
	a->vel.x = 0;
	a->vel.y = 0;
	a->vel.z = 0;
	a->n_bonds = 0;
	a->lim_bonds = INITIAL_BONDS;
	int c_buf = 0;
	while (name[c_buf] == ' ')
		c_buf++;
	strcpy(a->name, name+c_buf);
	a->check = FALSE;
	if (a->name[0] == 'H') {
		a->type = ATOM_H;
	} else if (a->name[0] == 'C'){
		a->type = ATOM_C;
	} else if (a->name[0] == 'O') {
		a->type = ATOM_O;
	} else {
		printf("Atom %s not added.\n", name);
		a->type = ATOM_UNKNOWN;
	}
	
	a->bonds = (struct Atom **) malloc(sizeof(struct Atom *) * INITIAL_BONDS);
	return 1;
}

/* add_bond()
 *  Creates bond between atom *a and *b. 
 * arg *a: Pointer to atom.
 * arg *b: Pointer to atom.
 * returns: 1 in all cases
 
 * Error modes;
 *  - Will seg fault if *a or *b do not exist, or if there is insufficient memory.
 */
int add_bond(struct Atom *a, struct Atom *b) {
	a->bonds[a->n_bonds] = b;
	b->bonds[b->n_bonds] = a;
	a->n_bonds++;
	b->n_bonds++;
	if (a->n_bonds >= a->lim_bonds) {
		a->lim_bonds *= 2;
		a->bonds = (struct Atom **) realloc(a->bonds, sizeof(struct Atom *) * a->lim_bonds);
	}
	if (b->n_bonds >= b->lim_bonds) {
		b->lim_bonds *= 2;
		b->bonds = (struct Atom **) realloc(b->bonds, sizeof(struct Atom *) * b->lim_bonds);
	}
	return 1;
}

void set_vector(struct Vector *a, double x, double y, double z) {
	a->x = x;
	a->y = y;
	a->z = z;
	return;
}

/* free_atoms()
 *  Frees bond arrays for each atom within molecule *m
 *
 * Error modes;
 *  - Will double free if any atom has not been initialized.
 */
void free_atoms(struct Molecule *m) {
	unsigned int i;
	for (i = 0; i < m->n_atoms; i++) {
		free(m->as[i].bonds);
	}
	return;
}

void crash(struct Molecule *m) {
	free_atoms(m);
	free(m->as);
	exit(-1);
}

/* print_moleculef()
 *  Prints molecular graph starting at atom n.
 *  Iterates in a functional manner, setting check to TRUE once an
 *  atom has been visited. Visits each atom once.
 * arg *a: Pointer to initial atom.
 * arg  n: Atom to begin graph at
 */
void print_moleculef(struct Atom * a, unsigned int n) {
	// if we have been visited, terminate execution.
	if (a->check == TRUE)
		return;
	a->check = TRUE;
	unsigned int i,j;
	for (j = 0; j < n+1; j++) 
		printf("> ");
	printf("Atom %s [%f, %f, %f] %d\n", a->name, a->v.x, a->v.y, a->v.z, a->n_bonds);
	
	// Recursively loop over, calling self for each neighbour.
	for (i = 0; i < a->n_bonds; i++) {
		print_moleculef(a->bonds[i], n+1);
	}
	return;
}

/* print_molecule()
 *  Prints molecular system, starting with one atom and looping to others
 */
void print_molecule(struct Atom * a) {
	if (a->check == TRUE)
		return;
	a->check = TRUE;
	unsigned int i;
	printf("%f, %f, %f\n", a->v.x, a->v.y, a->v.z);
	
	for (i = 0; i < a->n_bonds; i++) {
		print_molecule(a->bonds[i]);
	}
	return;
}


/* sub_vector()
 *  Subtracts vector *b from *a, and outputs into *res
 * arg *a: Vector
 * arg *b: Vector
 * arg *res: Vector, res = a - b
 */
void sub_vector(struct Vector *a, struct Vector *b, struct Vector *res) {
	// res = a - b
	res->x = a->x - b->x;
	res->y = a->y - b->y;
	res->z = a->z - b->z;
	return;
}

/* add_vector()
 *  Adds vector *b to *a, and outputs into *res
 * arg *a: Vector
 * arg *b: Vector
 * arg *res: Vector, res = a + b
 */
void add_vector(struct Vector *a, struct Vector *b, struct Vector *res) {
	// res = a + b
	res->x = a->x + b->x;
	res->y = a->y + b->y;
	res->z = a->z + b->z;
	return;
}


/* magnitude()
 *  Calculated magnitude of vector *v.
 * arg *v: Vector
 * returns: Magnitude
 */
float magnitude(struct Vector *v) {
	float n;
	n = powf(v->x, 2.) + powf(v->y, 2.) + powf(v->z, 2.);
	n = sqrtf(n);
	return n;
}

/* normalise()
 *  Normalise vector *v in situ
 * arg *v: Vector to be normalised.
 * 
 * Error modes
 *  - If vector has 0 magnitude, will do nothing.
 */
void normalise(struct Vector *v) {
	float n = magnitude(v);
	if (n == 0) {
		printf("Division by 0\n");
		return;
	}
	v->x /= n;
	v->y /= n;
	v->z /= n;
	return;
}

double dot(struct Vector *a, struct Vector *b) {
	return (a->x * b->x) + (a->y * b->y) + (a->z * b->z);
}

void print_vector(struct Vector *v) {
	printf("(%f, %f, %f)\n", v->x, v->y, v->z);
}

double calc_phi(struct Vector *a, struct Vector *b, struct Vector *c) {
	// calculates the ABC angle.
	struct Vector ab; 
	sub_vector(a, b, &ab);
	struct Vector bc;
	sub_vector(c, b, &bc);
	normalise(&ab);
	normalise(&bc);

	double v = dot(&ab, &bc);

	
	// because v is a double it can be -1.0000001
	if (v > 1)
		v = 1;
	if (v < -1)
		v = -1;
	return acos(v);
}

double calc_omega(struct Vector *a, struct Vector *b, struct Vector *c, struct Vector *d) {
	// calculates the ABC angle.
	struct Vector ab, cd;
	sub_vector(a, b, &ab);
	sub_vector(d, c, &cd);

	normalise(&ab);
	normalise(&cd);
	double v = dot(&ab, &cd);
	if (v > 1)
		v = 1;
	if (v < -1)
		v = -1;
	return acos(v);
}

double calc_distance(struct Vector *a, struct Vector *b) {
	double su = 0;
	su += powl(a->x - b->x, 2.);
	su += powl(a->y - b->y, 2.);
	su += powl(a->z - b->z, 2.);
	return sqrtl(su);
}


double bond_energy(struct Atom *a, struct Atom *b, struct Vector *apos, struct Vector *bpos, struct Model *m) {
	double req = m->eq_bond_length[a->type][b->type];
	double r = calc_distance(apos, bpos);
	double K = m->eq_bond_K[a->type][b->type];
	return K * powl(r - req, 2.);
}

double angle_energy(struct Atom *a, struct Atom *b, struct Atom *c, struct Vector *apos, struct Vector *bpos, struct Vector *cpos, struct Model *m) {
	// a-b-c
	double phi0;
	switch (b->hybridization) {
		case SP1: phi0 = 180. * (M_PI / 180.);break;
		case SP2: phi0 = 120. * (M_PI / 180.);break;
		case SP3: phi0 = 104. * (M_PI / 180.); break;
		default: phi0 = 104. * (M_PI / 180.); break;
	}
	double phi = calc_phi(apos, bpos, cpos);
	double rijeq = m->eq_bond_length[a->type][b->type];
	double rjkeq = m->eq_bond_length[b->type][c->type];
	
	double D = powl(rijeq - rjkeq, 2.) * powl(rijeq + rjkeq, -2.);
	double Kijk = 143.9 * m->angle_Z[a->type] * m->angle_C[b->type] * m->angle_Z[c->type];
	Kijk *= powl(rijeq + rjkeq, -1.);
	Kijk *= powl(phi0, -2.);
	Kijk *= expl(-2 * D);
	return Kijk * powl(phi - phi0, 2.);
}

double vdw_energy(struct Atom *a, struct Atom *b, struct Vector *apos, struct Vector *bpos, struct Model *m) {
	double EDEP = sqrtl(m->vdw_EDEP[a->type] * m->vdw_EDEP[b->type]);
	double sigma = sqrtl(m->vdwR[a->type] * m->vdwR[b->type]);
	double distance = calc_distance(apos, bpos);
	double term = sigma / distance;
	double term6 = term * term * term * term * term * term;
	double term12 = term6 * term6;
	
	return EDEP * (term12 - 2 * term6);
}

double dihedral_energy(struct Atom *a, struct Atom *b, struct Atom *c, struct Atom *d, struct Vector *apos, struct Vector *bpos, struct Vector *cpos, struct Vector *dpos, struct Model *m) {
    (void)a;
    (void)d;
    double omega = calc_omega(apos, bpos, cpos, dpos);
	
	// energy = pref * (1 + cos(n omega - gamma))
	double energy = m->dihedral_n[b->type][c->type] * omega - m->dihedral_gamma[b->type][c->type];
	energy = cos(energy);
	energy += 1;
	energy *= m->dihedral_pref[b->type][c->type];
	return energy;
}

double calc_energy(struct Molecule *m, unsigned int atom_offset, struct Vector * offset) {
	unsigned int i,j,k,l;
	struct Atom *a,*b,*c,*d;
	struct Vector zero;
	zero.x = 0; zero.y = 0; zero.z = 0;
	struct Vector apos, bpos, cpos, dpos;
	double energy = 0;
	double r;
	for (i = 0; i < m->n_atoms; i++) {
		a = &(m->as[i]);
		add_vector(&(a->v), (i == atom_offset)?offset:&(zero), &apos);
		for (j = 0; j < a->n_bonds; j++) {
			// Ebond
			b = a->bonds[j];
			add_vector(&(b->v), (b->i == atom_offset)?offset:&(zero), &bpos);
			energy += bond_energy(a, b, &apos, &bpos, m->model) / 2.; // divided by 2 as will be A-B and B-A
			for (k = j+1; k < a->n_bonds; k++) {
				// Eangle
				c = a->bonds[k];
				add_vector(&(c->v), (c->i == atom_offset)?offset:&(zero), &cpos);
				energy += angle_energy(b, a, c, &bpos, &apos, &cpos, m->model); // not divided as we only do this path once.
				for (l = 0; l < c->n_bonds; l++) {
					// Edihedral
					d = c->bonds[l];
					add_vector(&(d->v), (d->i == atom_offset)?offset:&(zero), &dpos);
					energy += dihedral_energy(d, c, a, b, &dpos, &cpos, &apos, &bpos, m->model)/2.; // ditto




				}
			}
		}
		for (j = i+1; j < m->n_atoms; j++) {
		
			if (i == j)
				continue;
		
			b = &(m->as[j]);
			add_vector(&(b->v), (b->i == atom_offset)?offset:&(zero), &bpos);
			
			
			//Evdw
			
			energy += vdw_energy(a, b, &apos, &bpos, m->model);
			
			
			
			// Eelectrostatic
			// Ees = ke q Q / r
			double q, Q;
			r = calc_distance(&apos, &bpos);
			q = m->model->charge[a->type];
			Q = m->model->charge[b->type];
			energy += COULOMB_CONSTANT * q * Q / r;

		}
	}
	return energy;
			
}

double kcal_to_atomic(double v) {

	return v * KCAL_ATOMIC;
}

double atomic_to_kcal(double v) {
	return v / KCAL_ATOMIC;
}

void set_temp_vel(struct Molecule *m, int atomid, double temp) {
	struct Vector v;

	v.x = ((rand()%100 - 49)/100.);
	v.y = ((rand()%100 - 49)/100.);
	v.z = ((rand()%100 - 49)/100.);
	if (v.x == 0 && v.y == 0 && v.z == 0)
		v.x = 1;
	normalise(&v);
	


	double mass = m->model->vdwM[m->as[atomid].type];

	double temp_factor = sqrtl(2 * BOLTZMANN_K * temp / mass);
	//if (atomid == 1)
		//printf("%lf, %lf, %lf :: %lf\n", v.x, v.y, v.z, temp_factor);
	m->as[atomid].vel.x = v.x * temp_factor;
	m->as[atomid].vel.y = v.y * temp_factor;
	m->as[atomid].vel.z = v.z * temp_factor;
}

double calc_temperature(struct Molecule *m) {
	double v;
	//T = m v^2 / (2k)
	unsigned int i;
	double sum = 0, mass;
	for (i = 0; i < m->n_atoms; i++) {
		v = magnitude(&(m->as[i].vel));
		mass = m->model->vdwM[m->as[i].type];
		sum += (mass * v * v) / (2 * BOLTZMANN_K);
	}
	//printf("Temp: %f K\n", sum/m->n_atoms);
	return sum / m->n_atoms;
}


void process_atom(int id, struct Molecule *m, double dn, double dt, double viscosity, double T) {

	double Tr = T; //calc_temperature(m);
	struct Vector o1, o2;
	struct Vector e_grad;
	double e1, e2;
	set_vector(&o1, dn, 0, 0);
	set_vector(&o2, -dn, 0, 0);
	e1 = kcal_to_atomic(calc_energy(m, id, &o1));
	e2 = kcal_to_atomic(calc_energy(m, id, &o2));
	e_grad.x = (e1 - e2) / (2 * dn);
	//printf("e1: %f, e2: %f, dn: %f, grad: %f\n", e1, e2, dn, e_grad.x);

	set_vector(&o1, 0, dn, 0);
	set_vector(&o2, 0, -dn, 0);
	e1 = kcal_to_atomic(calc_energy(m, id, &o1));
	e2 = kcal_to_atomic(calc_energy(m, id, &o2));
	e_grad.y = (e1 - e2) / (2 * dn);

	set_vector(&o1, 0, 0, dn);
	set_vector(&o2, 0, 0, -dn);
	e1 = kcal_to_atomic(calc_energy(m, id, &o1));
	e2 = kcal_to_atomic(calc_energy(m, id, &o2));
	e_grad.z = (e1 - e2) / (2 * dn);


	double mass = m->model->vdwM[m->as[id].type];

	struct Vector accel;
	// potential term
	double t3;
	//t1 = -e_grad.x;
	accel.x = -e_grad.x;
	accel.y = -e_grad.y;
	accel.z = -e_grad.z;
	// energy units are amu A^2 ps^-2, so units of gradient are amu A ps^-2.
	if (m->thermostat == LANGEVIN) {
		//t2 = -viscosity * mass * m->as[id].vel.x;
		accel.x -= viscosity * mass * m->as[id].vel.x;
		accel.y -= viscosity * mass * m->as[id].vel.y;
		accel.z -= viscosity * mass * m->as[id].vel.z;
		// not sure about viscosity units. reported to be ps^-1, but there's no mass term here.
		// so amu(?) ps^-1 A ps^-1 = amu A ps^-2

		double variance = 2 * mass * viscosity * BOLTZMANN_K * Tr / dt;
		// units are amu * ps^-1 * amu A^2 ps^-2 / K * K / ps
		//    = amu^2 * ps^-4 * A^2
		
		double std = sqrtl(variance);
		// so the units of the std are amu A ps^-2
		
		t3 = norm_rand(0, std);
		accel.x += t3; 
		accel.y += norm_rand(0, std);
		accel.z += norm_rand(0, std);
	}
	m->as[id].vel.x += dt * accel.x / mass;
	m->as[id].vel.y += dt * accel.y / mass;
	m->as[id].vel.z += dt * accel.z / mass;
	

	
	if (m->thermostat == ANDERSEN) {
		// with Andersen thermostat, the probability of a collision with bath is given as
		//  P(t) = v exp(-v t)
		// Giving a proportion of molecules which will have collided in that time. 
		// This proportion of atoms are then assigned random velocity based on the temperature.
		// Though v in this case is rate of collision, we reuse 'viscosity'.
		unsigned int k = 0;
		for (k = 0; k < m->n_atoms; k++)
			m->as[k].andersen_f = 0;
		double P = 1000. * viscosity * expl(-viscosity * dt);
		int Pr = (int) P;
		int ran = (rand()%1000);
		//printf("%d, %d -> %s\n", Pr, ran, (ran < Pr)?"yes":"no");
		if (ran < Pr) {
			//printf("%d, %d -> ", Pr, ran);
			set_temp_vel(m, id, T);
		}
		
		/*while (n_vel > 0) {
			k = (rand()%m->n_atoms);
			if (m->as[k].andersen_f == 1)
				continue;
			n_vel--;
			m->as[k].andersen_f = 1;
			printf("%d %d\n", k, n_vel);
			set_temp_vel(m, k, T);
		}*/
		
	}
	// amu A ps^-2 * ps / amu = A ps^-1
	return;
}
	
/* bond_xyz()
 *  XYZ files do not contain bond information; this function adds bonds
 *  between all atoms where the distance between the atoms is < bl.
 * arg *m: Molecule pointer
 * arg bl: Maximum bond length
 * TODO: Make bl able to vary depending on atoms involved 
 *       e.g., C-H bond shorter than C-C.
 */
void bond_xyz(struct Molecule *m, float bl) {
	unsigned int i, j;
	struct Vector temp;
	for (i = 0; i < m->n_atoms; i++) {
		for (j = i+1; j < m->n_atoms; j++) {
			// calculate magnitude of vector joining them
			sub_vector(&(m->as[i].v), &(m->as[j].v), &temp);
			if (magnitude(&temp) < bl) {
				// and if below bl, add a bond.
				add_bond(&(m->as[i]), &(m->as[j]));
				printf("Bonding %d - %d\n", i, j);
			}
		}
	}
	return;
}



/* read_xyz()
 *  Reads xyz file *filename, and generates atomic positions and puts
 *  them into molecule *m. 
 * arg *m: Molecule pointer
 * arg *filename: Filename.
 * returns: -1 if file does not exist or is wrong format.
 *           1 if okay.
 */
int read_xyz(struct Molecule *m, char *filename) {
	FILE * fp;
	char line[255];
	size_t len=255;
	printf("Reading xyz\n");
	fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("%s not found.\n", filename);
		return -1;
	}
	
	unsigned int c_line = 0;
	struct Atom *ca;
	while (fgets(line, len, fp)) {
		if (c_line == 0) {
			// first line in an xyz file gives the number of atoms.
			int k = sscanf(line, "%d\n", &(m->n_atoms));
			// if there is no number, close and exit.
			if (k != 1) {
				fclose(fp);
				return -1;
			} 
			// allocate memory for n_atoms atoms.
			m->as = (struct Atom *) malloc(sizeof(struct Atom) * m->n_atoms);
			if (m->as == NULL)
				crash(m);
		} else if (c_line >= 2) {
			// c_line = 1 is a comment line, which we ignore.
			// if the current atom (c_line - 2) is greater than n_atoms - 1, 
			// then exit. eg if there are 3 atoms, then as[0-2].
			if (c_line - 2 > m->n_atoms - 1) {
				fclose(fp);
				return -1;
			}
			// ca is pointer to the current atom being operated on.
			ca = &(m->as[c_line - 2]);
			ca->i = c_line-2;
			char *token = strtok(line, " \t");
			int i = 0;
			float x=0, y=0, z=0;
			char name[3];
			// while tokens are left...
			// i gives the column. First column is atom symbol, then x,y,z.
			while (token) {
				switch (i) {
					case 0: sprintf(name, "%2s", token); break;
					case 1: x = atof(token); break;
					case 2: y = atof(token); break;
					case 3: z = atof(token); break;
					// if there are more than 4 columns, crash and close.
					default: fclose(fp); return -1; break;
				}
				i++;
				token = strtok(NULL, " \t");
			}
			// add atom, and print out setup.
			add_atom(ca, x, y, z, name);
			//printf("%s :: (%f, %f, %f) %d\n", ca->name, ca->v.x, ca->v.y, ca->v.z, c_line-2);
		}
		c_line++;
	}
	fclose(fp);
	return 1;
}


/* save_xyz()
 *  Outputs *m as an xyz file in *filename. Tab delimited, x y z are of 
 *  format +00.00000.
 * returns: -1 if file does not open
 *           1 otherwise
 */
int save_xyz(struct Molecule *m, char *filename, char * mode) {
	FILE * fp;
	fp = fopen(filename, mode);
	if (fp == NULL)
		return -1;
	double T = calc_temperature(m);
	fprintf(fp, "%d\n%lf\n", m->n_atoms, T);

	unsigned int i;
	for (i = 0; i < m->n_atoms; i++) {
		fprintf(fp, "%2s\t%-2.6f\t%-2.6f\t%-2.6f\n", m->as[i].name, \
			m->as[i].v.x,\
			m->as[i].v.y,\
			m->as[i].v.z);
	}
	
	fclose(fp);
	return 1;
}

/* print_dir()
 *  Prints molecule to screen
 */
int print_dir(struct Molecule *m) {
	unsigned int i;
	for (i = 0; i < m->n_atoms; i++) {
		printf("%d\t%2s\t%-2.6f\t%-2.6f\t%-2.6f\n", i, m->as[i].name, \
			m->as[i].v.x,\
			m->as[i].v.y,\
			m->as[i].v.z);
	}
	return 0;
}

/*int assign_model(struct Molecule *m, char model_name[255]) {
	struct Model mod;
	if (strcmp(model_name, "MM3") == 0)
		setup_model(&mod, MOD_MM3);
	else {
		printf("Model does not exist.\n");
		crash(m);
	}
	
	m->model = &mod;
	return 1;
	
}*/

#define MINIMIZE 0
#define HEAT 1
#define PROD 2

void add_vector_sc(struct Vector *a, struct Vector *b, double dt) {
	a->x += (b->x * dt);
	a->y += (b->y * dt);
	a->z += (b->z * dt);
}

void propagate(struct Molecule *m, double dn, double dt, unsigned int steps, unsigned int output_step, char fn[500], double temp, double viscosity, int mode) {
	unsigned int i, k;

	//void process_atom(int id, struct Molecule *m, double dn, double dt, double viscosity, double T) {

	if (mode != HEAT)
		save_xyz(m, fn, "w");
	for (i = 0; i < steps; i++) {
		#pragma omp parallel for
		for (k = 0; k < m->n_atoms; k++) {
			if (mode != HEAT)
				process_atom(k, m, dn, dt, viscosity, temp);
			add_vector_sc(&(m->as[k].v), &(m->as[k].vel), dt);
			if (mode == MINIMIZE)
				set_vector(&(m->as[k].vel), 0, 0, 0);

		}
		if (i % output_step == 0)
			save_xyz(m, fn, "a");
	}
}

void minimize(struct Molecule *m, double dn, double dt, int steps, int output_step, char fn[500]) {
	propagate(m, dn, dt, steps, output_step, fn, 0, 0, MINIMIZE);
}


void iterate_once(struct Molecule *m, double dn, double dt, char fn[500], double temp, double viscosity) {
	propagate(m, dn, dt, 1, 1, fn, temp, viscosity, HEAT);
}


void heat(struct Molecule *m, double temp, double temp_step, int steps, double dn, double dt, double dnM, double dnT, double viscosity, char fn[500]) {
	unsigned int i;
	double T = 0;
	for (T = 0; T < temp + temp_step; T += temp_step) {
		for (i = 0; i < m->n_atoms; i++) {
			set_temp_vel(m, i, T);
			iterate_once(m, dn, dt, fn, T, viscosity);
			minimize(m, dnM, dnT, steps, 100, fn);
		}
	}
	return;
}

void production(struct Molecule *m, double dn, double dt, int steps, int output_step, char fn[500], double temp, double viscosity) {
	propagate(m, dn, dt, steps, output_step, fn, temp, viscosity, PROD);
}

/* run_script()
 *  Runs script in file *filename on molecule *m
 * Commands;
 *  - open FILENAME
 *     Opens xyz file at FILENAME
 *  - bond BL
 *     Generates bonds between atoms where distance < BL (angstroms)
 *  - rotate A B THETA
 *     Rotates atoms connected to B about the A-B axis an amount theta
 *     (radians)
 *  - output FILENAME
 *     Writes output as xyz file to FILENAME.
 * returns: -1 if file not available
 *           1 on success
 */
int run_script(char *filename, struct Molecule *m) {
	FILE * fp;
	char line[255];
	size_t len=255;
	fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("%s not found.\n", filename);
		return -1;
	}
	
	float bl;
	char command[255];
	int c;
	struct Vector zero;
	double dn, dt, temp, viscosity, temp_step, dnM, dnT;
	int steps, output_step;
	zero.x = 0; zero.y = 0; zero.z = 0;
	struct Model mod;
	while (fgets(line, len, fp)) {
		// reset system checks
		reset_check(m);
        // read in line, up to point c
		sscanf(line, "%255s %n", command, &c);
		if (command[0] == '%')
			continue;
		if (strcmp(command, "open") == 0) {
			// if there is no second argument, complain.
			if (sscanf(line+c, "%255s", command) != 1) {
				printf("error reading script\n%s", line);
				break;
			}
			read_xyz(m, command);
		} else if (strcmp(command, "bond") == 0) {
			if (sscanf(line+c, "%f", &bl) != 1) {
				printf("Error reading script\n%s", line);
				break;
			}
			bond_xyz(m, bl);
		} else if (strcmp(command, "output") == 0) {
			if (sscanf(line+c, "%255s", command) != 1) {
				printf("Error reading script\n%s", line);
				break;
			}
			save_xyz(m, command, "w");
		} else if (strcmp(command, "model") == 0) {
			if (sscanf(line + c, "%255s", command) != 1) {
				printf("Error reading script\n%s", line);
				break;
			}
			if (strcmp(command, "GAFF") == 0)
				setup_model(&mod, MOD_GAFF);
			m->model = &mod;
			
 		} else if (strcmp(command, "hybridize") == 0) {
 			//printf("%s\n", command);
 			if (sscanf(line+c, "%d %s", &steps, command) != 2) {
 				printf("Error reading script\n%s", line);
 				break;
 			}
 			//printf("%d steps\n", steps);
 			if (strcmp(command, "SP") == 0)
 				m->as[steps].hybridization = SP1;
 			else if (strcmp(command, "SP2") == 0)
 				m->as[steps].hybridization = SP2;
 			else if (strcmp(command, "SP3") == 0)
 				m->as[steps].hybridization = SP3;
 			else
 				printf("Hybridization not recognised %s\n", command);
 				
			
		} else if (strcmp(command, "energy") == 0) {
			if (m->model == NULL)
				printf("No model assigned\n");
			else
				printf("Energy: %f kcal/mol\n", calc_energy(m, 0, &zero));
		} else if (strcmp(command, "temperature") == 0) {
			printf("Temperature: %f K\n", calc_temperature(m));
			
		} else if (strcmp(command, "minimize") == 0) {
			if (sscanf(line + c, "%lf %lf %d %d %s", &dn, &dt, &steps, &output_step, command) != 5) {
				printf("Error reading script\n%s", line);
				break;
			}
			if (steps < 0)
				break;

			minimize(m, dn, dt, steps, output_step, command);
		} else if (strcmp(command, "heat") == 0) {
			if (sscanf(line + c, "%lf %lf %d %lf %lf %lf %lf %lf %s", &temp, &temp_step, &steps, &dn, &dt, &dnM, &dnT, &viscosity, command) != 9) {
				printf("Error reading script\n%s", line);
				break;
			}
			heat(m, temp, temp_step, steps, dn, dt, dnM, dnT, viscosity, command);
		} else if (strcmp(command, "iterate") == 0) {
			if (sscanf(line + c, "%lf %lf %lf %lf %s", &dn, &dt, &temp, &viscosity, command) != 5) {
				printf("Error reading script\n%s", line);
				break;
			}
			if (m->model == NULL) {
				printf("No model assigned.\n");
				break;
			} 
			if (m->thermostat == -1) {
				printf("No thermostat assigned.\n");
				break;
			}
			
			iterate_once(m, dn, dt, command, temp, viscosity);
		} else if (strcmp(command, "thermostat") == 0) {
			if (sscanf(line + c, "%s", command) != 1) {
				printf("Error reading thermostat\n");
				break;
			}
			if (strcmp(command, "LANGEVIN") == 0) {
				printf("Please note the Langevin thermostat is broken.\n");
				m->thermostat = LANGEVIN;
			} else if (strcmp(command, "ANDERSEN") == 0) {
				m->thermostat = ANDERSEN;
			}
		} else if (strcmp(command, "prod") == 0) {
			if (sscanf(line + c, "%lf %lf %lf %lf %d %d %s", &dn, &dt, &temp, &viscosity, &steps, &output_step, command) != 7) {
				printf("Error reading script\n%s", line);
				break;
			}
			if (steps < 0)
				break;
			if (m->model == NULL) {
				printf("No model assigned.\n");
				break;
			} 
			if (m->thermostat == -1) {
				printf("No thermostat assigned.\n");
				break;
			}
			production(m, dn, dt, steps, output_step, command, temp, viscosity);
		}
	}
	fclose(fp);
	return 1;
}

/* main()
 *  Runs system. If a script is passed, runs the script. Else,
 *  enters a do-while loop operating on commands in an interactive manner
 * args SCRIPT: can pass script filename which will be run in.
 */
int main(int argc, char *argv[]) {
	struct Molecule mol;
	mol.n_atoms = 0;
	mol.model = NULL;
	mol.thermostat = -1;
	
	char command[255];
	int n;

	char script_name[255];
	if (argc >= 2) {
		strcpy(script_name, argv[1]);
		run_script(script_name, &mol);
	} else {
		/* Interactive mode commands
		 * - open FILENAME
		 *    opens xyz file FILENAME
		 * - bond BL
		 *    sets up bonds between all atoms with distance < BL angstroms
		 * - graph N
		 *    prints out graph starting at atom N
		 * - print
		 *    prints out entire system
		 * - rotate A B THETA
		 *    applies a rotation to atoms connected to B about the A-B axis
		 *    an amount theta radians.
		 * - output FILENAME
		 *    outputs as XYZ file into FILENAME
		 * - run SCRIPT
		 *    runs script in SCRIPT.
		 * - exit
		 *    exits program
		 */
		do {
			reset_check(&mol);
			printf("> ");
			scanf("%255s", command);
			printf("%s\n", command);
			if (strcmp(command, "open") == 0) {
				scanf("%255s", command);
				read_xyz(&mol, command);
			} else if (strcmp(command, "bond") == 0) {
				printf(" max bond length (ang) > ");
				scanf("%255s", command);
				bond_xyz(&mol, atof(command));
			} else if (strcmp(command, "graph") == 0) {
				printf(" start atom > ");
				if (scanf("%d", &n) == 1)
					print_moleculef(&(mol.as[n]), 0);
			} else if (strcmp(command, "print") == 0) {
				print_dir(&mol);

			} else if (strcmp(command, "output") == 0) {
				printf(" filename > ");
				if (scanf("%255s", command) != 1)
					continue;
				save_xyz(&mol, command, "w");

			} else if (strcmp(command, "run") == 0) {
				printf(" script > ");
				if (scanf("%255s", command) != 1)
					continue;
				run_script(command, &mol);
			}
		} while (strcmp(command, "exit") != 0);
	}
	// clean up everything
	free_atoms(&mol);
	free(mol.as);
	return 1;
}
