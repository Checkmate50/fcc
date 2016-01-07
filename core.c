/*
 * Given a file with format:
 * m a
 * This programs generates a face-centered crystal
 * With m^3 cubes of side length a
 * Note that behavior is undefined for m, a <= 0
 *
 * Output is written to pos.pos[0]yz and is of the form:
 * N
 * H x_1 y_1 z_1
 * ...
 * H x_n y_n z_n
 *
 * Where N is the number of molecules
 * And x_i, y_i, and z_i are the coordinates of the i-th molecule
 *
 * Written by Dietrich Geisler
 * Last Modified 12/27/15
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define ARC4RANDOM_MAX     0x100000000

typedef struct {
    double* p; //position (3-d vector)
    double* f; //force (3-d vector)
    double* v; //velocity (3-d vector)
    double mass;
} Molecule;

void create_fcc(Molecule* data, int m, double a);
void zero_momentum(Molecule* molecules, int N);
double compute_en(Molecule* molecules, double eps, double sig, int index, int size, double L);
double compute_pair_en(Molecule p1, Molecule p2, double eps, double sig, double* f1, double* f2, double L);
double md_step(Molecule* molecules, double eps, double sig, int index, int size, double L, double delta_t);

int main(int argc, char* argv[]) {
    if (argc != 2) {
	printf("Please provide a file.\n");
	return 0;
    }
    int N, m;
    double a;
    FILE* ifile,* ofile,* efile;
    Molecule* molecules;
    double energy, kinetic, v;
    int i, j, k, count;
    srand(time(NULL));
    
    //Read in data
    ifile = fopen(argv[1], "r");
    fscanf(ifile, "%i %lf", &m, &a);
    fclose(ifile);
    
    N = 4*m*m*m; //4m^3
    molecules = (Molecule*) malloc(N*sizeof(Molecule));
    for (i = 0; i < N; i++) {
	molecules[i].p = malloc(3*sizeof(double));
	for (j = 0; j < 3; j++) {
	  molecules[i].p[j] = 0;
	}
	molecules[i].mass = 1;
    }
    create_fcc(molecules, m, a);

    //Output the result
    ofile = fopen("pos.xyz", "w");
    fprintf(ofile, "%i\n\n", N);
    for (i = 0; i < N; i++) {
	fprintf(ofile, "C\t%f\t%f\t%f\n", molecules[i].p[0], molecules[i].p[1], molecules[i].p[2]);
    }
    fclose(ofile);
    
    /*
    Molecule mol;
    ofile = fopen("data.txt", "w");
    for (i = 0; i < N; i++) {
	molecules[i].f = malloc(3*sizeof(double));
    }
    for (i = 0; i < N; i++) {
	energy = compute_en(molecules, 1, 1, i, N);
	mol = molecules[i];
	//fprintf(stderr, "%i:\t%f\t(%f, %f, %f)\n", i+1, energy, mol.f[0], mol.f[1], mol.f[2]);
	fprintf(stderr, "%f\t%f\n", a, energy);
	}*/
    
    //Writes the energy data as we vary a from .9 to 3 to file
    /*double energy;
    m = 5;
    N = 500;
    molecules = (Molecule*) malloc(N*sizeof(Molecule));
    ofile = fopen("data.txt", "w");
    for (a = .9; a < 3; a += .1) {
	create_fcc(molecules, m, a);
	energy = 0;
	for (i = 0; i < N; i++) {
	    molecules[i].f = malloc(3*sizeof(double));
	}
	for (i = 0; i < N; i++) {
	    energy += compute_en(molecules, 1, 1, i, N, a*m);
	}
	fprintf(ofile, "%f\t%f\n", a, energy/N);
    }*/

    //Writes the energy data as we vary N from 4 to 13500 to file
    /*
    a = 2;
    ofile = fopen("data.txt", "w");
    for (m = 1; m <= 15; m++) {
	N = 4*m*m*m;
	molecules = (Molecule*) malloc(N*sizeof(Molecule));
	create_fcc(molecules, m, a);
	energy = 0;
	for (i = 0; i < N; i++) {
	    molecules[i].f = malloc(3*sizeof(double));
	}
	for (i = 0; i < N; i++) {
	    energy += compute_en(molecules, 1, 1, i, N, a*m);
	}
	fprintf(ofile, "%i\t%f\n", N, energy/N);
	}*/

    for (i = 0; i < N; i++) {
	molecules[i].f = malloc(3*sizeof(double));
	molecules[i].v = malloc(3*sizeof(double));
	for (j = 0; j < 3; j++) {
	    //molecules[i].v[j] = 0;
	    molecules[i].v[j] = ((double)rand() / (RAND_MAX/2)) - 1; //Set velocity kinda randomly	    
	}
    }
    
    zero_momentum(molecules, N);
    for (i = 0; i < N; i++) {
	compute_en(molecules, 1, 1, i, N, a*m);
    }
    efile = fopen("energy.txt", "w");
    ofile = fopen("pos_t1.xyz", "w");
    char filename[20];
    for (count = 0; count < 10; count++) {
	snprintf(filename, 20, "pos_t%i.xyz", count+1);
	ofile = fopen(filename, "w");
	fprintf(ofile, "%i\n\n", N);
	for (i = 0; i < 1000; i++) {
	    energy = 0;
	    kinetic = 0;
	    for (j = 0; j < N; j++) {
		energy += md_step(molecules, 1, 1, j, N, a*m, .001);
		v = 0;
		for (k = 0; k < 3; k++) {
		    v += pow(molecules[j].v[k], 2);
		}
		kinetic += molecules[j].mass*v/2;
	    }
	    if (i % 100 == 0) {
		fprintf(efile, "%f\t%f\t%f\n", energy, kinetic, energy+kinetic);
	    }
	}
	for (j = 0; j < N; j++) {
	    fprintf(ofile, "C\t%f\t%f\t%f\n", molecules[j].p[0], molecules[j].p[1], molecules[j].p[2]);
	}
	printf("%i\n", count);
	fclose(ofile);
    }
    fclose(efile);
    for (i = 0; i < N; i++) {
	free(molecules[i].p);
	free(molecules[i].v);
	free(molecules[i].f);
    }
    free(molecules);
}

/*
 * Fills the given 4m^3 sized data array with molecular information
 * Assumes that molecules are arranged as an FCC
 */
void create_fcc(Molecule* data, int m, double a) {
    int i, j, k, I , J;
    double b = a / 2; //For equation simplicity
    int N = 4*m*m*m;
    //Below, we iterate through each "cell" in the crystal
    //Observe that we are we creating a cube of cells
    //The array access patterns allows us to have 4 atoms per cell
    //In each cell, the atoms are at the points given on page 2 of our document
    for (i = 0; i < m; i++) {
	for(j = 0; j < m; j++) {
	    for (k = 0; k < m; k++) {
		I = m*m*i;
		J = m*j;
		data[4*(I+J+k)].p[0] = i*a;
		data[4*(I+J+k)].p[1] = j*a;
		data[4*(I+J+k)].p[2] = k*a;
		data[4*(I+J+k) + 1].p[0] = b+i*a; 
		data[4*(I+J+k) + 1].p[1] = b+j*a; 
		data[4*(I+J+k) + 1].p[2] = k*a;
		data[4*(I+J+k) + 2].p[0] = b+i*a; 
		data[4*(I+J+k) + 2].p[1] = j*a;
		data[4*(I+J+k) + 2].p[2] = b+k*a;
		data[4*(I+J+k) + 3].p[0] = i*a; 
		data[4*(I+J+k) + 3].p[1] = b+j*a; 
		data[4*(I+J+k) + 3].p[2] = b+k*a;
	    }
	}
    }
}

/*
 * Sets the center of mass momentum for the given molecules to zero
 */
void zero_momentum(Molecule* molecules, int N) {
    int i, j;
    double M;
    double* VCM;
    M = 0;
    VCM = malloc(3*sizeof(double));
    for (i = 0; i < 3; i++) {
	VCM[i] = 0;
    }
    for (i = 0; i < N; i++) {
	M += molecules[i].mass;
	for (j = 0; j < 3; j++) {
	    VCM[j] += molecules[i].mass*molecules[i].v[j];
	}
    }
    for (i = 0; i < N; i++) {
	for (j = 0; j < 3; j++) {
	    molecules[i].v[j] -= M*VCM[j]/(N*molecules[i].mass);
	}
    }
}

/*
 * Computes the total energy and forces acting on the molecule given by index
 * Uses the list of molecules, excluding the given index
 * Stores the forces acting on the given molecule in the associated 'f'
 * Returns the total potential energy acting on the molecule
 */
double compute_en(Molecule* molecules, double eps, double sig, int index, int size, double L) {
    Molecule m = molecules[index];
    double* force;
    double* temp_force;
    double energy, e;
    int i, j;
    force = malloc(3*sizeof(double));
    temp_force = malloc(3*sizeof(double));
    energy = 0;
    for (i = 0; i < 3; i++) {
	m.f[i] = 0; //Reset force for calculation
    }
    for (i = 0; i < size; i++) {
	if (i == index) {
	    continue;
	}
	e = compute_pair_en(m, molecules[i], eps, sig, force, temp_force, L);
	if (i > index) {
	    energy += e;
	}
	for (j = 0; j < 3; j++) {
	    m.f[j] += force[j];
	}
    }
    return energy;
}

/*
 * Computes the energy and forces between the two given molecules
 * Stores the forces on each molecule within f1 and f2 respectively
 * Returns the total potential energy of this interaction
 */
double compute_pair_en(Molecule p1, Molecule p2, double eps, double sig, double* f1, double* f2, double L) {
    double* pos2 = malloc(3*sizeof(double));
    double r, X, Y, Z;
    double a, b;
    int i;
    for (i = 0; i < 3; i++) {
	pos2[i] = p2.p[i];
    }
    if (L > 0) {
	for (i = 0; i < 3; i++) {
	    if ((p1.p[i] - p2.p[i]) < -L/2) {
		pos2[i] -= L;
	    }
	    else if ((p1.p[i] - p2.p[i]) > L/2) {
		pos2[i] += L; 
	    }
	}
    }
    X = p1.p[0]-pos2[0];
    Y = p1.p[1]-pos2[1];
    Z = p1.p[2]-pos2[2];
    r = sqrt(X*X+Y*Y+Z*Z);
    

    //Fill in values for forces
    a = (12*pow(sig, 12)*(pos2[0]-p1.p[0]))/pow(r, 14);
    b = (6*pow(sig, 6)*(pos2[0]-p1.p[0]))/pow(r, 8);
    f1[0] = -4*eps*(a-b);
    f2[0] = -4*eps*(b-a);

    a = (12*pow(sig, 12)*(pos2[1]-p1.p[1]))/pow(r, 14);
    b = (6*pow(sig, 6)*(pos2[1]-p1.p[1]))/pow(r, 8);
    f1[1] = -4*eps*(a-b);
    f2[1] = -4*eps*(b-a);

    a = (12*pow(sig, 12)*(pos2[2]-p1.p[2]))/pow(r, 14);
    b = (6*pow(sig, 6)*(pos2[2]-p1.p[2]))/pow(r, 8);
    f1[2] = -4*eps*(a-b);
    f2[2] = -4*eps*(b-a);

    return 4*eps*(pow(sig/r, 12.0)-pow(sig/r, 6.0));
}

/*
 * Given a molecule from a list with a calculated force, velocity, and position,
 * Calculates the new position, velocity, and force of the given molecule
 * After the given time of delta_t has passed
 * These new values are stored in the molecule itself
 * Returns the potential energy of the given molecule at time t + delta_t
 */
double md_step(Molecule* molecules, double eps, double sig, int index, int size, double L, double delta_t) {
    Molecule m;
    double energy;
    int i;
    m = molecules[index];
    //Calculate the new position
    for (i = 0; i < 3; i++) {
	m.p[i] = m.p[i] + m.v[i]*delta_t + m.f[i]*(pow(delta_t, 2)/(2*m.mass));
	if (m.p[i] >= L) {
	    m.p[i] -= L;
	}
	else if (m.p[i] < 0) {
	    m.p[i] += L;
	}
    }
    //Calculate half-step velocity
    for (i = 0; i < 3; i++) {
	m.v[i] = m.v[i] + m.f[i]*(delta_t/(2*m.mass));
    }

    //Calculate the new forces
    energy = compute_en(molecules, eps, sig, index, size, L);

    //Calculate the new velocities
    //printf("index: %i\n\tPosition\tVelocity\tForce\n", index);
    for (i = 0; i < 3; i++) {
	m.v[i] = m.v[i] + m.f[i]*(delta_t/(2*m.mass));
//	printf("%c:\t%f\t%f\t%f\n", (char) ((int) 'X') + i, m.p[i], m.v[i], m.f[i]);
    }

    return energy;
}
