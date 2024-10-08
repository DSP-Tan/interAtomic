#include<stdio.h>

/*
struct Atoms{
int         Id;
int       Type;
double  charge;
double       x;
double       y;
double       z;
char     *name;
};
*/

struct Angles{
char atom1;				//These three chars are to identify the types of the atoms present.
char atom2;				//The extent to which to type of the atom affects the type of three body force that occurs between the three atoms
char atom3;				//is completely determined by the first letter of these atoms. (i.e. whether they are In, Ga, N, or Al);
double K2;				//These are the parameters for the class2 interactions, as describd on the angle_style class2 lammps document page.
double Theta0;
double M;
double N1;
double N2;
double r1;
double r2;
int present;
int type;
int at1;				//These integers are for if we want to completely reccord all the information for every threebody interaction in an
int at2;				//array of structs of this type.
int at3;
};

struct bond {
	int atType1;
	int atType2;
	int atID1;
	int atID2;
	double K1;
	double K2;
	double r0;
	int BType;
	int B_ID;
	int present;
};

// file parsing functions
int getCellSize( char datafile[] , double *cell );
void read_data ( char datafile[], double *x, double *y , double *z, int natoms, char **namespt);

double Distance ( double x1, double y1, double z1, double x2, double y2, double z2 );
int** atomtypes( int natoms,  char **namespt , int *type, double *charge, int *num_types );
struct bond* makebonds(char *source, int natoms, double *x_, double *y_, double *z_, int **types, int *type,int num_types, char **namespt ,double cutoff,double *cell,int *NumberOfBonds,struct bond *bondcoeffs, int *NumberOfBondTypes);
struct Angles* makeangles ( char *source,int natoms, double *x_,double *y_,double *z_, char **namespt, double cutoff,double *cell, struct bond *bonds,int NumBonds, int* NumberOfAngles, struct Angles *PosibAngs, int *NumberOfAngleTypes );
FILE* SkipTo( FILE *fptr , const char* command , char cmdline[] );
FILE* safe_open	( const char* filename , const char *mode );
