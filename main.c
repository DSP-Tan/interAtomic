/*This code reads in a gulp .gin input file and converts the pertinent atomic information
 * to lammps format. It also performs a few calculations of the initial values of certain thermodynamic
 * variables. The code also converts the gulp harmonic potentials with cutoffs to permanently bonded lammps
 * potentials */

//TO Do: I assume that in the initialisation of the anglecoeffs array below
//that 1.949, the equillibrium bond length for GaN, is only a temporary value,
//to be later replaced with the equillibrium bond length between whatever atoms
//happen to be in question. You must verify this, or there are serious problems.

//To Do: Put in safety measures to ensure that the supercell file and potential
//files have the correct formats so that the code doesn't read in bad numbers and
//then cause a crash somewhere.

//To Do: In read data, you are not assigning namespt in an efficient manner. Rectify this,
//you should increment the pointer, rather than start at the start of the pointer and count up from it.


//To Do: Check how the angle coefficients are written. Check their initialisation, try to understanding whatever it
//is that you were doing when you said print it if i==0. what's that about?

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "funcs.h"

// int getCellSize ( char datafile[] , double *cell );
// void read_data	( char datafile[] , double *x, double *y , double *z, int natoms , char **namespt );
int rand_lim    ( int limit) ;

int main( int argc, char** argv){

if (argc!=3)
{
printf("\nProgram Called Incorrectly");
printf("\n\nusage:\n%s [supercell file with cell parameters at top] [Potentials File] \n\n",argv[0]);
printf("The supercell file should have the format:\ncell\n");
printf("%%lf %%lf %%lf - These are the a b c lattice constants\n");
printf("frac\n");
printf("%%s %%lf %%lf %%lf\n");
printf("%%s %%lf %%lf %%lf\n");
printf("%%s %%lf %%lf %%lf\n");
printf("...\n");
printf("%%s %%lf %%lf %%lf - name x y z\n\n");
return -1;
}

printf("The supercell file should have the format:\ncell\n");
printf("%%lf %%lf %%lf - These are the a b c lattice constants\n");
printf("frac\n");
printf("%%s %%lf %%lf %%lf\n");
printf("%%s %%lf %%lf %%lf\n");
printf("%%s %%lf %%lf %%lf\n");
printf("...\n");
printf("%%s %%lf %%lf %%lf - name x y z\n\n");


FILE *fout1;
int natoms,i,l, **ptr, *type;			                          //natoms is number of atoms
double cutoff, cell[3], *charge;						                //cutoff is the Bonded potentials cutoff, this is used to determine bond and angle topology. It is 2.4 in the case of the model used here. The array cell contains the cell size in x,y and z. The pointer charge points to the memory location of the charges on each atom.
const char *dest,*inputCellFile, *Potentials;	              //This will be arrays containing the names of different files, two output files and two input files.
double *x, *z, *y, *x_, *y_, *z_; 						              // The x, y, z are the scaled lammps initial coords. The x_,y_,z_ are the unscaled coords obtained by multiplying the scaled coords by the lattice vectors.
int mol_ID = 1; 										                        //This is the molecule I.D needed in the lammps data file which will always be 1.
char **namespt,**useptt;

cutoff = 2.4;

/*bondptr is a pointer to an array of all bonds in the system
 *bondcoeffs[96] is an array which will contain all the bond types, 96 being
 *the maximimum number of different bonds in a particular configuration.
 *and coefficients which characterise them. angle ptr/coeffs are same
 *but for angles. */

struct bond *bondptr, bondcoeffs[96];
int NumberOfBonds;
int NumberOfBondTypes;
struct Angles *angleptr, anglecoeffs[9];
int NumberOfAngles;
int NumberOfAngleTypes;

for(i=0;i<96;i++){
   bondcoeffs[i].atType1 = 0;
   bondcoeffs[i].atType2 = 0;
   bondcoeffs[i].atID1 = 0;
   bondcoeffs[i].atID2 = 0;
   bondcoeffs[i].K1 = 0;
   bondcoeffs[i].K2 = 0;
   bondcoeffs[i].r0 = 0;
   bondcoeffs[i].BType = 0;
   bondcoeffs[i].B_ID = 0;
   bondcoeffs[i].present = 0;
   }
for(i=0;i<9;i++){
   anglecoeffs[i].atom1 ='\0';
   anglecoeffs[i].atom2 ='\0';
   anglecoeffs[i].atom3 ='\0';
   anglecoeffs[i].K2 =0;
   anglecoeffs[i].Theta0 = 109;
   anglecoeffs[i].M = 0;
   anglecoeffs[i].N1 = 0;
   anglecoeffs[i].N2 = 0;
   anglecoeffs[i].r1 = 1.949;
   anglecoeffs[i].r2 = 1.949;
   anglecoeffs[i].present = 0;
   anglecoeffs[i].at1 = 0;
   anglecoeffs[i].at2 = 0;
   anglecoeffs[i].at3 = 0;
   }

/*------------------------------------------------------------------------------------*/
/*-----------------OutPut and Input Files Definition----------------------------------*/
/*------------------------------------------------------------------------------------*/

inputCellFile=argv[1];       //The first argument to the funciton is the file containing all of the atomic positions and lattice constants.
Potentials   =argv[2];       //The second argument to the funciton is the file containing all the potential parameters.
dest         ="data.lammps"; //This is the lammps data file which will contain bond and angle topology and be input to lammps script.

puts(inputCellFile);
puts(dest);

/*------------------------------------------------------------------------------------------------------*/
/*----------Find number of atoms for memory allocation and variable declarations/definitions------------*/
/*------------------------------------------------------------------------------------------------------*/

natoms = getCellSize (inputCellFile , cell );

x      =( double *)malloc( natoms*sizeof(double) );
y      =( double *)malloc( natoms*sizeof(double) );
z      =( double *)malloc( natoms*sizeof(double) );
x_     =( double *)malloc( natoms*sizeof(double) );
y_     =( double *)malloc( natoms*sizeof(double) );
z_     =( double *)malloc( natoms*sizeof(double) );
type   =(   int*  )malloc( natoms*sizeof( int  ) );
charge =( double *)malloc( natoms*sizeof(double) );

useptt = namespt=(  char **)malloc( natoms*sizeof( char*) );
for(i=0;i<natoms;i++) {
   *useptt = (char*)malloc(8 * sizeof(char));
   memset(*useptt, '\0', 8);
   useptt++;
   }

printf("\n\nMemory Allocations Complete\n");

int num_types;								//This variable will hold the number of atom types.


/* With the number of atoms obtained by getCellSize, we can now see what size our arrays to contain the coords need to be
 * x,y and z are the fractional gulp coords which will be read in and x_,y_ and z_ are the absolute lammps coordinates which
 * will be output and used to calculate which atoms are bonded to which.
 */
read_data ( inputCellFile, x, y, z, natoms, namespt ) ;

//Get information on atom types, charges etc.
printf("\nAtom Types\n\n");
ptr = atomtypes( natoms, namespt , type, charge , &num_types );

//Get unscaled atomic coordinates
for(i=0;i<natoms;i++){
   x_[i] = x[i]*cell[0];
   y_[i] = y[i]*cell[1];
   z_[i] = z[i]*cell[2];
   }

//Get Bond and Angular information
printf("\n\nMakebonds\n\n");
bondptr = makebonds(Potentials,
                    natoms,x_,y_,z_,
                    ptr,type,num_types,namespt,cutoff,cell,
                    &NumberOfBonds,bondcoeffs,&NumberOfBondTypes);
printf("\n\nMake Angles \n\n");
angleptr = makeangles (Potentials,natoms,x_,y_,z_,namespt,cutoff,cell,bondptr,NumberOfBonds,&NumberOfAngles,anglecoeffs,&NumberOfAngleTypes);

//--------------Write General Information-------------------------------------------//

fout1 = fopen ( dest, "w" );
fprintf(fout1,"#Number of atoms, bonds and angles\n\n");
fprintf(fout1,"%d atoms\n",natoms);
fprintf(fout1,"%d bonds\n", NumberOfBonds);
fprintf(fout1,"%d angles\n\n", NumberOfAngles);

fprintf(fout1,"%d atom types\n",num_types);
fprintf(fout1,"%d bond types\n", NumberOfBondTypes);
fprintf(fout1,"%d angle types\n\n", NumberOfAngleTypes);

fprintf( fout1, "0.0\t%20.9lf xlo xhi\n", cell[0] );
fprintf( fout1, "0.0\t%20.9lf ylo yhi\n", cell[1] );
fprintf( fout1, "0.0\t%20.9lf zlo zhi\n\n", cell[2] );

fprintf( fout1, "Masses\n\n" );
   for(i=0;i<num_types;i++)
      fprintf( fout1, "%d\t1\n",i+1);

/*------------------------------Write Atomic Positions------------------------------------------------------*/
fprintf( fout1 ,"\n\nAtoms\t\n\n" );

for(i=0; i<natoms; i++)
   fprintf(fout1,"%3d %1d %2d\t %lf\t %20.9lf \t %20.9lf \t %20.9lf\n",i+1,mol_ID,type[i],charge[i], x_[i],y_[i],z_[i]);


/*------------------------------Write Bonding Information---------------------------------------------------*/
fprintf(fout1,"\n\nBonds\n\n");
for(i=0;i<NumberOfBonds;i++)
   fprintf( fout1, "%d\t%d\t%d\t%d\n", i+1, bondptr[i].BType, bondptr[i].atID1,bondptr[i].atID2);

/*------------------------------Write Angle Information And Perform Diagnostics-------------------------*/


fprintf(fout1,"\n\nAngles\n\n");
for(i=0;i<NumberOfAngles;i++)
   fprintf( fout1, "%d\t%d\t%d\t%d\t%d\n",i+1,angleptr[i].type,angleptr[i].at1,angleptr[i].at2,angleptr[i].at3 );

printf("\n\nThere are %d Angles\n\n", NumberOfAngles);



/*---------------------------------------------Write Coefficient Information---------------------*/
	fprintf(fout1,"\n\nBond Coeffs\n\n");
    for(i=0;i<NumberOfBondTypes;i++)
       fprintf( fout1, "%d\t%lf\t%lf\t%lf\t0.0\n" ,bondcoeffs[i].BType, bondcoeffs[i].r0, bondcoeffs[i].K2, bondcoeffs[i].K1);

	fprintf(fout1,"\n\nAngle Coeffs\n\n");
    l=anglecoeffs[0].type;
    for(i=0;i<9;i++){
       if(i==0 && anglecoeffs[0].present==1)
          fprintf(fout1,"%d\t%lf\t%lf\t%d\t%d\n",anglecoeffs[0].type,anglecoeffs[0].Theta0,anglecoeffs[0].K2/2,0,0);
       else if(l != anglecoeffs[i].type && anglecoeffs[i].present==1)
          fprintf(fout1,"%d\t%lf\t%lf\t%d\t%d\n",anglecoeffs[i].type,anglecoeffs[i].Theta0,anglecoeffs[i].K2/2,0,0);
       l = anglecoeffs[i].type;
       }

    fprintf(fout1,"\n\nBondBond Coeffs\n\n");
    l=anglecoeffs[0].type;
    for(i=0;i<9;i++){
       if(i==0 && anglecoeffs[0].present==1 )
          fprintf(fout1,"%d\t%lf\t%lf\t%lf\n",anglecoeffs[0].type,anglecoeffs[0].M,anglecoeffs[0].r1,anglecoeffs[0].r2);
       else if(l != anglecoeffs[i].type && anglecoeffs[i].present==1)
          fprintf(fout1,"%d\t%lf\t%lf\t%lf\n",anglecoeffs[i].type,anglecoeffs[i].M,anglecoeffs[i].r1,anglecoeffs[i].r2);
       l = anglecoeffs[i].type;
       }

    fprintf(fout1,"\n\nBondAngle Coeffs\n\n");
    l=anglecoeffs[0].type;
    for(i=0;i<9;i++){
       if(i==0 && anglecoeffs[0].present==1)
          fprintf(fout1,"%d\t%lf\t%lf\t%lf\t%lf\n",anglecoeffs[0].type,anglecoeffs[0].N1,anglecoeffs[0].N2,anglecoeffs[0].r1,anglecoeffs[0].r2);
       else if(l!=anglecoeffs[i].type && anglecoeffs[i].present==1)
          fprintf(fout1,"%d\t%lf\t%lf\t%lf\t%lf\n",anglecoeffs[i].type,anglecoeffs[i].N1,anglecoeffs[i].N2,anglecoeffs[i].r1,anglecoeffs[i].r2);
       l = anglecoeffs[i].type;
       }


return 0;
}








int rand_lim(int limit) {
/* return a random number between 0 and limit inclusive.
 */
//TO DO: Make it so that the random increment of coords is not just an increment, but a random adjustment between 1 and -1 % of the current coords.
//TO DO: It would be better to seed with the time so subsequent calls do not produce the same number.
	srand(2);
    int divisor = RAND_MAX/(limit+1);
    int retval;

    do {
        retval = rand() / divisor;
    } while (retval > limit);

    return retval;
}
