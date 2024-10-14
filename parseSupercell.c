#include "funcs.h"
#include <string.h>

/*This function gets all the parameters necessary to properly read the data from the gulp input
 * and to then properly create a bonded lammps input. It reads in the cell size in each dimension,
 * the spring constant for the harmonic pot, the equillibrium distance, and the cutoff for bonded atoms.
 */
int getCellSize(const char *datafile , double *cell )
{
	int natoms=0, i;
	char strang[80] ;
	double d ;

	FILE *fpoint = safe_open( datafile, "r" );

	//Get Lattice Constants
	SkipTo( fpoint, "cell", strang );
	fscanf(fpoint,"%lf %lf %lf",&cell[0],&cell[1],&cell[2]);

	//Determine Number of atoms.
	fpoint = SkipTo( fpoint, "frac",strang );
	while ( fscanf( fpoint, "%s %lf %lf %lf %d %d %d", strang,&d,&d,&d, &i,&i,&i ) == 7 )
		natoms++;

	fclose( fpoint );

  printf("natoms: %d\n",natoms);
  printf("cell:   %lf %lf %lf\n",cell[0],cell[1],cell[2]);

	return natoms;
}


void read_data (const char *datafile, double *x, double *y , double *z, int natoms, char **namespt)
{
int i, f;
char strang[80];

printf("\nread_data\n\n");

FILE *fpoint = safe_open( datafile, "r" );
fpoint = SkipTo( fpoint, "frac" , strang );

for(i=0;i<natoms;i++)
   fscanf( fpoint, "%s %lf %lf %lf %d %d %d", namespt[i],&x[i],&y[i],&z[i],&f,&f,&f );

printf("\n\n");
fclose( fpoint );

}
