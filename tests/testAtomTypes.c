#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "../funcs.h"

int cmpfunc (const void * a, const void * b);
int **atomtypes_test( int     natoms, char **namespt , int *type,
                      double *charge, int *num_types );

int main (){

int    natoms,i;
int    **types,**types_test, *type, *type_test;
double cell[3], *charge,*charge_test;
char   *inputCellFile, *potFile;
double *x, *z, *y;
char   **namespt,**useptt;

potFile       = "../Pots";

const char *inputCellFiles[] = {"10_10_20_1.txt",
                                "../supercell/sampleInputSupercell.txt",
                                "../supercell/bigInGaNSampleInput.txt"};


for(int i_test=0; i_test <3 ; i_test ++){

   natoms = getCellSize (inputCellFiles[i_test] , cell );

   x           =( double *)malloc( natoms*sizeof(double) );
   y           =( double *)malloc( natoms*sizeof(double) );
   z           =( double *)malloc( natoms*sizeof(double) );
   type        =(   int*  )malloc( natoms*sizeof( int  ) );
   type_test   =(   int*  )malloc( natoms*sizeof( int  ) );
   charge      =( double *)malloc( natoms*sizeof(double) );
   charge_test =( double *)malloc( natoms*sizeof(double) );

   useptt = namespt=(  char **)malloc( natoms*sizeof( char*) );
   for(i=0;i<natoms;i++) {
      *useptt = (char*)malloc(8 * sizeof(char));
      memset(*useptt, '\0', 8); useptt++;
      }

   read_data ( inputCellFiles[i_test], x, y, z, natoms, namespt ) ;

   int num_types, num_types_test;
   types = atomtypes( natoms, namespt , type, charge , &num_types );
   types_test = atomtypes_test( natoms, namespt , type_test, charge_test , &num_types_test );

   printf("Running test #%d on file:\n%s", i_test, inputCellFiles[i_test] );
   assert(num_types==num_types_test);
   for(i=0;i<natoms;i++){
     assert(type[i]==type_test[i]);
     assert(charge[i]==charge_test[i]);
   }
   printf("\nTest #%d passed.\n", i_test);
   }


return 0;
}


int** atomtypes_test( int natoms, char **namespt , int *type, double *charge , int *num_types ){

int *temp,i,j,k,**ptr;
int *elem_ptr, **usr_ptr;
*num_types = 1;
temp = (int *) malloc( natoms*( sizeof(int) ) ) ;

/*
 * The large if else if structure below assigns a number to each atom type
 * present in the Gulp input, and also gives the charge associated with this
 * atom type.
 *
 * However, this structure will result in type lists like [1,7,1,7,1,7]
 * Lammps must have types between 1 and n_types, where n_types is the total
 * number of distinct types. I.e. the above array should be:  [1,2,1,2,1] etc.
 *
 * The integer values for each atom type are reassigned  to be in this form
 * and a description of which type is which is printed.
 */

for( i=0; i<natoms; i++ )
{
	if ( strcmp(*(namespt+i),"Ga1") == 0 ){
		charge[i] = 0.859;
		type[i] = 1;
	}
	else if ( strcmp(*(namespt+i), "Ga2") == 0 ){
		type[i] = 2;
		charge[i] = 0.859;
	}
	else if ( strcmp(*(namespt+i), "In1") == 0 ){
		charge[i] = 1.084;
		type[i] = 3;
	}
	else if ( strcmp(*(namespt+i), "In2") == 0 ){
		charge[i] = 1.084;
		type[i] = 4;
	}
	else if ( strcmp(*(namespt+i), "Al1") == 0 ){
		charge[i] = 1.500;
		type[i] = 5;
	}
	else if ( strcmp(*(namespt+i), "Al2") == 0 ){
		charge[i] = 1.5;
		type[i] = 6;
	}
	else if ( strcmp(*(namespt+i), "N101") == 0 ){
		charge[i] = -0.859;
		type[i] = 7;
	}
	else if ( strcmp(*(namespt+i), "N201") == 0 ){
		charge[i] = -0.859;
		type[i] = 8;
	}
	else if ( strcmp(*(namespt+i), "N102") == 0 ){
		charge[i] = -1.01925;
		type[i] = 9;
	}
	else if ( strcmp(*(namespt+i), "N202") == 0 ){
		charge[i] = -1.01925;
		type[i] = 10;
	}
	else if ( strcmp(*(namespt+i), "N103") == 0 ){
		charge[i] = -0.91525;
		type[i] = 11;
	}
	else if ( strcmp(*(namespt+i), "N203") == 0 ){
		charge[i] = -0.91525;
		type[i] = 12;
	}
	else if ( strcmp(*(namespt+i), "N104") == 0 ){
		charge[i] = -1.1795;
		type[i] = 13;
	}
	else if ( strcmp(*(namespt+i), "N204") == 0 ){
		charge[i] = -1.1795;
		type[i] = 14;
	}
	else if ( strcmp(*(namespt+i), "N105") == 0 ){
		charge[i] = -1.0755;
		type[i] = 15;
	}
	else if ( strcmp(*(namespt+i), "N205") == 0 ){
		charge[i] = -1.0755;
		type[i] = 16;
	}
	else if ( strcmp(*(namespt+i), "N106") == 0 ){
		charge[i] = -0.9715;
		type[i] = 17;
	}
	else if ( strcmp(*(namespt+i), "N206") == 0 ){
		charge[i] = -0.9715;
		type[i] = 18;
	}
	else if ( strcmp(*(namespt+i), "N107") == 0 ){
		charge[i] = -1.33975;
		type[i] = 19;
	}
	else if ( strcmp(*(namespt+i), "N207") == 0 ){
		charge[i] = -1.33975;
		type[i] = 20;
	}
	else if ( strcmp(*(namespt+i), "N108") == 0 ){
		charge[i] = -1.23575;
		type[i] = 21;
	}
	else if ( strcmp(*(namespt+i), "N208") == 0 ){
		type[i] = 22;
		charge[i] = -1.23575;
	}
	else if ( strcmp(*(namespt+i), "N109") == 0 ){
		charge[i] = -1.13175;
		type[i] = 23;
	}
	else if ( strcmp(*(namespt+i), "N209") == 0 ){
		charge[i] = -1.13175;
		type[i] = 24;
	}
	else if ( strcmp(*(namespt+i), "N110") == 0 ){
		charge[i] = -1.02775;
		type[i] = 25;
	}
	else if ( strcmp(*(namespt+i), "N210") == 0 ){
		charge[i] = -1.02775;
		type[i] = 26;
	}
	else if ( strcmp(*(namespt+i), "N111") == 0 ){
		charge[i] = -1.5;
		type[i] = 27;
	}
	else if ( strcmp(*(namespt+i), "N211") == 0 ){
		charge[i] = -1.5;
		type[i] = 28;
	}
	else if ( strcmp(*(namespt+i), "N112") == 0 ){
		charge[i] = -1.396;
		type[i] = 29;
	}
	else if ( strcmp(*(namespt+i), "N212") == 0 ){
		charge[i] = -1.396;
		type[i] = 30;
	}
	else if ( strcmp(*(namespt+i), "N113") == 0 ){
		charge[i] = -1.292;
		type[i] = 31;
	}
	else if ( strcmp(*(namespt+i), "N213") == 0 ){
		charge[i] = -1.292;
		type[i] = 32;
	}
	else if ( strcmp(*(namespt+i), "N114") == 0 ){
		charge[i] = -1.188;
		type[i] = 33;
	}
	else if ( strcmp(*(namespt+i), "N214") == 0 ){
		charge[i] = -1.188;
		type[i] = 34;
	}
	else if ( strcmp(*(namespt+i), "N115") == 0 ){
		charge[i] = -1.084;
		type[i] = 35;
	}
	else if ( strcmp( *(namespt+i), "N215") == 0 ){
		charge[i] = -1.084;
		type[i] = 36;
	}
}

printf("Types Assigned\n");


/* Now these type names must be reassigned so that they are in ascending
 * order and a difference of 1 between each type. To do this we will make a copy
 * of type[], sort this array, and then count the inversions in this sorted array
 * to determine the number of distinct types. We will then rename the types such that
 * the first type is type one, the third type is type three, etc.
 */

//Make a copy of the types and sort it.

for( i=0;i<natoms;i++ )
	temp[i] = type[i];

clock_t start = clock(), diff;
int msec;

printf("Sort\n");
qsort( temp, natoms,sizeof(int), cmpfunc );
diff = clock() - start;

msec = diff * 1000 / CLOCKS_PER_SEC;
printf("Time taken to sort: %d seconds %d milliseconds\n", msec/1000, msec%1000);


//Find number of distinct atom types and record in variable count.
int count = 0;
for( i=1; i< natoms; i++ )
	if( temp[i-1] != temp[i] )
		count++;
*num_types = count+1;
printf("There are %d distinct atom types\n", (*num_types ) );
/*This counts the number of inversions, so it will be "number of distinct types" -1
 *Since the first inversion gives "count =1" when there are infact two distinct types involved
 *every subsequent inversion will however add the correct additional distinct atom type.
 */

/*Declare variable types which contains one of each type, and also an
 *example of an atom number with this type, so that we can see what prop
 *-erties are associated with this type.
 *types[i][0] contains type index of each distinct type present.
 *types[i][1] contains the index of the first atom seen with the type in
 *column 0. So this index can be used with the names array as follows:
 *names( types[i][1] ) will give you the chemical name for the type contained
 *in types[i][0].
 */


//Declare and Initialise array "types".

ptr = ( int**) malloc( (*num_types) * sizeof ( int * ) ) ;
usr_ptr = ptr;
for( i=0; i< (*num_types); i++){
	*usr_ptr = ( int* ) malloc( 2*sizeof( int ) );
	elem_ptr = *usr_ptr;
	*(elem_ptr) = 0;		//**usr_ptr = 0;
	*(elem_ptr+1) = 0;		//*( *usr_ptr +1 )
	usr_ptr++;
	}
usr_ptr = ptr;

printf("Types on Heap allocated");

int types[ (*num_types ) ][2];										    //Use malloc here so scope extends beyond function.

for (i =0; i< (*num_types); i++ ){
	types[i][0] = 0;
	types[i][1] = 0;
   }

printf("\nTypes on stack allocated\n");
/*Since the array is ordered, each inversion will correspond to the presence
 *of a new distinct atom type. This is used to put each distinct atom type
 *into the array types, along with information about the type, as well as
 *determining the total number of types and printing this along with info
 *about each type.
 */


count = 0;
for( i=1; i< natoms; i++ )
	if( temp[i-1] != temp[i] ){
		*(*(ptr+count))   = temp[i-1];
		*(*(ptr+count+1)) = temp[i];
		types[count][0]   = temp[i-1];
		types[count+1][0] = temp[i];
		for( k=0; k < natoms; k++){										//This loop is necessary to get information about which atom Ids have this particular type.
			if( temp[i-1] == type[k] ){
				types[count][1] = k;
				*(*(ptr+count)+1) = k;
			}
			else if( temp[i] == type[k] ){
				types[count+1][1] = k;
				*(*(ptr+count+1)+1) = k;
			}
		}
		count++;
		}
//Note: this is an inefficient way to access the bits of memory malloc'd earlier
//you should just make a copy of the pointer to the addresses and then incremement
//that to get the successive memory locatations.



count++; //Now count = num_types

/*Now we need only reassign the types so that instead of being something
 *like 1,7,1,7,47,1,7 it's 1,2,1,2,3,1,2
 */

for( i=0; i < *num_types; i++ )
	**(ptr+i) = i+1;


for( i=0; i<natoms; i++ )
	for(j=0; j< count; j++)
		if( type[i] == types[j][0] )
			type[i]  = j+1 ;


return ptr;
}
//So ptr points to num_atoms arrays of 2 integers. ptr+i gives you types[i].
// *(ptr+i) will give you the address pointed to by ptr+i. This will be the
// address of the 2 element arrays. So *(*(ptr+i)) = types[i][0] . Can probably
// do (*(*(ptr+i)))[1] = types[i][1].
// *(*(ptr+i)+1) = types[i][1]
