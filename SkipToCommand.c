#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"funcs.h"

/* This function 'Skip_To_Command' sets the file pointer to point at the part of the file
 * one line past where the specified command is found. I.e. at the first piece
 * of data associated with that command for commands like frac, or cell etc.
 * The parameter cmdline is to contain any parameters that may be specified on the line
 * on which the command appears. For example, this will give a string containing all the paramters
 * for each of the potentials. The parameters can then be extracted using sscanf
 */ 

// So in brief, this returns in fptr a pointer to the data which occurs in the line beneath where the command
// occurs in the file. And it returns in cmdline the line on which the command occurs, which may contain needed parameters.

FILE *SkipTo( FILE *fptr , const char *command , char cmdline[] )
{
/* The char array cmdline[] contains the last string read in by the function, which should be the string which
 * starts with the letters of "command" */		
int i,size,yes;
size = (int) ( strlen(command) );

while(fgets ( cmdline,80, fptr ) != NULL )
	{
	yes =0;
	for (i=0; i < size ;i++)
		if ( cmdline[i] == command[i] )
			yes++;
	if ( yes == size )
		break;
	}
	
return fptr;	
}

FILE *safe_open( const char *filename , const char *mode)
{
	/*This function is more or less just fopen with some added safety features
	* the mode here is the same mode in fopen, where r is for read, w is for write,
	* rb is for read binary, wb is for write binary. The arguments work the exact same
	* as fopen. */
	
	FILE *fptr;
	
	if( ( fptr = fopen ( filename, mode) ) == NULL )        
		{
		puts ( "Cannot open input file " );
		puts ( filename ) ;
		exit( 1 );
		}
return fptr;
}


double Distance ( double x1, double y1, double z1, double x2, double y2, double z2 ) //Function to calculate distances between two different atoms
{	
return sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );
}
