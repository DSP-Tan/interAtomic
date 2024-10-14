#include<stdio.h>
#include<string.h>
#include"funcs.h"
#include<math.h>
#include<malloc.h>
#include<stdlib.h>
#include<time.h>

int compare(const void *p1, const void *p2);

struct bond* makebonds(const char *source, int natoms, double *x_, double *y_, double *z_,int **types, int *type, int num_types,
char **namespt , double cutoff, double *cell , int *NumberOfBonds,struct bond *bondcoeffs,
int *NumberOfBondTypes )
{
double kGaN, kInN, kAlN;
//These are the force constants for the harmonic potential between Gan, InN and AlN respectively

double r0_GaN, r0_InN, r0_AlN,r0;

double xi_unwrapped, xj_unwrapped, yi_unwrapped, yj_unwrapped, zi_unwrapped, zj_unwrapped;
double r_ij,r_ij_unwrapped;
int Bond_ID;
FILE *ptr;
char string[80],  name1[10],  name2[10], name3[10];
double trash,K;
int NumBondTypes, NumDistinctBonds, NumDistinctPresentBonds,count,i,j,k,l;
//These variables all refer to bond types.

/*--------------------------------------------------------------------*/
/*------------Initialise Variables--------------------------------------
 * -------------------------------------------------------------------*/
kGaN = 9.0650;
kInN = 7.005;  //Note these are the GULP values/2 because lammps incorporates the hallf into the constant.
kAlN = 10.5;
r0_GaN = 1.949;
r0_InN = 2.155;
r0_AlN = 1.895;

count = 0;
for(i=0;i<10;i++)
	name1[i] = name2[i] = name3[i] = '\0';
for(i=0;i<80;i++)
	string[i] = '\0';


/*--------------------------------------------------------------------*/
/*----Read force constants from potentials file-----------------------*/
/*--------------------------------------------------------------------*/

ptr = safe_open( source , "r" );
SkipTo( ptr, "harm", string );

while(fscanf(ptr,"%s %s %lf %lf %s %d %d", name1, name2, &K, &r0, string,&j,&j )==7){
   printf("%s %s %lf %lf\n",name1,name2,K,r0);
   if(  (name1[0] == 'G' &&  name2[0] == 'N') || (name1[0] == 'N' &&  name2[0] == 'G') ){
      kGaN = K/2.0;
	   r0_GaN = r0;
	   }
   else if(  (name1[0] == 'I' &&  name2[0] == 'N') || (name1[0] == 'N' &&  name2[0] == 'I') ){
      kInN = K/2.0;
      r0_InN = r0;
	   }
   else if(  (name1[0] == 'A' &&  name2[0] == 'N') || (name1[0] == 'N' &&  name2[0] == 'A') ){
      kAlN = K/2.0;
	   r0_AlN = r0;
	   }
   }
rewind(ptr);


SkipTo( ptr, "polynomial", string );
fgets( string, 80 , ptr );				       //Skip line containing '1'


/*To get the maximum number of bond types, we go through the polynomial
 * section of the potentials file, as this uniquely determines each bond.
 * we check if the first atom type of a given line is present in our previously
 * determined atom types, and then we check if the second atom type is
 * also present.
*This gives the maximum possible number of bond types since it says if there
* is a Ga and a N present, there might be a bond between them. Though you could
* imagine a scenario where there are a Ga and an N but they are too far to bond. */

NumBondTypes = 0;
while( fscanf( ptr, "%s %s %lf %lf %lf %lf", name1, name2, &trash, &K, &trash, &trash ) == 6 ){
   count=0;
	for( i=0; i< num_types; i++ ){
		if( strcmp( name1, *(namespt + *(*(types+i)+1) )  ) == 0 ){ 			//Conditions like this can always be replaced by if ( ! condition ). If not condition. If not strcmp.
			count++;
         for( j=0; j< num_types; j++ )
            if( strcmp( name2, *(namespt + *(*(types+j)+1) )  ) == 0 ){
               count++;
               break;
               }
         }
		if( count == 2 ){
			NumBondTypes++;
			break;
		   }
	   }
   }

printf("\nThere are at most %d Bond Types\n",NumBondTypes);

//TO DO: Now that we have pared down number of bonds, fill each of these
//bonds with information about bond the two force constants, the atom types.

struct bond Bonds[NumBondTypes];								//Note: The bond ID is also stored in the element number here.
for( i=0; i< NumBondTypes; i++ )
	Bonds[i].present = 0;

rewind( ptr );
count = 0;

SkipTo( ptr, "polynomial", string );
fgets( string, 80 , ptr );				//Skip line containing '1'
//-----------NOTE!!!!!!!!:This method requires you to remove the #'s that are between the polynomial coeffs------------------


while( fscanf( ptr, "%s %s %lf %lf %lf %lf", name1, name2, &trash, &K, &trash, &trash ) == 6 )
	for( i=0; i< num_types; i++ )
		if( strcmp( name1, namespt[ *(*(types+i)+1) ] ) == 0  )
			for( j=0; j< num_types; j++ )
			   if( strcmp( name2, namespt[ *(*(types+j)+1) ] ) == 0  ){
			      Bonds[count].atType1 = **(types+i);
			      Bonds[count].atType2 = **(types+j);
			      Bonds[count].K1 = K;
			      Bonds[count].BType = count+1;
			      if(  namespt[ *(*(types+i)+1) ][0] == 'G' ||  namespt[ *(*(types+j)+1) ][0] == 'G' ){
				     Bonds[count].K2 = kGaN;
				     Bonds[count].r0 = r0_GaN;
				     }
				   else if(  namespt[ *(*(types+i)+1) ][0] == 'A' ||  namespt[ *(*(types+j)+1) ][0] == 'A' ){
				     Bonds[count].K2 = kAlN;
				     Bonds[count].r0 = r0_AlN;
				     }
				   else if(  namespt[ *(*(types+i)+1) ][0] == 'I' ||  namespt[ *(*(types+j)+1) ][0] == 'I' ){
				     Bonds[count].K2 = kInN;
				     Bonds[count].r0 = r0_InN;
				     }
				  count++;
				  }

fclose( ptr );

printf("count=%d, NumbondTypes = %d\n",count,NumBondTypes);
//for( i=0;i<count; i++ )
//	printf( "Atom1 %d Atom2 %d K1 %lf K2 %lf r0 %lf BondType %d\n" , Bonds[i].atType1, Bonds[i].atType2, Bonds[i].K1 , Bonds[i].K2 ,Bonds[i].r0, Bonds[i].BType );

//TO DO now you make sure that no two bond types give the same information
//If they do, set the bond type to be the same.

for( i=0;i<NumBondTypes;i++)
	for( j=0; j<NumBondTypes; j++ )
		if( Bonds[i].K1 == Bonds[j].K1 && Bonds[i].K2 == Bonds[j].K2 ){
			Bonds[i].BType = ( Bonds[i].BType > Bonds[j].BType ? Bonds[j].BType : Bonds[i].BType );
			Bonds[j].BType = ( Bonds[j].BType > Bonds[i].BType ? Bonds[i].BType : Bonds[j].BType );
			//Conditional operator assigns type to be printed. It goes Condition ? returnvalue1 : returnvalue 2 ; Turn this into a big case thing.
		   }

clock_t start = clock(), diff;
int msec;
qsort( Bonds, NumBondTypes, sizeof(struct bond), compare );
diff = clock() - start;

//msec = diff * 1000 / CLOCKS_PER_SEC;
//printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
//printf("\n\nThe sorted list of bond types is:\n");
//printf( "Type i\t Type j\t K_1\t\t K_2\t\t r_0\t\t Bond Type\n" );
//for( i=0;i<count; i++ )
//	printf( "%d\t %d\t %lf\t %lf\t %lf\t %d\n" , Bonds[i].atType1, Bonds[i].atType2, Bonds[i].K1 , Bonds[i].K2 ,Bonds[i].r0, Bonds[i].BType );



//--------------Write Bonding Information-------------------------------------------//

struct bond *bondptr;
struct bond *usrptr;


bondptr = ( struct bond * ) malloc( (natoms*2)  * sizeof (struct bond ) ) ;
usrptr = bondptr;
*NumberOfBonds = natoms*2;
Bond_ID = 1;

for (i=0; i<natoms; i++)
	for (j= i+1; j<natoms; j++)
		for( k=0; k< NumBondTypes; k++)
		   if( (type[i] == Bonds[k].atType1 && type[j] == Bonds[k].atType2 ) || (type[j] == Bonds[k].atType1 && type[i] == Bonds[k].atType2) )
		      {

		      r_ij = Distance( x_[i], y_[i], z_[i] , x_[j], y_[j], z_[j] );

		      if( r_ij < cutoff ){
               Bonds[k].present = 1;
               *usrptr = Bonds[k];
               usrptr->atID1 = i+1;
               usrptr->atID2 = j+1;
               usrptr->B_ID = Bond_ID;
               usrptr++;
               Bond_ID++;
               }

			  if (  (x_[i] > cell[0] - cutoff && x_[j] < cutoff) || (x_[j] > cell[0] -cutoff && x_[i] < cutoff)
              || (y_[i] > cell[1] - cutoff && y_[j] < cutoff) || (z_[i] > cell[2] - cutoff && z_[j] < cutoff)
              || (y_[j] > cell[1] - cutoff && y_[i] < cutoff) || (z_[j] > cell[2] - cutoff && z_[i] < cutoff) )
				{

			    xi_unwrapped = x_[i];
			    yi_unwrapped = y_[i];
			    zi_unwrapped = z_[i];

			    xj_unwrapped = x_[j];
			    yj_unwrapped = y_[j];
			    zj_unwrapped = z_[j];


				/*If atom i is near end of cell in any dimension, and atom j is near beginning of cell
				 * in this same dimension. Calculate the distance and make the bond between the image atoms
				 */

				if (x_[i] > cell[0] -cutoff && x_[j] < cutoff)
					xi_unwrapped = x_[i] - cell[0];
				if (y_[i] > cell[1] - cutoff && y_[j] < cutoff)
					yi_unwrapped = y_[i] - cell[1];
				if (z_[i] > cell[2] - cutoff && z_[j] < cutoff)
					zi_unwrapped = z_[i] - cell[2];


				/*If atom j is near end of cell in any dimension, and atom i is near beginning of cell
				 * in this same dimension. Calculate the distance and make the bond between the image atoms
				 */

				if (x_[j] > cell[0] -cutoff && x_[i] < cutoff)
					xj_unwrapped = x_[j] - cell[0];
				if (y_[j] > cell[1] - cutoff && y_[i] < cutoff)
					yj_unwrapped = y_[j] - cell[1];
				if (z_[j] > cell[2] - cutoff && z_[i] < cutoff)
					zj_unwrapped = z_[j] - cell[2];


				r_ij_unwrapped = Distance( xi_unwrapped, yi_unwrapped, zi_unwrapped, xj_unwrapped, yj_unwrapped, zj_unwrapped );


				if ( r_ij_unwrapped < cutoff )
					{
               Bonds[k].present = 1;
               *usrptr = Bonds[k];
               usrptr->atID1 = i+1;
               usrptr->atID2 = j+1;
               usrptr->B_ID = Bond_ID;
				   usrptr++;
					Bond_ID++;
					}

				if( r_ij_unwrapped < cutoff && r_ij < cutoff ){
				   printf("\n\nWARNING double bond created\n\n");
				   exit(-1);
               }
				}
			}

Bond_ID--;			//Bond_ID was incremented in last iteration of loop.


printf("There are %d pairs of atoms bonded to each other\n", Bond_ID );
printf("------------------------------------------------------------------------------------------------------------\n");



/*--------------------------------------Bond coefficients-------------------------------------------*/
// - Find distinct bond types
// - Make bond array of these
// - Use array to print Bond Coeffs part of read_data file.

//Find number of distinct bond types and record in variable count.
NumDistinctBonds = 0;
for( i=1; i< NumBondTypes; i++ )
	if( Bonds[i-1].BType != Bonds[i].BType )
		NumDistinctBonds++;
NumDistinctBonds++;
printf("There are %d distinct Bonds possible\n", NumDistinctBonds);

//Know there are NumDistinctBonds distinct bonds, make array and put appropriate information for bonds in it.
struct bond distinct_Bonds[NumDistinctBonds];
for( i=1,k=1; i< NumBondTypes; i++ )
   if(NumDistinctBonds == 1){                   //This method of bond assignment does not work if there is only one bond type.
      if (Bonds[i].present==1){
         distinct_Bonds[0] = Bonds[i];
         break;
         }
      else
         continue;
      }
	else if( Bonds[i-1].BType != Bonds[i].BType ){
		distinct_Bonds[k-1] = Bonds[i-1];
		distinct_Bonds[k] = Bonds[i];
		k++;
      }

printf("\n\nThe sorted list of distinct bond types is:\n");
printf( "Type i\t Type j\t K_1\t\t K_2\t\t r_0\t\t Bond Type\n" );
for( i=0;i<NumDistinctBonds; i++ )
	printf( "%d\t %d\t %lf\t %lf\t %lf\t %d\n" , distinct_Bonds[i].atType1,
   distinct_Bonds[i].atType2, distinct_Bonds[i].K1 , distinct_Bonds[i].K2
   ,distinct_Bonds[i].r0, distinct_Bonds[i].BType );


/*The for loop below fixes the fact that the distinct_Bonds array says that certain bonds are not present which
 *are present. This occurs when two bonds in Bonds which are the same but between different atom types have one of
 *these atom type pairs but not the other. So even though this bond is present, it takes it's information from the
 *element of Bonds[] corresponding to the atom pairs which aren't present.
*/

for( i=0; i<NumBondTypes; i++ )
	for (j =0; j< NumDistinctBonds; j++)
		if ( Bonds[i].K1 == distinct_Bonds[j].K1 && Bonds[i].K2 == distinct_Bonds[j].K2 && Bonds[i].present ==1 )
			distinct_Bonds[j].present =1;


printf("The number of total possible bond types is:  %d\n", NumBondTypes);
printf("The number of unique possible bond types is: %d\n",NumDistinctBonds);


/*So this array distinct_Bonds contains only one of each bond type, but does not contain the full
 *bonding information about which atom type bonds to which. You can cross reference the bond types of this
 *array with the full information of the other array if you need to.
 *
 *Before we print out the Bond_Coeffs, we will use the information of how many distinct bonds there are to write properly
 *and with every possible check the bonding information for each atom:
*/

NumDistinctPresentBonds=0;
for(j=0,count=0; j< NumDistinctBonds; j++)
	NumDistinctPresentBonds += distinct_Bonds[j].present;
printf("Out of %d possible distinct bonds, there are %d present\n", NumDistinctBonds, NumDistinctPresentBonds );
printf("Make sure the number of distinct present bond types is the same as that stated at the top of the data file\n");
*NumberOfBondTypes = NumDistinctPresentBonds;

struct bond Present_Bonds[NumDistinctPresentBonds];
i=0;
for(j=0;j<NumDistinctBonds;j++)
   if( distinct_Bonds[j].present == 1 ){
      Present_Bonds[i] = distinct_Bonds[j];
      Present_Bonds[i].BType = i+1;
      i++;
      }

/*Now make sure that the array of bonds with the full information about the bonding constants for particular pairs
 *of atom types has the right Bond Type associated with each of these bonds
 */
for( k=0; k<NumBondTypes; k++ )
   for( l=0; l< NumDistinctPresentBonds; l++)
      if ( Bonds[k].K1 == Present_Bonds[l].K1 && Bonds[k].K2 == Present_Bonds[l].K2 )
         Bonds[k].BType = Present_Bonds[l].BType;

usrptr = bondptr;
for(i=0;i<*NumberOfBonds;i++)
   for(j=0;j<NumDistinctPresentBonds;j++)
      if(usrptr->K1 == Present_Bonds[j].K1
      && usrptr->K2 == Present_Bonds[j].K2){
         usrptr->BType = Present_Bonds[j].BType;
         usrptr++;
         }




//------------------Write Bond Coeffs---------------------------------//

for(i=0;i<NumDistinctPresentBonds;i++)
   bondcoeffs[i] = Present_Bonds[i];

//Return array of structs containing all bonds:
return bondptr;


}



int compare(const void *p1, const void *p2)
{
   const struct bond *elem1 = p1;
   const struct bond *elem2 = p2;
   return (int)(elem1->BType - elem2->BType);
}

/*Note: When making these bonds, if two atoms which should be bonded are beyoned
 *the cutoff are distance from each other, then they will not be bonded. This means
 *there will never be any bonded interactions between the two atoms no matter how
 *close they come after initialisation of the minimisation. In GULP, there could be
 * two atoms between which there could be a bond, but which are further than cutoff
 * from each other. If these atoms come close during the simulation, they will
 * experience forces between each ohter.
 *
 * It is therefore imporant to ensure that your Gulp cutoffs and initial atomic positions
 * correspond properly to the bonding topology you would like to have for the whole run.
 */
