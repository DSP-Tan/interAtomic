#include<stdio.h>
#include"funcs.h"
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

int anglesort( int n, struct Angles *ptr );
int compaare(const void *p1, const void *p2);

struct Angles* makeangles(const char *source,int natoms, double *x_,double *y_,
double *z_, char **namespt, double cutoff,double *cell,
struct bond *bonds, int NumBonds, int *NumberOfAngles,
struct Angles *PosibAngs, int *NumberOfAngleTypes)
{
int i,j,k,l; //loop indices
char string[80],  name1[10],  name2[10], name3[10];
double trash;

//-------------------------------------------Angular Information-------------------------------//
/*NOTE doing the angular information is specific to how the .gin files are given to you. There is nothing general
 * about this. This assumes in different ways that the angular interactions are given in exactly the form they are given
 * in miguel's gulp scripts used to relax AlInGaN.
 *
 * It could be made more general to handle variations in the types of angular potentials, but there's probably no need.
 *
 */



//First we must find the number of distinct angle types present:

//To find the present angles, we say that any two bonds which have an atom in common form an angle.



/*First all the possible types of angle interaction and the parameters
 *for these types are recorded from the GULP input script. There are 9 of
 *these in the InAlGaN input scripts. This program is therefore specific
 *to that particular script, and cannot be used on others.
 */

FILE *ptr;

ptr = safe_open( source , "r" );
ptr = SkipTo( ptr, "three" , string );
for( i=0; i< 9; i++ ){
fscanf(ptr,"%s %s %s %lf %lf %lf %lf %s %d %d",name1,name2,name3, &PosibAngs[i].K2, &PosibAngs[i].Theta0,&trash,&trash,string,&l,&l);
PosibAngs[i].atom1 = name1[0];
PosibAngs[i].atom2 = name2[0];
PosibAngs[i].atom3 = name3[0];
}

ptr = SkipTo( ptr, "bacross" , string );
for( i=0; i< 9; i++ )
   fscanf( ptr, "%s %s %s %lf %lf %lf %lf %lf %lf %lf %s %d %d %d %d", name1,name2,name3, &PosibAngs[i].N1, &PosibAngs[i].N2,&PosibAngs[i].r1,&PosibAngs[i].r2,&trash, &trash,&trash, name1,&l,&l,&l,&l);

ptr = SkipTo( ptr, "bcross" , string );
for( i=0; i< 9; i++ )
   fscanf( ptr, "%s %s %s %lf %lf %lf %lf %lf %s %d %d %d", name1,name2,name3, &PosibAngs[i].M, &trash,&trash,&trash,&trash, name1,&l,&l,&l );


fclose( ptr );

int NumAngles;							      // This is the total number of angular interactions.
int Ang1, Ang2, Ang3;				      // These variables contain, for a given angle, the ID's of the three atoms participating in this angle
int NumAngleTypes;						   // This is the number of types of angle.
int AngleAtom1, AngleAtom2,AngleAtom3;


NumAngles = natoms*6;                  // The number of angles in wurzite GaN is 6 times the number of atoms.
*NumberOfAngles = NumAngles;
struct Angles *angleptr,*usrptr;

angleptr = ( struct Angles * ) malloc(NumAngles * sizeof (struct Angles ) ) ;
usrptr = angleptr;


NumAngles=0;
printf("-------------------------------------------------\n\n");


for(i=0;i<NumBonds;i++)
   for(j=i+1;j< NumBonds; j++ )
      if(bonds[i].atID1 == bonds[j].atID1 || bonds[i].atID1 == bonds[j].atID2 || bonds[i].atID2 == bonds[j].atID1 || bonds[i].atID2 == bonds[j].atID2 ){
         if(bonds[i].atID1 == bonds[j].atID1 || bonds[i].atID1 == bonds[j].atID2){
			   Ang2 = bonds[i].atID1-1;
			   Ang1 = ( bonds[i].atID1 == bonds[j].atID1 ? bonds[j].atID2 : bonds[j].atID1 ) -1 ;
			   Ang3 = bonds[i].atID2 -1;
            }
         else if(bonds[i].atID2 == bonds[j].atID1 || bonds[i].atID2 == bonds[j].atID2){
			   Ang2 = bonds[i].atID2 -1;
			   Ang1 = ( bonds[i].atID2 == bonds[j].atID1 ? bonds[j].atID2 : bonds[j].atID1 ) -1 ;
			   Ang3 = bonds[i].atID1 -1;
            }

		 //printf("Ang1 = %d Ang2 = %d Ang3 = %d\n\n",Ang1,Ang2,Ang3);
		 /*The purpose of the -1's is so that we can use these atom ID's to index into the array of atom names
		  *which starts at 0, and use this to find the species information about the atoms with these IDs
		  */

		 for(k=0;k<9;k++)
		    if(PosibAngs[k].atom1 == **(namespt+Ang2))
		       if( (PosibAngs[k].atom2 == **(namespt+Ang1) && PosibAngs[k].atom3 == **(namespt+Ang3))
              || (PosibAngs[k].atom2 == **(namespt+Ang3) && PosibAngs[k].atom3 == **(namespt+Ang1)) )
		          {
                AngleAtom1 = (PosibAngs[k].atom2== **(namespt+Ang1) ? Ang1 :Ang3);
                AngleAtom2 = Ang2;
                AngleAtom3 = (PosibAngs[k].atom3== **(namespt+Ang3) ? Ang3 :Ang1);
                PosibAngs[k].present = 1;
                *usrptr = PosibAngs[k];
                usrptr->at1 = AngleAtom1+1;
                usrptr->at2 = AngleAtom2+1;
                usrptr->at3 = AngleAtom3+1;
                usrptr++;
                NumAngles++;
		          }
      }


l=1;
for(i=0;i<9;i++)
	if ( PosibAngs[i].present == 1 ){
		PosibAngs[i].type = l;
		if( PosibAngs[i+1].present==1
      && PosibAngs[i+1].K2 == PosibAngs[i].K2
      && PosibAngs[i+1].M == PosibAngs[i].M
      && PosibAngs[i+1].N1 == PosibAngs[i].N1
      && PosibAngs[i+1].N2 == PosibAngs[i].N2
      && PosibAngs[i+1].r1 == PosibAngs[i].r1
      && PosibAngs[i+1].r2 == PosibAngs[i].r2 )
			continue;
		l++;
	   }

printf("Let's look at the angle types:\n\n");
for(i=0;i<9;i++)
   printf("Type: %d,Present:%d, K2: %lf, M: %lf, N1: %lf, N2: %lf,r1: %lf, r2: %lf, %c,%c,%c\n", PosibAngs[i].type,PosibAngs[i].present,PosibAngs[i].K2,PosibAngs[i].M,PosibAngs[i].N1,PosibAngs[i].N2,PosibAngs[i].r1,PosibAngs[i].r2,PosibAngs[i].atom1,PosibAngs[i].atom2,PosibAngs[i].atom3);

printf("\nThe number of angles is : %d\nnatomsx6 = %d\n",NumAngles,natoms*6);


//Now you must set the angle types of the angles in the supercell.
usrptr = angleptr;
for(i=0;i<NumAngles;i++)
   for(j=0;j<9;j++)
      if( usrptr->K2 == PosibAngs[j].K2 ){    //This one parameter is enough to specify an angle
         usrptr->type = PosibAngs[j].type;
         usrptr++;
         break;
         }

NumAngleTypes = l-1;
*NumberOfAngleTypes = NumAngleTypes;
printf("There are %d Angle Types\n", NumAngleTypes);


/*Before we write these coefficients, we may want to eliminate repeated
 *Angle types from our array of PosibAngs[k]
 */
int repeat = 0;
for(i=1;i<9;i++)
	if( PosibAngs[i].present == 1 && PosibAngs[i-1].present ==1)
	   if(PosibAngs[i].type == PosibAngs[i-1].type)
	      repeat++;

printf("\nThere are %d repeated angle types\n\n",repeat);
qsort( angleptr, *NumberOfAngles,  sizeof(struct Angles), compaare );
printf("\nSort again?\n");
anglesort(*NumberOfAngles,angleptr);
printf("Sorted\n");

return angleptr;

}



int anglesort( int n, struct Angles *ptr ){
   int status, i;
   struct Angles a;
   status = 0;
   for( i=0; i< (n-1);i++ ){
      if ( ptr[i].type > ptr[i+1].type ){
         status++;
         a = ptr[i];
         ptr[i] = ptr[i+1];
         ptr[i+1] = a;
	     }
    }

	if( status == 0 )
	   return status;
	else
	   anglesort( n, ptr );

	return status;
}

int compaare(const void *p1, const void *p2)
{
   const struct Angles *elem1 = p1;
   const struct Angles *elem2 = p2;
   return (int)(elem1->type - elem2->type);
}
