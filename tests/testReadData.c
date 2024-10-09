#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "../funcs.h"

int main(){

char datafile[] = "2_2_4_1.txt";
int natoms = 32;
double testX[natoms], testY[natoms], testZ[natoms];
double x[natoms], y[natoms], z[natoms];
char** tName;
char** names;

int i;

tName=(  char **)malloc( natoms*sizeof( char*) );
names=(  char **)malloc( natoms*sizeof( char*) );
for(i=0;i<natoms;i++){
   tName[i] = ( char* )malloc( 8*sizeof( char ) );
   names[i] = ( char* )malloc( 8*sizeof( char ) );
   memset(tName[i], '\0', 8);
   memset(names[i], '\0', 8);
   }


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

tName[0] ="Ga1";   testX[0] = 0.000000;testY[0] = 0.000000; testZ[0] = 0.000000;
tName[1] ="N101";  testX[1] = 0.250000;testY[1] = 0.166667; testZ[1] = 0.058985;
tName[2] ="Ga1";   testX[2] = 0.500000;testY[2] = 0.000000; testZ[2] = 0.000000;
tName[3] ="N101";  testX[3] = 0.750000;testY[3] = 0.166667; testZ[3] = 0.058985;
tName[4] ="Ga1";   testX[4] = 0.250000;testY[4] = 0.500000; testZ[4] = 0.000000;
tName[5] ="N101";  testX[5] = 0.000000;testY[5] = 0.666667; testZ[5] = 0.058985;
tName[6] ="Ga1";   testX[6] = 0.750000;testY[6] = 0.500000; testZ[6] = 0.000000;
tName[7] ="N101";  testX[7] = 0.500000;testY[7] = 0.666667; testZ[7] = 0.058985;
tName[8] ="Ga2";   testX[8] = 0.250000;testY[8] = 0.166667; testZ[8] = 0.250000;
tName[9] ="N201";  testX[9] = 0.000000;testY[9] = 0.000000; testZ[9] = 0.308985;
tName[10] ="Ga2";  testX[10]= 0.750000;testY[10]= 0.166667; testZ[10]= 0.250000;
tName[11] ="N201"; testX[11]= 0.500000;testY[11]= 0.000000; testZ[11]= 0.308985;
tName[12] ="Ga2";  testX[12]= 0.000000;testY[12]= 0.666667; testZ[12]= 0.250000;
tName[13] ="N201"; testX[13]= 0.250000;testY[13]= 0.500000; testZ[13]= 0.308985;
tName[14] ="Ga2";  testX[14]= 0.500000;testY[14]= 0.666667; testZ[14]= 0.250000;
tName[15] ="N201"; testX[15]= 0.750000;testY[15]= 0.500000; testZ[15]= 0.308985;
tName[16] ="Ga1";  testX[16]= 0.000000;testY[16]= 0.000000; testZ[16]= 0.500000;
tName[17] ="N101"; testX[17]= 0.250000;testY[17]= 0.166667; testZ[17]= 0.558985;
tName[18] ="Ga1";  testX[18]= 0.500000;testY[18]= 0.000000; testZ[18]= 0.500000;
tName[19] ="N101"; testX[19]= 0.750000;testY[19]= 0.166667; testZ[19]= 0.558985;
tName[20] ="Ga1";  testX[20]= 0.250000;testY[20]= 0.500000; testZ[20]= 0.500000;
tName[21] ="N101"; testX[21]= 0.000000;testY[21]= 0.666667; testZ[21]= 0.558985;
tName[22] ="Ga1";  testX[22]= 0.750000;testY[22]= 0.500000; testZ[22]= 0.500000;
tName[23] ="N101"; testX[23]= 0.500000;testY[23]= 0.666667; testZ[23]= 0.558985;
tName[24] ="Ga2";  testX[24]= 0.250000;testY[24]= 0.166667; testZ[24]= 0.750000;
tName[25] ="N201"; testX[25]= 0.000000;testY[25]= 0.000000; testZ[25]= 0.808985;
tName[26] ="Ga2";  testX[26]= 0.750000;testY[26]= 0.166667; testZ[26]= 0.750000;
tName[27] ="N201"; testX[27]= 0.500000;testY[27]= 0.000000; testZ[27]= 0.808985;
tName[28] ="Ga2";  testX[28]= 0.000000;testY[28]= 0.666667; testZ[28]= 0.750000;
tName[29] ="N201"; testX[29]= 0.250000;testY[29]= 0.500000; testZ[29]= 0.808985;
tName[30] ="Ga2";  testX[30]= 0.500000;testY[30]= 0.666667; testZ[30]= 0.750000;
tName[31] ="N201"; testX[31]= 0.750000;testY[31]= 0.500000; testZ[31]= 0.808985;

read_data ( datafile, x, y , z, natoms, names);

double e = 1e-8;
for(i=0;i<natoms;i++){
  assert(x[i] >= testX[i]-e && x[i] <= testX[i]+e);
  assert(y[i] >= testY[i]-e && y[i] <= testY[i]+e);
  assert(z[i] >= testZ[i]-e && z[i] <= testZ[i]+e);
  assert( strcmp(tName[i],names[i])==0 );
  }

return 0;

}
