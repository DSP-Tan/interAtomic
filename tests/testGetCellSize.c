#include <stdio.h>
#include<assert.h>
#include "../funcs.h"

int getCellSize( char datafile[] , double *cell );

// To Do: Test 3 sample input cells, with the values to test and the filenames
// in arrays so that you can just look over them.

// Consider in the make file actually creating these cells if they do not exist.
// They should not really be in the repo as they would be a bit large.

// Here we will test getCellSize using 3 sample input cells of a known size.
int main(){

double cell[3];
int test_a_1, test_b_1, test_c_1, test_natoms_1;
int test_a_2, test_b_2, test_c_2, test_natoms_2;
int test_a_3, test_b_3, test_c_3, test_natoms_3;

char datafile1[80]= "10_10_20_1.txt";
char datafile2[80]= "2_2_4_1.txt";
char datafile3[80]= "3_7_12_1.txt";

int natoms = getCellSize(datafile1,cell);
test_natoms_1 = natoms==4000;
printf("natoms: %d\n",test_natoms_1);

double e = 1e-8;
test_a_1 = 31.889999389648438 -e<= cell[0] && cell[0]<= 31.889999389648438 +e;
test_b_1 = 27.616739273071289 -e<= cell[1] && cell[1]<= 27.616739273071289 +e;
test_c_1 = 51.849998474121094 -e<= cell[2] && cell[2]<= 51.849998474121094 +e;
printf("cell: %d %d %d\n",test_a_1 ,test_b_1,test_c_1);

assert(test_natoms_1);
assert(test_a_1);
assert(test_b_1);
assert(test_c_1);

return test_natoms_1 && test_a_1 && test_b_1 && test_c_1;

}
