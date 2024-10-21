We will test all functions on 3 GaN, 3In, and 3 real InGaN supercells.

For the larger cells we will just have to test the input and output using frozen older 
hard coded version of the functions we are testing, which we know work. 

This way we can do whatever changes we want and be sure that they still give the correct
behaviour.

In fact you can do it for all of the tests becaues the other way is a pain in the Hole.

We will start with supercell files which are contained in the ../supercell folder. However,
as time goes on we will want to do tests which do not rely on the presence of such files, 
but which rather generate them using a conerted to C version of the fortran files in the
supercell folder.


Now: find out how to make a mixed supercell of a more reasonable size.

We will do 3 AlInGaN cells with 10% In and 10% Al of such a max size that does
not take too long to make. 

Or else you should have a pure GaN, a pure InN, a pure AlN, an InGaN, an AlGaN, and an AlInGaN

The good news is that to completely re-create this you only need to reproduce supercell.f in C.
That should not be too hard to do.
