******** Introduction ******************


Self-modifying CGP is a form of CGP which includes self-modifying operation. This means that
evolved programs include instructions that change the program. Thus they can be iterated (indefinitely) and
a possibly infinite sequence of programs can be evolved!

For further details consult the tutorials on the CGP web site or buy the CGP book.

The smcgp programs have implemented a variant of the standard form described in publications. I have renamed
some primitive functions and also I have implemented output functions difefrently. In my version
the output functions write to registers which are read for fitness assessment. This is a cleaner way to do things
that the output methods as published. 


The only description of the new methods is published in a small workshop in Paris in October 2011. This is the only account
so far of my implementation of SMCGP. So if you end up publishing anything using my smcgp program I think you should cite

  
@INPROCEEDINGS {miller2011devleann,
    title = {Neuro-centric and holocentric approaches to the evolution of developmental neural networks},
    author = {J.F.~Miller},
    editor = {T.~Kowaliw and N.~Bredeche and R.~Doursat},
    title = {Proceedings of DevLeaNN: A Workshop on Development and Learning in Artificial Neural Networks},
    address = {Paris, France},
    pages = {18--25},
    year = {2011},
    url = {http://devleann.iscpif.fr}
}

The paper is downloadable from the CGP web site (see my publications by going to http://www.cartesiangp.co.uk/jfm-publications.html)


******** SMCGP program *********************

The self-modifying Cartesian GP program consists of the following files:

smcgp.c
smcgp-functions.c
node-function.c
smcgp.h

smcgp.c contains the main function.
smcgp-functions.c contains almost all of the function implementations
node-function.c contains the implementation of the primitive function set
smcgp.h contains all defined constants, global variables and function prototypes


 The program allows two kinds of problems to be tackled. Boolean, where
all the primitive functions relate to binary data and real-valued, where the primitive functions are mathematical functions opearting
on real-valued data.

If you want to work on the program for boolean data, you need to make sure one of the defined symbols is uncommented in smcgp.h.


/* */
#define DATA_IS_UNSIGNED_INT
/* */


/* comment this in if you want to work with integer data */
#define DATA_IS_INT
*/

/* comment this in if you want to work with real-value data */
/*
#define DATA_IS_DOUBLE
*/

So here DATA_IS_UNSIGNED_INT is uncommented. Under these circumstances, when the code is compiled it will generate an executable
program that will read in specifications of boolean functions in compressed 32-bit format (see later). However, if instead we
have 

/* 
#define DATA_IS_UNSIGNED_INT
*/


/* comment this in if you want to work with integer data */
#define DATA_IS_INT
*/

/* comment this in if you want to work with real-value data */
/* */
#define DATA_IS_DOUBLE
/* */

now DATA_IS_DOUBLE is uncommented so when the program is compiled it will generate an executable that reads floating point data.

Currently I have not implemented anything particular to integer data (so for now keep DATA_IS_INT commented out)

The different data types come with their own .par files. The parameter files (.par) define all the parameters required for
the program (i.e. evolutionary parameters, control of runs, what files the program generates etc).



I provide two .par files: smcgp-boolean.par and smcgp-double.par. 


Most of the parameters in the .par files are standard to CGP (note standard CGP can be downloaded from the cgp web site)


The SMCGP code has defined constants that define datatypes that will allow it to work on three types of data

compressed boolean (as in the boolean experiments). All truth tables are compressed into 32 bit words and the
boolean functions operate on these words. So the truth tables epar3.plu, epar4.plu etc have been compressed
from truth tables in binary.


***************** RUNNING THE SMCGP PROGRAMS *******************

To run smcgp you need to say

smcgp smcgp.par 

where smcgp.par is file with a valid list of parameters.


***************** EXAMINING the .par files *********************

The smcgp-boolean.par file look like this:


5				 population_size
1.0				 per_cent_connection_mutate
1.0				 per_cent_function_mutate
1.0				 per_cent_argument_mutate
5000000  			 num_generations
1				 num_runs_total
100				 num_nodes
100				 levels_back
1				 progress_report
50000			         report_interval
173156789		         global_seed
0				 save_best_chrom
0				 run_from_chrom
targets.txt                      targetfiles
1				 ToDoLength
3				 num_iterations
1                                create_dot_files
0                                create_i0_files
0                                display_junk
 1  0    INCI-IO                   0
 0  0    DECI-IO                   1
 0  0    SKPI-IO                   2
 1  1    INCO-IO                   3
 0  1    DECO-IO                   4
 0  1    SKPO-IO                   5
 0  1    DEL-SM                    6
 0  1    ADD-SM                    7
 0  1    MOV-SM                    8
 0  1    OVR-SM                    9
 1  1    DUP-SM                   10
 0  1    CRP-SM                   11
 1  1    CHC-SM                   12
 0  1    CHF-SM                   13
 0  1    CHA-SM                   14
 0  0    CON-BOOL                 15
 0  0	 ZERO-BOOL		  16
 0  0	 ONE-BOOL	          17
 0  1	 WIRA-BOOL	          18
 0  1	 WIRB-BOOL		  19
 0  1	 NOTA-BOOL		  20
 0  1    NOTB-BOOL		  21
 1  2	 AND1-BOOL		  22
 0  2	 AND2-BOOL                23
 0  2	 AND3-BOOL                24
 1  2	 NOR-BOOL                 25
 0  2	 XOR-BOOL                 26
 0  2	 XNOR-BOOL                27
 1  2	 OR-BOOL                  28
 0  2	 OR1-BOOL                 29
 0  2	 OR2-BOOL                 30
 1  2	 NAND-BOOL                31
 0  3	 MUX1-BOOL                32
 0  3	 MUX2-BOOL                33
 0  3	 MUX3-BOOL                34
 0  3	 MUX4-BOOL                35


Each node in smcgp has the following genes: one function, a number of connections, three arguments

The number of function genes is the same as the number of nodes (num_nodes above). The number of connection genes
depends on the maximum arity of your chosen functions. Thelist of 36 functions are given in the .par file. You choose
which functions will be used using the first number in each function row. For example in the above the following primitive functions
have been chosen: INCI, INCO, DUP, CHC, AND1, NOR, OR, NAND. All the other possible node functions are disallowed. The second number
in each function primitive row is the arity of the functions (how many inputs they have). After the text description of the function
(which must have a dash in it to delimit the name, see later when I talk about how to generate graph pictures of the programs) there
is a number. This is the function gene.

So in the case above the maximum arity of each node is 2, so the program assumes that each node has two function genes. The computation performed by a node only uses the individual nodes arity. So technically IO functions ignore their connection genes.

Thus in the above case per_cent_connection_mutate is 1.0%. This means that when a genotype is mutated, only 2 connection genes
in the genotype will be mutated. In the above per_cent_function_mutate and per_cent_argument_mutate are also 1.0%, so this means
that only 1 function gene will be mutated. Since each node is assumed to have three aruments, 1% means 3 argument genes will be mutated.

Most of the parameters are standard to CGP. However the smcgp program requires some additional parameters. 

You need to have a file (here targets.txt) which says what problems you want solved at each iteration.
targets.txt              targetfiles  

For example targets could look like:

6        num_target_files
epar3.plu
epar4.plu
epar5.plu
epar6.plu
epar7.plu
epar8.plu

In this case num_iterations (see below) will have to be set to num_target_files -1 (first iteration is 0)

the target files (whatever the datatype used - see smcgp.h) themselves should have the format

.i 3
.o 1
.p n
n lines like this
inputnumber1 inputnumber2 inputnumber3 outputnumber1
.e

where if the datatype is either unsigned int, int or double. In the case of Boolean it is unsigned int.

the number after ".i" is how many inputs the target problem has
the number after ".o" is how many outputs the target problem has
the number after ".p" is how many inpit-output conditions the target problem has

ToDo is your chosen length for the ToDo self-modification list

1				 ToDo

num_iterations is how many iterations of smcgp you want to do. 

5				 num_iterations

(i.e. for 6 target files you need 5 iterations - the first one is assumed to be iteration 0)

You can set num_iterations to any integer less than num_target_files.

1                                create_dot_files
0                                create_i0_files
0                                display_junk

When create_dot_files is set to 1 (as above) the program will create a graph description in the form of a .dot file which can be read
into a free software package called graphviz (http://www.graphviz.org/). The program will create a .dot file at the end of
each evolutionary run. Once you convert it to pdf with graphviz you can view a .pdf file which will show the evolved graph.

When create_io_files is set to 1 the program will generate .dat files showing you what the evolved program outputs for each
data input. This is useful to check that the program is working properly and in case you want to plot what the output looks like
(e.g. in the case of evolving symbolic expressions).

display_junk allows you to decide whether you want the unconnected junk nodes to be included in the graph picture or not. If included
junk nodes are displayed in grey.

Input gathering and output producing functions

 1  0    INCI                   0
 0  0    DECI                   1
 0  0    SKPI                   2
 1  1    INCO                   3
 0  1    DECO                   4
 0  1    SKPO                   5


Self-modification functions

 0  1    DEL                    6
 0  1    ADD                    7
 0  1    MOV                    8
 0  1    OVR                    9
 1  1    DUP                   10
 0  1    CRP                   11
 1  1    CHC                   12
 0  1    CHF                   13
 0  1    CHA                   14

Spare function that can can return a constant
 0  0    CONST                 15

Boolean functions

 0  0	   ZERO			 16
 0  0	   ONE			 17
 0  1	   WIREA	         18
 0  1	   WIREB		 19
 0  1	   NOTA			 20
 0  1      NOTB		         21
 1  2	   AND1			 22
 0  2	   AND2                  23
 0  2	   AND3                  24
 1  2	   NOR                   25
 0  2	   XOR                   26
 0  2	   EXNOR                 27
 1  2	   OR                    28
 0  2	   OR1                   29
 0  2	   OR2                   30
 1  2	   NAND                  31
 0  3	   MUX1                  32
 0  3	   MUX2                  33
 0  3	   MUX3                  34
 0  3	   MUX4                  35

The symbolic regression program works with a .par file that looks like this:

        
5				 population_size
1.0				 per_cent_connection_mutate
1.0				 per_cent_function_mutate
1.0				 per_cent_argument_mutate
200000			 num_generations
1				 num_runs_total
100				 num_nodes
100				 levels_back
1				 progress_report
10000			       report_interval
173156789		       global_seed
0				 save_best_chrom
0				 run_from_chrom
symb-reg-targets.txt      targetfiles
1				 ToDo
1				 num_iterations
1                        create_dot_files
0                        create_i0_files
0                        display_junk
 1  0   INCI-IO               0
 0  0   DECI-IO               1         
 0  0   SKPI-IO               2
 1  1   INCO-IO               3
 0  1   DECO-IO               4
 0  1   SKPO-IO               5
 1  1   DEL-SM                6
 0  1   ADD-SM                7
 0  1   MOV-SM                8
 0  1   OVR-SM                9
 1  1   DUP-SM               10
 0  1   CRP-SM               11
 0  1   CHC-SM               12
 0  1   CHF-SM               13
 0  1   CHP-SM               14
 0  0   CON-MATH             15
 0  1	  MOD-MATH             16
 0  1	  SQRT-MATH            17
 0  1	  REC-MATH             18
 0  1	  SIN-MATH             19
 0  1	  COS-MATH             20
 0  1	  TAN-MATH             21
 0  1	  EXP-MATH             22
 0  1	  SINH-MATH            23
 0  1	  COSH-MATH            24
 0  1	  TANH-MATH            25
 0  1	  LOGE-MATH            26
 0  1	  LG10-MATH            27
 0  2	  SIN+-MATH            28
 0  2	  COS+-MATH            29
 0  2	  HYP-MATH             30
 0  2	  POW-MATH             31
 1  2	  ADD-MATH             32
 1  2	  SUB-MATH             33
 1  2	  MUL-MATH             34
 1  2	  DIV-MATH             35

Most of it is the same, but the primitive functions are different.

The file symb-reg-targets.txt looks like

2        num_target_files
quintic-polynomial.dat
sextic-polynomial.dat

So when run with the program, this requires two data files 


quintic-polynomial.dat and sextic-polynomial.dat

These data files are sampled from quintic and sextic polynomial functions
described by John Koza. They are just examples. However using the smcgp program
you can try to evolve programs that when iterated computer whatever functions you want.

So in principle this could be a long sequence of different functions!

The smcgp program is just a starting point and it could be adapted to solve many other problems.

For instance in our published smcgp papers we have got it to computer pi and e, find formulae
for Fibonacci series and many other problems.

It is worth noting that if the number of iteration is set to 0, and all the self-modifying functions
are set to zero (not used) the program effecitively becomes a form of standard cgp (albeit using
different ways of geting inputs and providing outputs).

Julian Miller, August 5th, 2013
University of York