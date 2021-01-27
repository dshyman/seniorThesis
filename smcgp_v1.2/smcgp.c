/* smcgp.c
   Julian F. Miller (c), 2013
   version 1.2:
   Dept. of Electronics, University of York, UK
*/

#include <stdio.h>
#include "smcgp.h"


int main(int argc,char* argv[])
{
	char parfile[MAX_NUM_LETTERS];

	validate_command_line(argc,argv,parfile);

	get_parameters(parfile);

	write_cgp_info(argv[0]);

	run_EA(num_runs_total);

	puts("");
	puts("*********    SMCGP COMPLETED      *********");

   return 0;
}

