/* smcgp-functions.c
   Julian F. Miller (c) 2013
   version 1.2
   Dept. of Electronics, University of York, UK

   The genotype is organized into nodes (structures)
   The connection genes are now relative (i.e. they
   count back from the node) unlike CGP where they are absolute
   There are no inputs or outputs now, instead there are
   input and output functions.

   IMPORTANT: For Boolean unsigned problems
   program outputs are arranged most significant on the left
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "smcgp.h"


/* validate cgp program command line and create file strings for .par and data files */
void validate_command_line(int argc, char* argv[], char parfile[])
{
	puts("");
	puts("*********    WELCOME TO SELF_MODIFYING CARTESIAN GENETIC PROGRAMMING      *********");
	puts("********* Validating command line arguments to smcgp program *********");

	if (argc!=2)
	{
		puts("INCORRECT NUMBER OF ARGUMENTS");
		puts("Type smcgp <file.par> then return");
		exit(1);
	}
	if (strlen(argv[1])>(MAX_NUM_LETTERS-1))
	{
   		puts("filename for parameter file is too long");
		printf("It should be less than %d characters\n",MAX_NUM_LETTERS);
		exit(2);
	}
	strcpy(parfile,argv[1]);

}


/* calculates how many mutations to do per chromosome */
int get_num_mutant(int num_genes, double per_cent_mutate)
{
	return (int)((double)num_genes*per_cent_mutate/100.0);
}

/* determines the number of mutations required
   for functions, connections and arguments
*/
void get_number_of_mutations(void)
{
	num_connection_mutations = get_num_mutant(num_connection_genes,per_cent_connection_mutate);
	num_function_mutations   = get_num_mutant(num_nodes,per_cent_function_mutate);
	num_argument_mutations   = get_num_mutant(num_argument_genes,per_cent_argument_mutate);
}

/* error checks the globals parameters read in
   get_parameters
*/
void check_parameter_error(void)
{
	FILE *fp;

	if (population_size > MAX_NUM_CHROMOSOMES)
	{
		printf("Too large a population size (<= %d)\n",MAX_NUM_CHROMOSOMES);
		exit(0);
	}

	if (num_connection_mutations <= 0)
	{
		printf("Your chosen mutation rate of %7.3lf for per_cent_connection_mutate\n",per_cent_connection_mutate);
		printf("means that zero connection mutations will take place.\n");
		printf("This is probably an error so the program is shutting down\n");
		exit(0);
	}

	if (num_function_mutations <= 0)
	{
		printf("Your chosen mutation rate of %7.3lf for per_cent_function_mutate\n",per_cent_function_mutate);
		printf("means that zero function mutations will take place.\n");
		printf("This is probably an error so the program is shutting down\n");
		exit(0);
	}

	if (num_argument_mutations <= 0)
	{
		printf("Your chosen mutation rate of %7.3lf for per_cent_argument_mutate\n",per_cent_argument_mutate);
		printf("means that zero argument mutations will take place.\n");
		printf("This is probably an error so the program is shutting down\n");
		exit(0);
	}


	if (num_runs_total < 1)
	{
		puts("Number of runs of EA must be at least 1");
		exit(0);
	}
	else if (num_runs_total > MAX_NUM_RUNS)
	{
		printf("Number of runs of EA must be less than %d\n", MAX_NUM_RUNS);
		exit(0);
	}

	if (levels_back > num_nodes)
	{
		printf("ERROR. levels-back exceeeds the number of nodes\n");
		exit(0);

	}

	if ((progress_report< 0) || (progress_report > 1))
	{
		puts("Progress report parameter must be 0 or 1");
		exit(0);
	}

	fp = fopen(targetfiles,"r");
	if (!fp)
	{
		printf("Cannot find file %s which lists the target problems\n",targetfiles);
		exit(0);
	}
	fclose(fp);



    if (ToDoLength > 10)
	{
		printf("ToDoLength of %d too big (> 10). Please reduce and re-run\n",ToDoLength);
		exit(0);

	}

	if (num_functions < 1)
	{
		puts("You have not selected any functions");
		exit(0);
	}

	if (num_input_functions == 0)
	{
		puts("You have not selected any input functions");
		exit(0);

	}

	if (num_output_functions == 0)
	{
		puts("You have not selected any output functions");
		exit(0);

	}
}

/* gives name to the functions that are nice for printing in a graph (dot).
   This takes the names in the .par file and looks for the letters on the
   left of the delimiter (dash) */
void extract_node_function_string(char string[MAX_NUM_LETTERS])
{
    int i = 0;
    char new_string[MAX_NUM_LETTERS];
    char delimiter = '-';

    do
    {
        new_string[i] = string[i];
        i++;
        if (string[i] == delimiter) /* jump out when encountering the delimiter in node name */
            break;

    }
    while (i <= (int) strlen(string));
    new_string[i]='\0';

    strcpy(string, new_string);
}

/* read from the parameter file all the global parameters */
void get_parameters(char parfile[MAX_NUM_LETTERS])
{
	int		i;
	int     arity, max_arity;
	char	dummy[50];
	FILE*	fp;

	printf("\n********* Reading parameters defined in %s *********\n",parfile);
	fp=fopen(parfile,"r");
	if (!fp)
	{
		printf("Missing file: %s\n",parfile);
		exit(1);
	}
	fscanf(fp,"%d %s", &population_size,dummy);
	fscanf(fp,"%lf %s",&per_cent_connection_mutate,dummy);
	fscanf(fp,"%lf %s",&per_cent_function_mutate,dummy);
	fscanf(fp,"%lf %s",&per_cent_argument_mutate,dummy);
	fscanf(fp,"%d %s", &num_generations,dummy);
	fscanf(fp,"%d %s", &num_runs_total,dummy);
	fscanf(fp,"%d %s", &num_nodes,dummy);
	fscanf(fp,"%d %s", &levels_back,dummy);
	fscanf(fp,"%d %s", &progress_report,dummy);
	fscanf(fp,"%d %s", &report_interval,dummy);
	fscanf(fp,"%u %s", &global_seed,dummy);
	fscanf(fp,"%d %s", &save_best_chrom,dummy);
	fscanf(fp,"%d %s", &run_from_chrom,dummy);
	fscanf(fp,"%s %s", targetfiles,dummy);
	fscanf(fp,"%d %s", &ToDoLength,dummy);
	fscanf(fp,"%d %s", &num_iterations,dummy);
	fscanf(fp,"%d %s", &create_dot_files,dummy);
	fscanf(fp,"%d %s", &create_io_files,dummy);
	fscanf(fp,"%d %s", &display_junk,dummy);

	num_functions = 0;
	num_output_functions = 0;
	num_input_functions = 0;
	max_arity = 0;

	for (i = 0; i < MAX_NUM_FUNCTIONS; i++)
	{
		/* number[] holds whether the function is used or not */
		fscanf(fp,"%d%d%s%s",&number[i],&arity,node_types[i],dummy);
		extract_node_function_string(node_types[i]);
		if (number[i])
		{
			allowed_functions[num_functions][0]=i;
			allowed_functions[num_functions][1]=arity;
			num_functions++;
			if (i<= 2)
			{
				num_input_functions++;
			}
			else if ((i >=3) && (i <= 5))
			{
				allowed_output_functions[num_output_functions] = i;
				num_output_functions++;
			}

			if (arity > max_arity)
				max_arity = arity;

		}
	}
	fclose(fp);

	/* each node is assigned max_arity connection genes */
	num_connection_genes_per_node = max_arity;
	num_connection_genes = num_nodes*num_connection_genes_per_node;
	num_argument_genes = num_nodes*MAX_NUM_ARGUMENTS;
    num_genes_per_node = num_connection_genes_per_node+MAX_NUM_ARGUMENTS+1;
    num_genes = num_nodes*num_genes_per_node;

    get_number_of_mutations();

	check_parameter_error();

	/* get input data size and gets target_file_names */
	read_target_file_info(&max_num_inputs,
					      &max_num_outputs,
		                  &max_num_tests);

	/* allocate memory for data_inputs and outputs */
	data_inputs = create_3d_datatype_space(max_num_inputs,
										   max_num_tests,num_iterations+1);


	data_outputs = create_3d_datatype_space(max_num_outputs,
										    max_num_tests,
										    num_iterations+1);

	/* get data from all files in targetfiles.txt */
	for (i = 0; i <= num_iterations; i++)
	{
		read_data(target_file_names[i], i);
	}
	/* calculate the perfect score i.e. all test cases in each target problem
	   are correct
    */
	perfect = define_perfect();

	srand(global_seed);

	puts("********* Parameters read *********");
	puts("********* Beginning execution *********");
}

/* write out parameter values in results file */
void write_cgp_info(char command[])
{
	int		i;
	FILE*	fp;

	fp=fopen("smcgp.txt","w");
	fprintf(fp,"The program is                   %s\n",command);
	fprintf(fp,"population_size is               %d\n",population_size);
	fprintf(fp,"connection mutation rate is      %6.2lf\n",per_cent_connection_mutate);
	fprintf(fp,"function mutation rate is        %6.2lf\n",per_cent_function_mutate);
	fprintf(fp,"argument mutation rate is        %6.2lf\n",per_cent_argument_mutate);
	fprintf(fp,"num_generations is               %d\n",num_generations);
	fprintf(fp,"num_runs is                      %d\n",num_runs_total);
	fprintf(fp,"num_nodes is                     %d\n",num_nodes);
	fprintf(fp,"levels_back is                   %d\n",levels_back);
	fprintf(fp,"progress report is               %d\n",progress_report);
	fprintf(fp,"report interval is               %d\n",report_interval);
	fprintf(fp,"global_seed is                   %u\n",global_seed);
	fprintf(fp,"save_best_chrom is               %d\n",save_best_chrom);
	fprintf(fp,"run_from_chrom is                %d\n",run_from_chrom);
	fprintf(fp,"target_files is				     %s\n",targetfiles);
	fprintf(fp,"ToDoLength list length is        %d\n",ToDoLength);
	fprintf(fp,"num_iterations is                %d\n",num_iterations);
	fprintf(fp,"create_dot_files is              %d\n",create_dot_files);
	fprintf(fp,"create_io_files is               %d\n",create_io_files);
	fprintf(fp,"display_junk is                  %d\n",display_junk);

	for (i=0;i<MAX_NUM_FUNCTIONS;i++)
	{
		fprintf(fp,"%d %s %d\n",number[i],node_types[i],i);
	}
	fprintf(fp,"\nHere are the Results\n");
	fclose(fp);
}

/*  returns a random integer between 0 and range-1 */
int newrand(int range)
{
    int temp;

    temp=rand() % range;
    return(temp);
}



/* read an item of data from problem
   specification file. Format of data
   depends on defined data type
   This will need to be edited if the user defines
   a new data type
*/
data_type myfscanf(FILE* fp)
{
	data_type datum_read;

	#ifdef DATA_IS_UNSIGNED_INT
		fscanf(fp,"%lu", &datum_read);
	#endif

	/* this converts input data immediately to double */
	#ifdef DATA_IS_INT
		fscanf(fp,"%lf", &datum_read);
	#endif

	#ifdef DATA_IS_DOUBLE
		fscanf(fp,"%lf", &datum_read);
	#endif

	return datum_read;
}

/* print an item of data of type data_type */
void myfprintf(FILE* fp, data_type datum, int linefeed)
{

    if (!linefeed)
    {
        #ifdef DATA_IS_UNSIGNED_INT
            fprintf(fp,"%lu\t", datum);
        #endif

        /* this converts input data immediately to double */
        #ifdef DATA_IS_INT
            fprintf(fp,"%lf\t", datum);
        #endif

        #ifdef DATA_IS_DOUBLE
            fprintf(fp,"%lf\t", datum);
        #endif
    }
    else
        fprintf(fp,"\n");
}


/* gets the number of inputs and outputs and tests
   from a data file and stores them in arrays.
*/
void read_data_header(char datafile[MAX_NUM_LETTERS],
					  int iteration)
{
    char	dummy[MAX_NUM_LETTERS];
    FILE*	fp;

    fp=fopen(datafile,"r");
    if (!fp)
    {
        puts("ERROR. Missing input data file (e.g. .plu (compressed Boolean) or .dat (data)");
        exit(1);
    }
    else
    {
        /* store the number of inputs and outputs and data points in arrays */
		fscanf(fp,"%s %d",dummy,&num_inputs[iteration]);
		fscanf(fp,"%s %d",dummy,&num_outputs[iteration]);
		fscanf(fp,"%s %d",dummy,&num_tests[iteration]);
		fclose(fp);
    }
}

/* find the maximum value in an integer array */
int get_max_array(int array[], int num_items)
{
	int i, max;

	max = array[0];
	for (i = 0; i < num_items; i++)
	{
		if (array[i] > max)
			max = array[i];
	}

	return max;

}

/* this function reads a file (e.g. targetfiles.txt)
   which lists the target problems. It analyzes
   these problems and determines how many data files are
   required and their number of inputs and outputs and test cases
   it calculates the maximum number of inputs, outputs and test cases
   so that the arrays data_inputs and data_outputs
   can be dimensioned. This is wasteful but it will do for now.
*/
void read_target_file_info(int* max_num_inputs,
						   int* max_num_outputs,
						   int* max_num_tests)
{
	int i, num_target_files;
	char dummy[MAX_NUM_LETTERS];
	FILE *fp;

	/* the string targetfiles is a filename string which is a
	   global read from the .par file */
	fp = fopen(targetfiles,"r");

	if (fp == NULL)
	{
		printf("ERROR. Cannot find file %s. This defines the target problems.\n",targetfiles);
		exit(0);
	}
	else
	{
		/* read targetfiles */

		fscanf(fp,"%d%s", &num_target_files,dummy);

		i = 0;
		do
		{
			fscanf(fp, "%s", target_file_names[i]);
			i++;

		}
		while (!feof(fp));

		i--;
		if (i != num_target_files)
		{
			printf("\nERROR. The number of files names %s in %d does not equal the number of target files\n",targetfiles,num_target_files);
			exit(0);

		}

		if (num_iterations+1 > num_target_files)
		{
			printf("\nERROR. The number of files names in %s must be one more than the number of iterations (%d)\n",targetfiles,num_iterations);
			exit(0);
		}

		/* find the number of inputs, outputs and tests
		   the arrays num_inputs, num_outputs and num_tests are
		   globals */
		for (i = 0; i <= num_iterations; i++)
			read_data_header(target_file_names[i], i);

		/* determine the maximum values */
		*max_num_inputs  =  get_max_array(num_inputs, num_iterations+1);
		*max_num_outputs =  get_max_array(num_outputs, num_iterations+1);
		*max_num_tests   =  get_max_array(num_tests, num_iterations+1);
	}
}


/* For each iteration reads input and output data from a file(s)
   and stores it in data_inputs and data_outputs
*/
void read_data(char datafile[MAX_NUM_LETTERS], int iteration)
{
    int		i,j, temp;
    char	dummy[MAX_NUM_LETTERS];
    FILE*	fp;

    fp=fopen(datafile,"r");

    if (!fp)
    {
        puts("ERROR. Missing input data file (e.g. .plu (compressed Boolean) or .dat (data)");
        exit(1);
    }
    else
    {
		fscanf(fp,"%s %d",dummy,&temp);
		fscanf(fp,"%s %d",dummy,&temp);
		fscanf(fp,"%s %d",dummy,&temp);

		for (i = 0; i < num_tests[iteration]; i++)
		{
			for(j = 0; j < num_inputs[iteration]; j++)
				data_inputs[iteration][i][j] = myfscanf(fp);
			for(j = 0; j < num_outputs[iteration]; j++)
				data_outputs[iteration][i][j] = myfscanf(fp);
		}
		fclose(fp);
    }
}

/* calculates a perfect score (if any)
   */
double define_perfect(void)
{
	int i;
	double perfect = 0.0;

	#ifdef DATA_IS_UNSIGNED_INT /* Boolean case: using compressed truth tables */
		for (i = 0; i <= num_iterations; i++)
		{
			perfect_at_it[i] = (double) pow2(num_inputs[i])*num_outputs[i];
            perfect = perfect + perfect_at_it[i];
		}
	#endif

	#ifdef DATA_IS_INT
		for (i = 0; i <= num_iterations; i++)
		{
			perfect_at_it[i] = (double) num_tests[i];
            perfect = perfect + perfect_at_it[i];
		}

	#endif

	#ifdef DATA_IS_DOUBLE
		for (i = 0; i <= num_iterations; i++)
		{
            perfect_at_it[i] = (double) num_tests[i]*(1-ERROR_THRESHOLD);
            perfect = perfect + perfect_at_it[i];
		}

	#endif


	return perfect;
}


/* prints a chromosome to a file
   when append is 1, the function appends the information to the file
   when append is 0, the function creates a new file
*/
void fprint_a_chromosome(node* chromosome,
						 int   num_nodes_chromosome,
						 char  name[],
						 int   append)
{
	int		i, j;
	FILE*	fp;

	if (append)
		fp=fopen(name,"a");
	else
	   fp=fopen(name,"w");

	fprintf(fp,"\ngenotype is\n");

	for (i = 0; i < num_nodes_chromosome; i++)
	{
		fprintf(fp,"[%s]",node_types[chromosome[i].function]);

		fprintf(fp,"(");

		for (j = 0; j < num_connection_genes_per_node-1; j++)
			fprintf(fp,"%d,",chromosome[i].connection[j]);

		fprintf(fp,"%d",chromosome[i].connection[num_connection_genes_per_node-1]);

		fprintf(fp,": ");

		for (j = 0; j < MAX_NUM_ARGUMENTS-1; j++)
			fprintf(fp,"%d,",chromosome[i].argument[j]);

		fprintf(fp,"%d",chromosome[i].argument[MAX_NUM_ARGUMENTS-1]);

		fprintf(fp,"):%d\n",i);
   }

   fprintf(fp,"\n\n");
   fclose(fp);
}

/* prints a phenotype to a file
   when append is 1, the function appends the information to the file
   when append is 0, the function creates a new file
*/
void fprint_a_phenotype(node phenotype[MAX_NUM_NODES],
						int   num_nodes_phenotype,
						char  name[],
						int   append)
{
	int		i, j;
	FILE*	fp;

	if (append)
		fp=fopen(name,"a");
	else
	   fp=fopen(name,"w");

	fprintf(fp,"\n\nAll nodes in phenotype (active and inactive) are given below.\n");
	fprintf(fp,"Format is [node function](connections : arguments): node position\n");

	for (i = 0; i < num_nodes_phenotype; i++)
	{
		fprintf(fp,"[%s]",node_types[phenotype[i].function]);

		fprintf(fp,"(");

		for (j = 0; j < num_connection_genes_per_node-1; j++)
			fprintf(fp,"%d,",phenotype[i].connection[j]);

		fprintf(fp,"%d",phenotype[i].connection[num_connection_genes_per_node-1]);

		fprintf(fp,": ");

		for (j = 0; j < MAX_NUM_ARGUMENTS-1; j++)
			fprintf(fp,"%d,",phenotype[i].argument[j]);

		fprintf(fp,"%d",phenotype[i].argument[MAX_NUM_ARGUMENTS-1]);

		fprintf(fp,"):%d\n",i);
   }

   fprintf(fp,"\n\n");
   fclose(fp);
}

/* prints a chromosome to the screen */
void print_a_chromosome(node* chromosome, int num_nodes_chromosome)
{
   int i, j;

   printf("\n");

   for (i = 0; i < num_nodes_chromosome; i++)
   {

	  printf("%d   ",chromosome[i].function);

      for (j = 0; j< num_connection_genes_per_node; j++)
         printf("%d ",chromosome[i].connection[j]);

	  printf("   ");

      for (j = 0; j< MAX_NUM_ARGUMENTS; j++)
         printf("%d ",chromosome[i].argument[j]);


      printf("\n");

   }

   printf("\n");

}

/* prints a chromosome to the screen */
void print_a_phenotype(node phenotype[MAX_NUM_NODES],
					   int num_nodes_phenotype)
{
   int i, j;

   printf("\n");

   for (i = 0; i < num_nodes_phenotype; i++)
   {

	  printf("%d   ",phenotype[i].function);

      for (j = 0; j< num_connection_genes_per_node; j++)
         printf("%d ",phenotype[i].connection[j]);

	  printf("   ");

      for (j = 0; j< MAX_NUM_ARGUMENTS; j++)
         printf("%d ",phenotype[i].argument[j]);


      printf("\n");

   }

   printf("\n");

}



/* prints a chromosome to file in raw format,
  so that it can be read by the program
  the format is: function,  connections,  arguments etc
*/
void fprint_a_raw_chromosome(node* chromosome,
							 int   num_nodes_chromosome,
							 char  name[MAX_NUM_LETTERS], int append)
{
   int i, j;
   FILE* fp;

   if (append)
	   fp = fopen(name, "a");
   else
	   fp = fopen(name, "w");

   fprintf(fp,"\n");

   for (i = 0; i < num_nodes_chromosome; i++)
   {

	  fprintf(fp,"%d   ",chromosome[i].function);

      for (j = 0; j< num_connection_genes_per_node; j++)
         fprintf(fp,"%d ",chromosome[i].connection[j]);

	  fprintf(fp,"   ");

      for (j = 0; j< MAX_NUM_ARGUMENTS; j++)
         fprintf(fp,"%d ",chromosome[i].argument[j]);


      fprintf(fp,"\n");

   }

   fprintf(fp,"\n");

   fclose(fp);
}

/* determines from the function gene what its arity is.
   It uses the global arrays allowed_functions */
int get_arity(int function_gene)
{
	int i, arity = 0;

	for (i = 0; i < num_functions; i++)
		if (function_gene == allowed_functions[i][0])
			arity = allowed_functions[i][1];

	return arity;
}


/* determine which nodes in the chromosome will act as program outputs.
   This is important as it is these which determine the chain of nodes
   that define the phenotype. However, it is possible that some
   of the output nodes are still redundant as whatever they write to the
   output register gets overwritten by another function. Strictly speaking
   we should determine which output nodes last wrote to each position
   in the output register. It is very complicated as output nodes
   write and then change the output register writing head, so even
   if they wrote a value that was overwritten they could affect subsequent
   output writing functions because they change the position of the
   writing head!
   The output addresses hold the addresses in the chromosome
   where the outputs are taken from (like cgp output genes).
*/
int get_output_nodes_genotype(node* chromosome,
					          int   num_nodes_chromosome,
					          int   output_nodes[MAX_NUM_NODES])
{
	int i;
	int num_output_nodes = 0;

	/* find out how many active OUTPUT nodes there are */
	/* read the nodes from the right */
	for (i = num_nodes_chromosome - 1; i >= 0; i--)
	{
		if (function_type_is_output(chromosome[i].function))  /* the function is an output */
		{
			output_nodes[num_output_nodes] = i; /*store output node addresses from right */
			num_output_nodes++;
		}
	}

	return num_output_nodes;
}

/* like get_output_nodes_genotype but works on phenotype */
int get_output_nodes_phenotype(int   num_nodes_phenotype,
					           int   output_nodes[MAX_NUM_NODES])
{
	int i;
	int num_output_nodes = 0;

	/* find out how many active OUTPUT nodes there are */
	/* read the nodes from the right */
	for (i = num_nodes_phenotype - 1; i >= 0; i--)
	{
		if (function_type_is_output(phenotype[i].function))  /* the function is an output */
		{
			output_nodes[num_output_nodes] = i; /*store output node addresses from right */
			num_output_nodes++;
		}
	}

	return num_output_nodes;
}


/* calculates a boolean array: node_used that indicates
   whether a node in the chromosome is active or not */
void get_node_used_in_genotype(node* chromosome,
				               int   num_nodes_chromosome,
				               int   node_used[MAX_NUM_NODES])
{
	int		i,j, num_output_nodes;
	int		node_genes[MAX_NUM_CONNECTIONS_PLUS_FUNCTION];
	int		output_nodes[MAX_NUM_NODES];


    /* analyze chromosome and determine where the outputs are coming from */
	num_output_nodes = get_output_nodes_genotype(chromosome,
		                                         num_nodes_chromosome,
		                                         output_nodes);
	/* initialise the node_used array */
	for (i = 0; i < num_nodes_chromosome; i++)
		node_used[i] = 0;

	/* all num_output_nodes are active */
	for (i = 0; i < num_output_nodes; i++)
		node_used[output_nodes[i]]=1;

	for (i = num_nodes_chromosome - 1; i >= 0; i--)
	{
		if (node_used[i])
		{
			/* write connection genes for node into array node_genes */
			for (j = 0; j < num_connection_genes_per_node;j++)
				node_genes[j] = chromosome[i].connection[j];

			node_genes[num_connection_genes_per_node] = chromosome[i].function;

			/* each function has an arity stored in
			   allowed_functions[][1].
			   Find the nodes whose data is used
			*/
			for (j = 0; j < get_arity(node_genes[num_connection_genes_per_node]); j++)
				if (node_genes[j] <= i) /* if nodes point within graph */
					node_used[i-node_genes[j]]=1;
		}
	}
}

/* calculates a boolean array: node_used that indicates
   whether a node in the phenotype is active or not */
void get_node_used_in_phenotype(int   num_nodes_phenotype,
				                int   node_used[MAX_NUM_NODES])
{
	int		i,j, num_output_nodes;
	int		node_genes[MAX_NUM_CONNECTIONS_PLUS_FUNCTION];
	int		output_nodes[MAX_NUM_NODES];


    /* analyze chromosome and determine where the outputs are coming from */
	num_output_nodes = get_output_nodes_phenotype(num_nodes_phenotype,
		                                          output_nodes);
	/* initialise the node_used array */
	for (i = 0; i < num_nodes_phenotype; i++)
		node_used[i] = 0;

	/* all num_output_nodes are active */
	for (i = 0; i < num_output_nodes; i++)
		node_used[output_nodes[i]]=1;

	for (i = num_nodes_phenotype - 1; i >= 0; i--)
	{
		if (node_used[i])
		{
			/* write connection genes for node into array node_genes */
			for (j = 0; j < num_connection_genes_per_node;j++)
				node_genes[j] = phenotype[i].connection[j];

			node_genes[num_connection_genes_per_node] = phenotype[i].function;

			/* each function has an arity stored in
			   allowed_functions[][1].
			   Find the nodes whose data is used
			*/
			/* What if relative address beyond left end of graph? */
			/* in this case the node input connects to zero, so we do nothing */
			for (j = 0; j < get_arity(node_genes[num_connection_genes_per_node]); j++)
				if (node_genes[j] <= i) /* if nodes point within graph */
					node_used[i-node_genes[j]]=1;
		}
	}
}

/* Boolean function that returns 1 or 0 depending on whether
   the_node of a chromosome is an sm function or not */
int node_is_sm(int the_node)
{
	if ((phenotype[the_node].function >= 6) && (phenotype[the_node].function <= 14))
		return 1;
	else
		return 0;
}

/* print active chromosome genes to file */
void fprint_active_genes(node* chromosome,
						 int   num_nodes_chromosome,
						 char  name[MAX_NUM_LETTERS])
{
	int		i,j;
	int		node_used[MAX_NUM_NODES];

	FILE*	fp;

	fp=fopen(name,"a");

	fprintf(fp,"\nActive genes at are:\n");

	get_node_used_in_genotype(chromosome,
		                     num_nodes_chromosome,
		                     node_used);

	fprintf(fp,"\nnode_used is\n");

	for (i = 0; i < num_nodes_chromosome; i++)
		fprintf(fp,"%2d: %d\n",i,node_used[i]);

	fprintf(fp,"\n");

	for (i = 0; i < num_nodes_chromosome; i++)
	{
		if (node_used[i])
		{

			fprintf(fp,"[%s]",node_types[chromosome[i].function]);

			fprintf(fp,"(");

			for (j = 0; j < num_connection_genes_per_node-1; j++)
				fprintf(fp,"%d,",chromosome[i].connection[j]);

			fprintf(fp,"%d",chromosome[i].connection[num_connection_genes_per_node-1]);

			fprintf(fp,": ");

			for (j = 0; j < MAX_NUM_ARGUMENTS-1; j++)
				fprintf(fp,"%d,",chromosome[i].argument[j]);

			fprintf(fp,"%d",chromosome[i].argument[MAX_NUM_ARGUMENTS-1]);

			fprintf(fp,"):%d\n",i);
		}
   }

   fprintf(fp,"\n\n");
   fclose(fp);
}


/* print active phenotype nodes (I call then phenes) to file */
void fprint_active_phenes(int   num_nodes_phenotype,
						  char  name[MAX_NUM_LETTERS],
						  int   iteration)
{
	int		i,j;
	int		node_used[MAX_NUM_NODES];

	FILE*	fp;

	fp=fopen(name,"a");

	fprintf(fp,"\nActive phenotype nodes at iteration %d are:\n",iteration);

	get_node_used_in_phenotype(num_nodes_phenotype,
		                       node_used);

	for (i = 0; i < num_nodes_phenotype; i++)
	{
		if (node_used[i])
		{

			fprintf(fp,"[%s]",node_types[phenotype[i].function]);

			fprintf(fp,"(");

			for (j = 0; j < num_connection_genes_per_node-1; j++)
				fprintf(fp,"%d,",phenotype[i].connection[j]);

			fprintf(fp,"%d",phenotype[i].connection[num_connection_genes_per_node-1]);

			fprintf(fp,": ");

			for (j = 0; j < MAX_NUM_ARGUMENTS-1; j++)
				fprintf(fp,"%d,",phenotype[i].argument[j]);

			fprintf(fp,"%d",phenotype[i].argument[MAX_NUM_ARGUMENTS-1]);

			fprintf(fp,"):%d\n",i);
		}
   }

   fprintf(fp,"\n\n");
   fclose(fp);
}

/* creates a dot file of the phenotype ready for graphviz dot program */
void fprint_dot_graph(int num_nodes_phenotype, int iteration, int run, char name[50])
{
	int		i,j;
    int		node_used[MAX_NUM_NODES];
    char    basename[50];
	char	dotfile[MAX_NUM_LETTERS];
	char    it[10] ={"_it"}, itstring[5];
	char    r[10] = {"_run"}, runstring[5];
	/* char   pdffile[MAX_NUM_LETTERS], cmd[1024]; */
	char    argstring[MAX_NUM_LETTERS],arg[5];
	FILE*	fp;

    for (i = 0; i < 5; i++)
        basename[i] = name[i];
    basename[5] = '\0'; /* remove .txt from name */

    sprintf(itstring, "%d", iteration);
    strcat(it, itstring);
    sprintf(runstring, "%d", run);
    strcat(r, runstring);

    /* Create filename for .dot file */
	strcpy(dotfile, basename);
	strcat(dotfile, r);
	/* put the iteration on the end */
	strcat(dotfile, it);
	/* strcpy(pdffile, name); */
	strcat(dotfile, ".dot");

	/* strcat(pdffile, ".pdf"); */

    /* to draw the graphs we need the array node_used */
	get_node_used_in_phenotype(num_nodes_phenotype,node_used);

	/* Open .dot file */
	fp=fopen(dotfile,"w");

	/* print start of graph */
    fprintf(fp,"digraph G\n{\n   rankdir = LR;\n");

	for (i = 0 ;i < num_nodes_phenotype; i++)
	{
        /* make argument string */
        strcpy(argstring, "");
        for (j = 0; j < MAX_NUM_ARGUMENTS; j++)
        {
            sprintf(arg, "%d ", phenotype[i].argument[j]);
            strcat(argstring,arg); /* TODO: Place these under SM nodes - how */
        }
        /* draw the function node */
        if (node_used[i])
            fprintf(fp,"\n\n   f%d [shape=circle  style=bold label=\"%s %d\"];",i, node_types[phenotype[i].function],  i);
        else
            if (display_junk)
                fprintf(fp,"\n\n   f%d [shape=circle color=gray style=bold label=\"%s %d\"];",i, node_types[phenotype[i].function], i);

        /* Look up arity of function and check the nodes connections */
		for (j = 0; j < get_arity(phenotype[i].function); j++)
        {
            if (i-phenotype[i].connection[j] < 0) /* input to node is a zero */
            {
                if (node_used[i])
                {
                    fprintf(fp,"\n   0 [shape=box label=\"0\"];");
                    fprintf(fp,"\n   0 -> f%d;", i);
                }
                else
                {
                    if (display_junk)
                    {
                        fprintf(fp,"\n   0 [shape=box color=gray label=\"0\"];");
                        fprintf(fp,"\n   0 -> f%d [color = gray];", i);
                    }

                }
            }
            else /* input to node is another node */
            {
                if (node_used[i])
                {
                    fprintf(fp,"\n   f%d -> f%d [headlabel = %d labeldistance = 2.0 fontsize = 10];", i-phenotype[i].connection[j], i,  i-phenotype[i].connection[j]);
                }
                else
                {
                    if (display_junk)
                        fprintf(fp,"\n   f%d -> f%d [headlabel = %d labeldistance = 2.0 color=gray fontsize = 10];", i-phenotype[i].connection[j], i, i-phenotype[i].connection[j]);
                }
            }
        }
	}

   /* print end of graph */
   fprintf(fp, "}\n");

   fclose(fp);
}

/* print out to a file the input-output mapping produced by smcgp at a given iteration */
void print_evolved_input_output_table(int num_nodes_to_process,
                                      int nodes_to_process[MAX_NUM_NODES],
                                      int run,
                                      int iteration)
{
	int			i, fitness_test;
	data_type	cgp_outputs[MAX_NUM_OUTPUTS];
	char        runstring[5], itstring[5];
	char        name[MAX_NUM_LETTERS];
	FILE        *fp;


	#ifdef DATA_IS_UNSIGNED_INT
		get_bit_width(iteration);
	#endif

	sprintf(runstring, "%d", run);
	sprintf(itstring, "%d", iteration);
	strcpy(name, "smcgp_io_run");
	strcat(name, runstring);
	strcat(name, "_it");
	strcat(name, itstring);
	strcat(name, ".dat");

	fp = fopen(name, "w");

    if (fp)
    {
        for (fitness_test = 0; fitness_test < num_tests[iteration]; fitness_test++)
        {
            decode_cgp(cgp_outputs,num_nodes_to_process, nodes_to_process,
                    fitness_test, iteration);

            for (i = 0; i < num_inputs[iteration]; i++)
            {
                myfprintf(fp, data_inputs[iteration][fitness_test][i], 0);
            }
            for (i = 0; i < num_outputs[iteration]; i++)
            {
                 /* mask off irrrelevant bits in small input problems */
                #ifdef DATA_IS_UNSIGNED_INT
                if (num_inputs[iteration] == 2)
                    cgp_outputs[i] = cgp_outputs[i] & 15;
                if (num_inputs[iteration] == 3)
                    cgp_outputs[i] = cgp_outputs[i] & 255;
                else if (num_inputs[iteration] == 4)
                    cgp_outputs[i] = cgp_outputs[i] & 65535;
                #endif
                myfprintf(fp, cgp_outputs[i], 0);
            }

            myfprintf(fp, 0, 1);
        }
    }
    else
    {
        printf("\nIn print_evolved_input_output_table: could not open %s for writing", name);
        exit(0);
    }
}

/* prints the genotype and the active genes and then
   generates all the phenotypes and prints out the
   phenes and the active phenes. If withdot = 0 then run
   is ignored */
void fprint_all_phenotypes(node* genotype,
						   char name[MAX_NUM_LETTERS],
						   int withdot, int run, int print_io)
{
	int     i, j;
	int     num_nodes_phenotype;
	ToDo    ToDoList[MAX_TODO];
	int     num_sm_nodes_to_process;
	int		num_nodes_to_process = 0;
	int     nodes_to_process[MAX_NUM_NODES];
	int     iteration = 0;
	double	fit = 0.0;
	double  fit_at_it;

	FILE*   fp;

	num_nodes_phenotype = num_nodes;

    /* make the initial phenotype a copy of the genotype */
    copy_chrom_to_phenotype(genotype,num_nodes_phenotype);

    do
    {
        num_nodes_to_process = get_nodes_to_process(num_nodes_phenotype,
                                                    ToDoList,
                                                    &num_sm_nodes_to_process,
                                                    nodes_to_process,
                                                    iteration);

        /* create dot files to draw graphs of phenotypes -using package graphviz */
        if (withdot)
        {
            fprint_dot_graph(num_nodes_phenotype, iteration, run, name);
        }



        if (iteration != num_iterations)
        {
            fp = fopen(name, "a");
            fprintf(fp, "\niteration is %d\n", iteration);
            fprintf(fp, "\nnum_nodes_phenotype is %d\n", num_nodes_phenotype);
            fprintf(fp, "\nnum_nodes_to_process is %d\n", num_nodes_to_process);
            fprintf(fp, "\nThese are the nodes_to_process \n");
            for (i = 0; i < num_nodes_to_process; i++)
                fprintf(fp, " %d ", nodes_to_process[i]);
            fprintf(fp,"\nnum_sm_nodes_to_process is %d\n", num_sm_nodes_to_process);
            if (num_sm_nodes_to_process != 0)
            {
                fprintf(fp,"\nThe ToDo list is\n");
                for (i = 0; i < num_sm_nodes_to_process; i++)
                {
                    fprintf(fp," (pos: %d: args: ", ToDoList[i].position);
                    for (j = 0; j < MAX_NUM_ARGUMENTS; j++)
                        fprintf(fp,"%d ", ToDoList[i].argument[j]);
                    fprintf(fp,"funct: %d)", ToDoList[i].function);
                }
            }
            fclose(fp);
        }

        fprint_a_phenotype(phenotype, num_nodes_phenotype, name, 1);
        fprint_active_phenes(num_nodes_phenotype, name, iteration);


        /* calculate the fitness of this phenotype */
        fit_at_it = fitness_at_iteration(num_nodes_to_process,
                                         nodes_to_process,
                                         iteration);

        fit = fit + fit_at_it;

        if (print_io)
        {
            print_evolved_input_output_table(num_nodes_to_process,nodes_to_process, run, iteration);
        }

        /* print out the fitness */
        fp = fopen(name, "a");
        fprintf(fp,"\nfit_at_it is %4.2lf and fit is %4.2lf\n", fit_at_it, fit);
        fclose(fp);

        iteration++;
        if (iteration != (num_iterations+1))
        {
            /* generate the new phenotype */
            if (num_sm_nodes_to_process  != 0)
            {
                generate_new_phenotype_verbose(&num_nodes_phenotype,
                                               ToDoList,
                                                num_sm_nodes_to_process,
                                                name,iteration);
            }
        }
    } /* end of iteration loop */
    while (iteration <= num_iterations);
}

/* copy a node */
node copy_node(node source)
{
	int  j;
	node the_copy;

	for (j = 0; j < num_connection_genes_per_node; j++)
		the_copy.connection[j] = source.connection[j];

	for (j = 0; j < MAX_NUM_ARGUMENTS; j++)
        the_copy.argument[j] = source.argument[j];

    the_copy.function = source.function;

	return the_copy;
}

/* copy a chromosome */
void copy_chrom(node* source,
				node* destination,
				int   number_of_nodes)
{
	int i;

	/* copy all nodes */
	for (i = 0; i < number_of_nodes; i++)
	   destination[i] = copy_node(source[i]);
}

/* copy a chromosome */
void copy_chrom_to_phenotype(node* source,
				             int   number_of_nodes)
{
	int i;

	/* copy all nodes */
	if (number_of_nodes > MAX_NUM_NODES)
	{
	    printf("\nIn copy_chrom_to_phenotype");
	    printf("\nnumber_of_nodes(%d) exceeds allowed maximum(%d)", number_of_nodes, MAX_NUM_NODES);
	    printf("\nEither there is a bug or MAX_NUM_NODES needs to be increased");
	    exit(0);
	}
	else
	{
        for (i = 0; i < number_of_nodes; i++)
            phenotype[i] = copy_node(source[i]);
	}

}


/* copy a phenotype */
void copy_phen(node source[MAX_NUM_NODES],
			   node destination[MAX_NUM_NODES],
			   int  number_of_nodes)
{
	int i;

	/* copy all nodes */
	if (number_of_nodes > MAX_NUM_NODES)
	{
	    printf("\nIn copy_phen");
	    printf("\nnumber_of_nodes(%d) exceeds allowed maximum(%d)", number_of_nodes, MAX_NUM_NODES);
	    printf("\nEither there is a bug or MAX_NUM_NODES needs to be increased");
	    exit(0);
	}
	else
	{
        for (i = 0; i < number_of_nodes; i++)
            destination[i] = copy_node(source[i]);
	}

}


/* copy a chromosome from the population into a chromsome */
void copy_pop_to_chrom(node** chromosomes,
					   node*  chromosome,
					   int    pop_member)
{
	int i;

	/* copy all nodes */
	for (i = 0; i < num_nodes; i++)
		chromosome[i] = copy_node(chromosomes[pop_member][i]);
}


/* copy a chromosome into the population */
void copy_chrom_to_pop(node** chromosomes,
					   node*  chromosome,
					   int    pop_member)
{
	int i;

	/* copy all nodes */
	for (i = 0; i < num_nodes; i++)
		chromosomes[pop_member][i] = copy_node(chromosome[i]);
}

/* generate a starting population
   from a chromosome read from a file (smcgp.chr)
*/
void read_from_chrom(node** chromosomes)
{
	int	  i, j, k, the_node;
	int   gene;
	FILE* fp;


    fp = fopen("smcgp.chr","r");
    if (!fp)
    {
		puts("Missing file smcgp.chr (contains a chromosome)");
        exit(1);
    }
    else
	{
        /* make starting population copies of loaded chromosome */
		for (j = 0; j < population_size; j++)
        {
			if (j==0)
            {
				i = 0;
				the_node = 0;
				k = 0;
				do
				{
					fscanf(fp,"%d",&gene);

					/* function gene is the first gene of a node */
					if (i % num_genes_per_node == 0)
					{
						chromosomes[j][the_node].function = gene;
						k = 0;
					}
					else
					{
						if (k < num_connection_genes_per_node) /* read connections */
						{
							chromosomes[j][the_node].connection[k] = gene;
						}
						else /* read arguments */
						{
							chromosomes[j][the_node].argument[k-num_connection_genes_per_node] = gene;
						}
						k++;
						if (k == num_connection_genes_per_node + MAX_NUM_ARGUMENTS)
							the_node++;

					}
					i++;
					if (i == num_genes)
						break;
				}
				while (!feof(fp));

				if (the_node!=num_nodes)
				{
					puts("ERROR. Number of nodes in smcgp.chr does not match the expected number");
					printf("\nnum_nodes required is %d, num_nodes read is %d\n",num_nodes, the_node);
					puts("Check the number of nodes in the .par file");
					exit(0);
				}
            }
            else
				copy_chrom_to_pop(chromosomes, chromosomes[0],j);
        }
		fclose(fp);
	}
}

/* This function analyzes the input node function and writes
   the data inputs of the node into the inputs array (in[0]).
   It also updates the input_pointer */
void in_nodes(int       function_type,
			  int*      input_pointer,
			  int       the_node,
			  int       nodes_to_process[MAX_NUM_NODES],
			  data_type in[MAX_NUM_CONNECTION_GENES_PER_NODE],
              int       iteration,
			  int       fitness_test)
{
	int inp;

	inp = *input_pointer;

	/* get input data from place indicated by inp */
	in[0] = data_inputs[iteration][fitness_test][inp];

	switch(function_type)
	{
		case 0:  /*   INCI (aka INP): increment input_pointer  */
			inp++;
			if (inp == num_inputs[iteration]) /* wrap around */
				inp = 0;
			break;
		case 1:  /*   DECI (aka INPP): decrement input_pointer  */
			inp--;
			if (inp == -1) /* wrap around */
				inp = num_inputs[iteration] - 1;
			break;
		case 2: /* SKPI (aka SKIPP): move input_pointer by argument[0] */
			inp = inp + phenotype[nodes_to_process[the_node]].argument[0];
            inp = inp % num_inputs[iteration]; /* make sure inp is valid */
			if (inp < 0) /* make inp valid */
				inp = inp + num_inputs[iteration];
	}
	*input_pointer = inp;
}

/* This function analyzes the output node function and writes
   the data inputs of the node into the node_outputs array and
   also the array in[] which stores the inputs to the node */
   /* output[] stores the outputs of nodes */
void out_nodes(int       function_type,
			   int*      output_pointer,
			   data_type output_register[MAX_NUM_OUTPUTS],
			   int       the_node_to_process,
			   int       nodes_to_process[MAX_NUM_NODES],
			   data_type in[MAX_NUM_CONNECTION_GENES_PER_NODE],
			   data_type node_output[MAX_NUM_NODES],
               int       iteration)
{
	int data_address;
	int the_node_we_are_on;
	int oup;

	oup = *output_pointer;

	the_node_we_are_on = nodes_to_process[the_node_to_process];
	/* go back this many nodes from where we are: relative addressing */
	/* the number of nodes to count back is stored in the first node connection */
	data_address = the_node_we_are_on - phenotype[the_node_we_are_on].connection[0];

	if (data_address < 0) /* left of first node */
		in[0] = 0;
	else
		in[0] = node_output[data_address];

	/* write first data input to output register */
	/* that is what output nodes do just write their input
	   to the output register */
	output_register[oup] = in[0];

    /* update output_pointer */
	switch(function_type)
	{
		case 3: /* INCO: increment output_pointer */
			oup++;
			if (oup == num_outputs[iteration])
				oup = 0;
			break;
		case 4:  /* DECO: decrement output_pointer */
			oup--;
			if (oup == -1)
				oup = num_outputs[iteration]-1;
			break;
		case 5:  /* SKPO: move output_pointer by argument[0] */
			oup = oup + phenotype[nodes_to_process[the_node_to_process]].argument[0];
			oup = oup % num_outputs[iteration];
			if (oup < 0)
				oup = oup + num_outputs[iteration];
	}

	*output_pointer = oup;
}

/* This function analyzes the sm node function and writes
   the data inputs of the node into the node inputs array */
void sm_nodes(int       the_node_to_process,
			  int       nodes_to_process[MAX_NUM_NODES],
			  data_type in[MAX_NUM_CONNECTION_GENES_PER_NODE],
			  data_type	output[MAX_NUM_NODES])
{
	int data_address;

    /* use data coming from first connection to sm node */
	/* this is not data dependent. If data dependency is wanted
	   then we need to analyze the actual data arriving
	   at the connections to the sm node and choose which
	   data will be passed on. */
	/* NOTE TO MYSELF: why not make all nodes data dependent?
	   this could be done by all nodes having two functions,
	   so math functions would have an if (data condition) choose
	   math1 else choose math 2. I need to think about this
	*/
	data_address = nodes_to_process[the_node_to_process]-phenotype[nodes_to_process[the_node_to_process]].connection[0];
	if (data_address < 0) /* left of first node */
		in[0] = 0;
	else
		in[0]=output[data_address];
}

/* This function analyzes the comp node function and writes
   the data inputs of the node into the node inputs array */
void comp_nodes(int       function_type,
				int       the_node,
			    int       nodes_to_process[MAX_NUM_NODES],
			    data_type in[MAX_NUM_CONNECTION_GENES_PER_NODE],
				data_type output[MAX_NUM_NODES])
{
	int i, data_address;

	if (function_type == 15) /*CONST*/
	{
		in[0] = phenotype[nodes_to_process[the_node]].argument[0];
		#ifdef DATA_IS_UNSIGNED_INT
			in[0] = in[0] % 2; /* either 0 or 1 */
		#endif
	}
	else
	{
		for (i = 0; i < num_connection_genes_per_node; i++)
		{   /* get input data to node by relative address */
			data_address = nodes_to_process[the_node]-phenotype[nodes_to_process[the_node]].connection[i];
			if (data_address < 0) /* left of first node */
				in[i] = 0;
			else
				in[i]=output[data_address];
		}
	}
}

/* Boolean function that returns 1 when the function_gene refers to
   an input gathering function */
int function_type_is_input(int function_gene)
{
    return (function_gene <= 2);
}

/* Boolean function that returns 1 when the function_gene refers to
   an input gathering function */
int function_type_is_output(int function_gene)
{
    return (function_gene >= 3) && (function_gene <= 5);
}

/* Boolean function that returns 1 when the function_gene refers to
   an sm function */
int function_type_is_sm(int function_gene)
{
    return (function_gene >= 6) && (function_gene <= 14);
}

/* Boolean function that returns 1 when the function_gene refers to
   a computational function */
int function_type_is_comp(int function_gene)
{
    return (function_gene >= 15);
}

/* this decodes the cgp chromosome.
   It is given data_inputs corresponding to a single
   test case and it calculates what the cgp phenotype gives
   for the data outputs (cgp_outputs)
   It only processes nodes that are used (nodes_to_process)
*/
void decode_cgp(data_type cgp_outputs[MAX_NUM_OUTPUTS],
				int num_nodes_to_process,
				int nodes_to_process[MAX_NUM_NODES],
				int fitness_test,
				int iteration)
{
	int			i, the_node;
	int			function_type;
	int         input_pointer = 0;
	int         output_pointer = 0;
	data_type	in[MAX_NUM_CONNECTION_GENES_PER_NODE];
	data_type	output[MAX_NUM_NODES]={0ul};
	data_type   output_register[MAX_NUM_OUTPUTS]={0ul};

	#ifdef DEBUG
		printf("\nEntering decode_cgp:\n");
	#endif

	/* only process nodes that are used */
    for (the_node = 0; the_node < num_nodes_to_process; the_node++)
    {
		/* get node function */
        function_type = phenotype[nodes_to_process[the_node]].function;

		if (function_type_is_input(function_type)) /* Input function */
		{
			in_nodes(function_type, &input_pointer,
					 the_node, nodes_to_process,
			         in,iteration,fitness_test);
		}
		else if (function_type_is_output(function_type)) /*Output function */
		{
			out_nodes(function_type,
				      &output_pointer, output_register,
					  the_node, nodes_to_process,
			          in, output, iteration);
		}
		else if (function_type_is_sm(function_type)) /* SM function */
		{
			sm_nodes(the_node, nodes_to_process,
					 in, output);
		}
		else if (function_type_is_comp(function_type)) /* comp function */
		{
			comp_nodes(function_type,
				       the_node, nodes_to_process,
					   in, output);
		}

        output[nodes_to_process[the_node]] = node_type(in,function_type);
	}

	/* write cgp outputs */
	/* this just reads them from an output register that
	   the output producing functions write to */
    for (i = 0; i < num_outputs[iteration]; i++)
	{
		 cgp_outputs[i]=output_register[i];
	}

	#ifdef DEBUG
		printf("\nLeaving decode_cgp:\n");
	#endif
}


#ifdef DATA_IS_UNSIGNED_INT
	/* counts how many bits the cgp calcualted 32-bit output
	   has in common with compressed truth table output read from read_data()
	*/
	double correctness_test(data_type data_output,
						 data_type cgp_output,
						 int iteration)
	{
		register int		i;
		int					result = 0;
		data_type			temp;

		temp=(~data_output)^cgp_output;              /* xnor the calculated output with the desired */

		for (i = 0; i < bit_width[iteration]; i++) /* examine the bit_width bits */
			result = result + getbit(temp,i);

		return (double)result;
	}
#endif

#ifdef DATA_IS_INT

	/* checks whether the output integer cgp produces
	   is the same as the desired output integer
	*/
	double correctness_test(data_type data_output,
						    data_type cgp_output,
						    int iteration)
	{
		double	result = 0.0;
		double	x;

		x = fabs(cgp_output-data_output);

		/* hits based fitness */
		#ifdef HITS_BASED_FITNESS
			if (x == 0)
				result = 1;
		#else
			result = 1.0/(1.0 + x);  /* error based fitness */
		#endif

		return result;
	}
#endif

#ifdef DATA_IS_DOUBLE

	/* checks how close the output double cgp produces
	   is to the desired output
   */
	double correctness_test(data_type data_output,
							data_type cgp_output,
							int iteration)
	{
		double	result = 0.0;
		double x;

		x = fabs(cgp_output-data_output);

		/* hits based fitness */
		#ifdef HITS_BASED_FITNESS
			if (x <= ERROR_THRESHOLD)
				result = 1.0;
		#else                         /* error based fitness */
				result = 1.0 - x;     /* an alternative is result = 1.0/(1.0 + x); */
		#endif

		return result;
	}

#endif

/* evaluate the fitness at a single test point */
double evaluate_cgp_outputs(data_type cgp_outputs[MAX_NUM_OUTPUTS],
							int test, int iteration)
{
	int i;
	double fit = 0.0, x;

	#ifdef DEBUG
		printf("\nEntering evaluate_cgp_outputs:\n");
	#endif

	for (i = 0; i < num_outputs[iteration]; i++)
	{
		x = correctness_test(data_outputs[iteration][test][i],cgp_outputs[i], iteration);
		fit = fit + x;
	}

	#ifdef DEBUG
		printf("\nLeaving evaluate_cgp_outputs:\n");
	#endif

	return fit;
}

/* determines if the fitness is in the perfection interval */
int reached_perfect_score(double fit, double perfect)
{

	if (fit >= perfect)
		return 1;
	else
		return 0;
}

/* used with Boolean data types */
void get_bit_width(int iteration)
{
	/* */
	if (num_inputs[iteration] == 2)
		bit_width[iteration] = 4;
	else if (num_inputs[iteration] == 3)
		bit_width[iteration] = 8;
	else if (num_inputs[iteration] == 4)
		bit_width[iteration] = 16;
	else
		bit_width[iteration] = 32;
}

/* determine the fitness of a phenotype at an iteration */
double fitness_at_iteration(int num_nodes_to_process,
							int nodes_to_process[MAX_NUM_NODES],
							int iteration)
{
	int			fitness_test;
	double		fit = 0.0, fitness_increment;
	data_type	cgp_outputs[MAX_NUM_OUTPUTS];

	#ifdef DEBUG
		printf("\nEntering fitness_at_iteration:\n");
	#endif

	#ifdef DATA_IS_UNSIGNED_INT
		get_bit_width(iteration);
	#endif

	for (fitness_test = 0; fitness_test < num_tests[iteration]; fitness_test++)
	{

		decode_cgp(cgp_outputs,num_nodes_to_process, nodes_to_process,
				  fitness_test, iteration);

		fitness_increment = evaluate_cgp_outputs(cgp_outputs,fitness_test, iteration);

		fit = fit + fitness_increment;

	}

	#ifdef DEBUG
		printf("\nLeaving fitness_at_iteration:\n");
	#endif

	return fit;
}

/* sort a one-dimensional arrays using insertion sort.
   Insertion sort is fast with small list sizes */
void insertionsort(int* array, int num_items)
{
    int i, j, value;

    for (i = 1; i < num_items ; i++)
    {
        value = array[i];
        j = i - 1;
        while ((j >= 0) && (array[j] > value))
        {
            array[j+1] = array[j];
            j = j - 1;
        }
        array[j+1] = value;
    }
}



/* the positions and arguments must be always be valid
   Note this function generates three arguments
   and sorts them in numerical order. Do we need this
   for all functions that change the size? ADD for
   instance only used the first argument.

   The problem that this function solves, is what to do when
   sm_arguments are no longer valid - because the phenotype
   has changed length. Simon does not repair he just discards
   the operation if it is invalid. Hmmm, this might be simpler.
*/
void make_arguments_and_position_valid(int  num_nodes_phenotype,
									   int* sm_position,
									   int  sm_arguments[MAX_NUM_ARGUMENTS],
			                           int  true_sm_arguments[MAX_NUM_ARGUMENTS])
{
	int i, true_sm_position;

	true_sm_position = *sm_position;
	/* We never said in any papers how this was handled:
	 Simon does not run an SM instruction if it isn't valid.
	 In other words it becomes a junk node.
	 */

	/* sm_positions must always be within the phenotype size */
	if (num_nodes_phenotype > 1)
		true_sm_position = true_sm_position % num_nodes_phenotype;

	/* valid sm_arguments must always be within phenotype size. Here the
	validated sm_position (x in papers) is added to the arguments and the modulus
	taken. Here they operate from the current sm position,
	is this a good idea or should they just be absolute?   */
	for (i = 0; i < MAX_NUM_ARGUMENTS; i++)
		true_sm_arguments[i] = (sm_arguments[i]+true_sm_position) % num_nodes_phenotype;

	insertionsort(true_sm_arguments,MAX_NUM_ARGUMENTS);

	*sm_position = true_sm_position;
}


/* calculate new phenotype size after SM operations
*/
int  get_size(int   num_nodes_phenotype,
			  int   true_sm_arguments[MAX_NUM_ARGUMENTS],
			  int   sm_operation)
{
	int    how_many_nodes = 0;
	int    num_nodes_to_delete;
	int    num_nodes_to_add;
	int    num_nodes_after;
	int    num_nodes_to_overwrite;

	switch (sm_operation)
	{
		case 6: /* DEL */
			num_nodes_to_delete = true_sm_arguments[1] - true_sm_arguments[0];
			how_many_nodes = num_nodes_phenotype - num_nodes_to_delete;
            /* note how_many_nodes cannot be zero as num_nodes_to_delete < num_nodes_phenotype */
		break;
		case 7:   /* Add (ADD) Add P0 new random nodes after (x) */
			num_nodes_to_add = true_sm_arguments[0];
            how_many_nodes = num_nodes_phenotype + num_nodes_to_add;
			if (how_many_nodes > MAX_NUM_NODES)
			{
				printf("\nhow_many_nodes exceeds  the max allowed (%d) nodes\n", MAX_NUM_NODES);
				printf("\nRe-run the experiments with a larger value of MAX_NUM_NODES\n");
				exit(0);
			}
		break;
		case 9:  /* OVR */
			num_nodes_after = num_nodes_phenotype - true_sm_arguments[2];
			num_nodes_to_overwrite = true_sm_arguments[1] - true_sm_arguments[0];
			if (num_nodes_after <= num_nodes_to_overwrite)
				how_many_nodes = true_sm_arguments[2] + num_nodes_to_overwrite + 1;
			else
				how_many_nodes = num_nodes_phenotype;
		break;
		case 10: /* DUP */
			num_nodes_to_add = true_sm_arguments[1] - true_sm_arguments[0];
            how_many_nodes = num_nodes_phenotype + num_nodes_to_add;
		break;
		case 11: /* CRP */
			how_many_nodes = true_sm_arguments[1] - true_sm_arguments[0];
		break;
		default: /* no change */
			how_many_nodes = num_nodes_phenotype;
	}

	if (how_many_nodes > MAX_NUM_NODES)
	{
		printf("\nIn get_size: sm_operation is %d, too many nodes\n", sm_operation);
		exit(0);
	}

	return how_many_nodes;
}



/* The DEL operation deletes the nodes between (P0 + x) and (P1 + x) */
/* returns the new phenotype and the number of nodes in the
   phenotype. */
/* Simon's version did this
Delete the nodes between (P0 + x) and (P0 + x + P1)
*/
void  DEL(int*  num_nodes_phenotype,
		  int sm_position, int   true_sm_arguments[MAX_NUM_ARGUMENTS])
{
	int	   i, j;
	int    how_many_nodes = 0;
	int    num_nodes_p;
	node   copy_of_phenotype[MAX_NUM_NODES];


	num_nodes_p = *num_nodes_phenotype;

	copy_phen(phenotype,copy_of_phenotype,num_nodes_p);

	how_many_nodes = get_size(num_nodes_p,true_sm_arguments, 6);

	/* copy over undeleted sections */
	j = 0;
	for (i = 0; i < true_sm_arguments[0]; i++)
	{
		phenotype[j] = copy_node(copy_of_phenotype[i]);
		j++;
	}
	for (i = true_sm_arguments[1]; i < num_nodes_p; i++)
	{
		phenotype[j] = copy_node(copy_of_phenotype[i]);
		j++;
	}

	*num_nodes_phenotype = how_many_nodes;
}

/* Add (ADD) Add P0 new random nodes after (x) */
void  ADD(int*  num_nodes_phenotype,
		  int sm_position, int   true_sm_arguments[MAX_NUM_ARGUMENTS])
{
	int	   i, j;
	int    how_many_nodes = 0;
	int    num_nodes_p;
	node   rand_node;
	node  copy_of_phenotype[MAX_NUM_NODES];

	num_nodes_p = *num_nodes_phenotype;

	copy_phen(phenotype,copy_of_phenotype,num_nodes_p);


	how_many_nodes = get_size(num_nodes_p,true_sm_arguments, 7);

	j = 0;
	for (i = 0; i < sm_position; i++)
	{
		phenotype[j] = copy_node(copy_of_phenotype[i]);
		j++;
	}

	/* generate random nodes to be inserted */
	for (i = 0; i < true_sm_arguments[0]; i++)
	{
		rand_node = generate_a_random_node(j);
		phenotype[j] = copy_node(rand_node);
		j++;
	}

	/* copy in remaining part of original phenotype */
	for (i = sm_position; i < num_nodes_p; i++)
	{
		phenotype[j] = copy_node(copy_of_phenotype[i]);
		j++;
	}

	*num_nodes_phenotype = how_many_nodes;
}


/* The MOV operation moves the nodes between (P0 + x) and (P1 + x)
   and inserts them after (P2 + x)
   The function returns the new phenotype and the number of nodes in the
   phenotype.
   Simon's version did this
   Move the nodes between (P0 + x) and (P0 + x + P1) and insert after
   (P0 + x + P2).
*/

void  MOV(int*  num_nodes_phenotype,
		  int sm_position, int   true_sm_arguments[MAX_NUM_ARGUMENTS])
{
	int	   i, j;
	int    num_nodes_p;
	node  copy_of_phenotype[MAX_NUM_NODES];

	/* write the number of nodes in the phenotype tp num_nodes_p */
	num_nodes_p = *num_nodes_phenotype;

    /* should this be memcopy for speed of copying? */
	copy_phen(phenotype,copy_of_phenotype,num_nodes_p);

	/* copy up to P0 + x */
	j = 0;
	for (i = 0; i < true_sm_arguments[0]; i++)
	{
		phenotype[j] = copy_node(copy_of_phenotype[i]);
		j++;
	}
	/* copy from P1 + x to P2 +x */
	for (i = true_sm_arguments[1]; i < true_sm_arguments[2]; i++)
	{
		phenotype[j] = copy_node(copy_of_phenotype[i]);
		j++;
	}
	/* copy from P0 + x to P1 + x and insert after P2 + x */
	for (i = true_sm_arguments[0]; i < true_sm_arguments[1]; i++)
	{
		phenotype[j] = copy_node(copy_of_phenotype[i]);
		j++;
	}
	/* copy from P2 + x to end of phenotype */
	for (i = true_sm_arguments[2]; i < num_nodes_p; i++)
	{
		phenotype[j] = copy_node(copy_of_phenotype[i]);
		j++;
	}

	*num_nodes_phenotype = num_nodes_p;
}

/*
   Overwrite (OVR): Copy the nodes between (P0 + x) and (P1 + x) to position
   (P2 + x), replacing existing nodes in the target position
   The function returns the new phenotype and the number of nodes in the
   phenotype.
   Simon's version did this
   Copy the nodes between (P0 + x) and (P0 + x + P1) and insert after
   (P0 + x + P2) overwriting the nodes there in the process.
*/

void  OVR(int*  num_nodes_phenotype,
		  int sm_position, int   true_sm_arguments[MAX_NUM_ARGUMENTS])
{
	int	   i, j;
	int    num_nodes_p;
	int    how_many_nodes;
	node  copy_of_phenotype[MAX_NUM_NODES];

	num_nodes_p = *num_nodes_phenotype;

	copy_phen(phenotype,copy_of_phenotype,num_nodes_p);

	how_many_nodes = get_size(num_nodes_p,true_sm_arguments, 9);

	/* copy up to P0 + x */
	j = 0;
	for (i = 0; i < true_sm_arguments[0]; i++)
	{
		phenotype[j] = copy_node(copy_of_phenotype[i]);
		j++;
	}
	/* copy from P1 + x to P2 +x */
	for (i = true_sm_arguments[1]; i < true_sm_arguments[2]; i++)
	{
		phenotype[j] = copy_node(copy_of_phenotype[i]);
		j++;
	}
	/* copy from P0 + x to P1 + x and insert after P2 + x */
	for (i = true_sm_arguments[0]; i < true_sm_arguments[1]; i++)
	{
		phenotype[j] = copy_node(copy_of_phenotype[i]);
		j++;
	}
	/* copy from P2 + x to end of phenotype */
	for (i = true_sm_arguments[2]; i < num_nodes_p; i++)
	{
		phenotype[j] = copy_node(copy_of_phenotype[i]);
		j++;
	}

	*num_nodes_phenotype = how_many_nodes;
}

/*
   DUP copies the nodes between (P0 + x) and (P1 + x) and inserts them
   after position (P2 + x)
   The function returns the new phenotype and the number of nodes in the
   phenotype.
   Simon's version did this
   Duplication (DUP) Copy the nodes between (P0 + x) and (P0 + x + P1) and insert after
   (P0 + x + P2)
*/

void  DUP(int*  num_nodes_phenotype,
		  int   sm_position, int   true_sm_arguments[MAX_NUM_ARGUMENTS])
{
	int	   i, j;
	int    num_nodes_p;
	int    how_many_nodes;
	node   copy_of_phenotype[MAX_NUM_NODES];

	num_nodes_p = *num_nodes_phenotype;

	copy_phen(phenotype,copy_of_phenotype,num_nodes_p);

	how_many_nodes = get_size(num_nodes_p,true_sm_arguments, 10);

	/* copy up to P0 + x */
	j = 0;
	for (i = 0; i < true_sm_arguments[2]; i++)
	{
		phenotype[j] = copy_node(copy_of_phenotype[i]);
		j++;
	}
	/* copy from P0 + x to P1 + x and insert after P2 + x */
	for (i = true_sm_arguments[0]; i < true_sm_arguments[1]; i++)
	{
		phenotype[j] = copy_node(copy_of_phenotype[i]);
		j++;
	}
	/* copy from P2 + x to end of phenotype */
	for (i = true_sm_arguments[2]; i < num_nodes_p; i++)
	{
		phenotype[j] = copy_node(copy_of_phenotype[i]);
		j++;
	}

	*num_nodes_phenotype = how_many_nodes;
}

/* The CRP operation standing for CROP deletes the nodes outside of (P0 + x) and (P1 + x) */
/* returns the new phenotype and the number of nodes in the
   phenotype. */
void  CRP(int*  num_nodes_phenotype,
		  int sm_position, int   true_sm_arguments[MAX_NUM_ARGUMENTS])
{
	int	   i, j;
	int    how_many_nodes = 0;
	int    num_nodes_p;
	node   copy_of_phenotype[MAX_NUM_NODES];

	num_nodes_p = *num_nodes_phenotype;

	copy_phen(phenotype,copy_of_phenotype,num_nodes_p);

	how_many_nodes = get_size(num_nodes_p,true_sm_arguments, 11);

	/* copy from P0 + x to P1 + x */
	j = 0;
	for (i = true_sm_arguments[0]; i < true_sm_arguments[1]; i++)
	{
		phenotype[j] = copy_node(copy_of_phenotype[i]);
		j++;
	}

	*num_nodes_phenotype = how_many_nodes;
}


/* Change the connections of node P0 to those of node P1 */
void  CHC(int*  num_nodes_phenotype,
		  int sm_position, int   true_sm_arguments[MAX_NUM_ARGUMENTS])
{
	int	   j;
	int    num_nodes_p;
	node   copy_of_phenotype[MAX_NUM_NODES];

	num_nodes_p = *num_nodes_phenotype;

	copy_phen(phenotype,copy_of_phenotype,num_nodes_p);

	/* change connections of P0 node to those of P1 node */
	for (j = 0; j < num_connection_genes_per_node; j++)
		phenotype[true_sm_arguments[0]].connection[j] = copy_of_phenotype[true_sm_arguments[1]].connection[j];

	/* number of nodes is unchanged */
	*num_nodes_phenotype = num_nodes_p;
}


/* Change the function of node P0 to that of node P1 */
void  CHF(int*  num_nodes_phenotype,
		  int sm_position, int   true_sm_arguments[MAX_NUM_ARGUMENTS])
{
	int    num_nodes_p;
	node   copy_of_phenotype[MAX_NUM_NODES];

	num_nodes_p = *num_nodes_phenotype;

	copy_phen(phenotype,copy_of_phenotype,num_nodes_p);

	/* change function of P0 node to that of P1 node */
	phenotype[true_sm_arguments[0]].function = copy_of_phenotype[true_sm_arguments[1]].function;

	/* number of nodes is unchanged */
	*num_nodes_phenotype = num_nodes_p;
}


/* Change the arguments of node P0 to those of node P1 */
void  CHA(int*  num_nodes_phenotype,
		  int sm_position, int   true_sm_arguments[MAX_NUM_ARGUMENTS])
{
	int	   j;
	int    num_nodes_p;
	node   copy_of_phenotype[MAX_NUM_NODES];

	num_nodes_p = *num_nodes_phenotype;

	copy_phen(phenotype,copy_of_phenotype,num_nodes_p);

	/* change arguments of P0 node to those of P1 node */
	for (j = 0; j < MAX_NUM_ARGUMENTS; j++)
		phenotype[true_sm_arguments[0]].argument[j] = copy_of_phenotype[true_sm_arguments[1]].argument[j];

	/* number of nodes is unchanged */
	*num_nodes_phenotype = num_nodes_p;
}

/* from the ToDo list and which item on the list (which_sm)
   get the position, arguments and function
*/
void get_sm_arguments(ToDo  ToDoList[MAX_TODO],
                      int   which_sm,
					  int*  sm_position,
					  int   sm_arguments[MAX_NUM_ARGUMENTS],
					  int*  sm_operation)
{
	int j;

	*sm_position = ToDoList[which_sm].position;
	for (j = 0; j < MAX_NUM_ARGUMENTS; j++)
		sm_arguments[j] = ToDoList[which_sm].argument[j];
	*sm_operation = ToDoList[which_sm].function;
}

/* Apply the operations on ToDo list to the phenotype
   to generate the new phenotype */
void generate_new_phenotype(int*   num_nodes_phenotype,
                            ToDo   ToDoList[MAX_TODO],
						    int    num_sm_nodes_to_process)
{
	int   i;
	int   sm_operation, sm_position;
	int   sm_arguments[MAX_NUM_ARGUMENTS];
	int   true_sm_arguments[MAX_NUM_ARGUMENTS]={0};
	int   num_nodes_p;

	num_nodes_p = *num_nodes_phenotype;

	/* get sm parameters */
	for (i = 0; i < num_sm_nodes_to_process; i++)
	{
		/* extract the details of the sm operations from the
		   ToDoList */
		get_sm_arguments(ToDoList,i,
			             &sm_position,
						 sm_arguments,
						 &sm_operation);

		/* make sure the sm_position (x in papers) fits in phenotype
		   and the arguments also fit and are sorted
		   in numerically ascending order
		*/
		make_arguments_and_position_valid(num_nodes_p,
									      &sm_position,
									      sm_arguments,
			                              true_sm_arguments);
		switch (sm_operation)
		{
			case 6: /* DEL-SM-6 */
			DEL(&num_nodes_p,sm_position,true_sm_arguments);
			break;
			case 7: /* ADD_SM-7 */
			ADD(&num_nodes_p,sm_position,true_sm_arguments);
			break;
			case 8: /* MOV-SM-8 */
			MOV(&num_nodes_p,sm_position,true_sm_arguments);
			break;
			case 9: /* OVR-SM-9 */
			OVR(&num_nodes_p,sm_position,true_sm_arguments);
			break;
			case 10:/* DUP-SM-10*/
			DUP(&num_nodes_p,sm_position,true_sm_arguments);
			break;
			case 11:/* CRP-SM-11*/
			CRP(&num_nodes_p,sm_position,true_sm_arguments);
			break;
			case 12:/* CHC-SM-12*/
			CHC(&num_nodes_p,sm_position,true_sm_arguments);
			break;
			case 13:/* CHF-SM-13*/
			CHF(&num_nodes_p,sm_position,true_sm_arguments);
			break;
			case 14:/* CHA-SM-14*/
			CHA(&num_nodes_p,sm_position,true_sm_arguments);
			break;
		}
	}

	*num_nodes_phenotype = num_nodes_p;
}



void generate_new_phenotype_verbose(int*   num_nodes_phenotype,
									ToDo   ToDoList[MAX_TODO],
						            int    num_sm_nodes_to_process,
									char   name[MAX_NUM_LETTERS],
									int    iteration)
{
	int   i, j;
	int   sm_operation, sm_position;
	int   sm_arguments[MAX_NUM_ARGUMENTS];
	int   true_sm_arguments[MAX_NUM_ARGUMENTS]={0};
	int   num_nodes_p;
	FILE* fp;

	num_nodes_p = *num_nodes_phenotype;

	/* get sm parameters */
    fp = fopen(name, "a");
	fprintf(fp, "\nThe sm nodes to process are: \n");
	fclose(fp);
	for (i = 0; i < num_sm_nodes_to_process; i++)
	{
		/* extract the details of the sm operations from the
		   ToDoList */
		get_sm_arguments(ToDoList,i,
			             &sm_position,
						 sm_arguments,
						 &sm_operation);

		fp = fopen(name, "a");
		fprintf(fp, "[%s](..:",node_types[sm_operation]);
		for (j = 0; j < MAX_NUM_ARGUMENTS; j++)
			fprintf(fp, "%d,",sm_arguments[j]);
		fprintf(fp, "):%d", sm_position);
		fclose(fp);

		/* make sure the sm_position (x in papers) fits in phenotype
		   and the arguments also fit and are sorted
		   in numerically ascending order
		*/
		make_arguments_and_position_valid(num_nodes_p,
									      &sm_position,
									      sm_arguments,
			                              true_sm_arguments);

		fp = fopen(name, "a");
		fprintf(fp, "\nAfter make_arguments_and_position_valid\n");
		fprintf(fp, "[%s](..:",node_types[sm_operation]);
		for (j = 0; j < MAX_NUM_ARGUMENTS; j++)
			fprintf(fp, "%d,",true_sm_arguments[j]);
		fprintf(fp, "):%d", sm_position);
		fclose(fp);

		switch (sm_operation)
		{
			case 6: /* DEL-SM-6 */
			DEL(&num_nodes_p,sm_position,true_sm_arguments);
			break;
			case 7: /* ADD_SM-7 */
			ADD(&num_nodes_p,sm_position,true_sm_arguments);
			break;
			case 8: /* MOV-SM-8 */
			MOV(&num_nodes_p,sm_position,true_sm_arguments);
			break;
			case 9: /* OVR-SM-9 */
			OVR(&num_nodes_p,sm_position,true_sm_arguments);
			break;
			case 10:/* DUP-SM-10*/
			DUP(&num_nodes_p,sm_position,true_sm_arguments);
			break;
			case 11:/* CRP-SM-11*/
			CRP(&num_nodes_p,sm_position,true_sm_arguments);
			break;
			case 12:/* CHC-SM-12*/
			CHC(&num_nodes_p,sm_position,true_sm_arguments);
			break;
			case 13:/* CHF-SM-13*/
			CHF(&num_nodes_p,sm_position,true_sm_arguments);
			break;
			case 14:/* CHA-SM-14*/
			CHA(&num_nodes_p,sm_position,true_sm_arguments);
			break;
		}
	}

	fp = fopen(name,"a");
	fprintf(fp,"\nnew phenotype is\n");
	fclose(fp);

	fprint_a_phenotype(phenotype,num_nodes_p, name, 1);

	fprint_active_phenes(num_nodes_p, name, iteration);

	*num_nodes_phenotype = num_nodes_p;
}

/* Calculates the addresses of nodes used
   and stores them in nodes_to_process.
   Returns the number of nodes used.
   Also determines the ToDoList and how many sm_nodes need
   to be processed in it (it has a maximum of ToDoLength).
   Iteration is needed as the nodes that need to be
   processed depend on how many outputs are required
*/
int get_nodes_to_process(int   num_nodes_phenotype,
						 ToDo  ToDoList[MAX_TODO],
						 int*  num_sm_nodes_to_process,
						 int   nodes_to_process[MAX_NUM_NODES],
						 int   iteration)
{
	int		i, j;
	int		num_nodes_to_process = 0;
    int		node_used[MAX_NUM_NODES];
	int     num_sm_nodes_used = 0;

	get_node_used_in_phenotype(num_nodes_phenotype,
				               node_used);

	/* find number of used nodes */
	for (i = 0; i < num_nodes_phenotype; i++)
	{
		if (node_used[i])
		{
			nodes_to_process[num_nodes_to_process] = i;
			num_nodes_to_process++;
			if (node_is_sm(i)) /* it is used and is sm */
			{
				if (num_sm_nodes_used < ToDoLength) /* it is within budget */
				{
					ToDoList[num_sm_nodes_used].position = i;
					for (j = 0; j < MAX_NUM_ARGUMENTS; j++)
						ToDoList[num_sm_nodes_used].argument[j] = phenotype[i].argument[j];
					ToDoList[num_sm_nodes_used].function = phenotype[i].function;
					num_sm_nodes_used++;
				}
			}
		}
	}

	*num_sm_nodes_to_process= num_sm_nodes_used;

	return num_nodes_to_process;
}

/* generates the new phenotypes from iteration 1 onwards and
   accumulates the fitness
*/
double iterate_phenotype_and_accumulate_fitness(node* genotype, int* last_perfect_iteration)
{
	int     num_nodes_phenotype;
	ToDo    ToDoList[MAX_TODO];
	int     num_sm_nodes_to_process;
	int		num_nodes_to_process = 0;
	int     nodes_to_process[MAX_NUM_NODES];
	int     iteration = 0;
	int     quit;
	double	fit = 0.0;
	double  fit_at_it;

	num_nodes_phenotype = num_nodes;
	copy_chrom_to_phenotype(genotype,num_nodes_phenotype);
	do
	{
		quit = 0;
		num_nodes_to_process = get_nodes_to_process(num_nodes_phenotype,
                                                    ToDoList,
                                                    &num_sm_nodes_to_process,
                                                    nodes_to_process, iteration);

		fit_at_it = fitness_at_iteration(num_nodes_to_process,nodes_to_process,iteration);
		fit = fit + fit_at_it;

        /* if fit_at_it is not perfect then stop iterating */
		if (fit_at_it < perfect_at_it[iteration])
		{
            quit = 1;
            *last_perfect_iteration = iteration-1;
		}
        else /* fit_at_it >= perfect_at_it - it can be greater for symbolic regression */
        {
            if  (iteration != num_iterations)
                if (num_sm_nodes_to_process  != 0)
                    generate_new_phenotype(&num_nodes_phenotype,ToDoList,num_sm_nodes_to_process);
            *last_perfect_iteration = iteration;
            iteration++;
        }
	}
	while ((iteration <= num_iterations) && (quit == 0));

	return fit;
}


double iterate_phenotype_and_accumulate_fitness_with_comments(node* genotype)
{
	int     num_nodes_phenotype;
	ToDo    ToDoList[MAX_TODO];
	int     num_sm_nodes_to_process;
	int		num_nodes_to_process = 0;
	int     nodes_to_process[MAX_NUM_NODES];
	int     iteration = 0;
	int     quit;
	double	fit = 0.0;
	double  fit_at_it;

	num_nodes_phenotype = num_nodes;

	/* make the initial phenotype a copy of the genotype */
	copy_chrom_to_phenotype(genotype,num_nodes_phenotype);

	do
	{
		quit = 0;
		/* analyze phenotype and find out which nodes we
		   need to process and also while we are at it
		   determine the ToDoList */

		num_nodes_to_process = get_nodes_to_process(num_nodes_phenotype,
			                                        ToDoList,
			                                        &num_sm_nodes_to_process,
			                                        nodes_to_process,
			                                        iteration);

		/* calculate the fitness of this phenotype */
		fit_at_it = fitness_at_iteration(num_nodes_to_process,
			                             nodes_to_process,
										 iteration);
		fit = fit + fit_at_it;
		iteration++;

		if (fit_at_it != perfect_at_it[iteration])
			quit = 1;
		else /* got perfect fitness at iteration */
			quit = 0;

		if  ( (iteration != (num_iterations+1)) && (quit == 0) )
		{
			/* generate the new phenotype */
			if (num_sm_nodes_to_process  != 0)
			{

				generate_new_phenotype(&num_nodes_phenotype,
			                           ToDoList,
								       num_sm_nodes_to_process);
			}
		}

	} /* end of iteration loop */
	while ((iteration <= num_iterations) && (quit == 0));

	return fit;
}


/* this is the EA fitness function
*/
double fitness(node* genotype, int* last_perfect_iteration)
{
	double		fit;

	#ifdef DEBUG
		printf("\nEntering fitness:\n");
	#endif

	/* iterate phenotype */
	fit = iterate_phenotype_and_accumulate_fitness(genotype, last_perfect_iteration);

	#ifdef DEBUG
		printf("\nLeaving fitness:\n");
	#endif

	return fit;
}


/* returns a random valid connection allele that
   obeys the constraints imposed by levels_back.
*/
int get_connection_allele(int which_node)
{
	int gene;

	if (which_node > 0)
	{
		if (levels_back > which_node)
			gene = newrand(which_node)+1;
		else
			gene = newrand(levels_back)+1;
	}
	else
		gene = 0;

	return gene;
}


/* returns a random valid argument allele
   shouldn't one provide how many nodes? */
int get_argument_allele(void)
{
	int gene;

	gene = newrand(num_nodes);

	return gene;
}

/* returns a random valid function allele */
int get_function_allele(void)
{
	return allowed_functions[newrand(num_functions)][0];
}


/* returns a random allowed output function allele */
int get_output_function_allele(void)
{
	int which_output_function;

	if (num_output_functions == 1)
		return allowed_output_functions[0];
	else
	{
		which_output_function = newrand(num_output_functions);
		return allowed_output_functions[which_output_function];
	}
}

/* creates a randomly generated node */
node generate_a_random_node(int position)
{
	int		i;
	node	rand_node;

	/* get random connection genes */
	for (i = 0; i < num_connection_genes_per_node; i++)
		rand_node.connection[i] = get_connection_allele(position);

	/* get random argument genes */
	for (i = 0; i < MAX_NUM_ARGUMENTS; i++)
		rand_node.argument[i] = get_argument_allele();

	rand_node.function = get_function_allele();

	return rand_node;

}


/* generates a random chromosome. Used by initialise */
void generate_a_random_chromosome(node* chromosome)
{
	int the_node;

	for (the_node = 0; the_node < num_nodes; the_node++)
		chromosome[the_node] = generate_a_random_node(the_node);

}


/* creates initial population of chromosomes
   either having been generated from a single
   chromosome from a file or by generating
   an entire random population
*/
void initialise(node** chromosomes)
{
    int  i;

    if (run_from_chrom)
		read_from_chrom(chromosomes);
    else  /* generate random population */
		for (i = 0; i < population_size; i++)
			generate_a_random_chromosome(chromosomes[i]);
}

/* calculate best population fitness and the best chromosome */
double  get_best_chromosome(node** chromosomes,
                            node*  best_chromosome,
						    double previous_best_fitness,
						    int    gen,
						    int*   last_perfect_iteration)
{
	int		i, last_p_it;
    double	fitness_max, fit;
    int		best_member;

	fitness_max = -1.0;
    best_member = 0;

	#ifdef DEBUG
			printf("\nEntering get_best_chromosome:\n");
	#endif

    for (i = 0; i < population_size; i++)
    {
        last_p_it = 0;

		if ((i == population_size -1) && (gen > 1))
			fit = previous_best_fitness;
		else
		{
			fit = fitness(chromosomes[i], &last_p_it);
		}

		if (fit > fitness_max)
		{
			fitness_max = fit;
			best_member = i;
			*last_perfect_iteration = last_p_it;
		}
		/* break out of this as soon as we get a perfect score */
		if ( reached_perfect_score(fit, perfect))
			break;
	}

	/* store the best chromosome */
    copy_pop_to_chrom(chromosomes, best_chromosome, best_member);

	#ifdef DEBUG
			printf("\nLeaving get_best_chromosome:\n");
	#endif

	return fitness_max;
}


/* get a random connection gene */
int get_valid_connection_gene(int num_connection_genes)
{
	int new_gene;

	new_gene = newrand(num_connection_genes);

	return new_gene;
}

/* get a random function gene */
int get_valid_function_gene(int num_functions)
{
	int new_gene;

	new_gene = newrand(num_functions);

	return new_gene;
}

/* get a random argument gene */
int get_valid_argument_gene(int num_argument_genes)
{
	int new_gene;

	new_gene = newrand(num_argument_genes);

	return new_gene;
}

/* mutate a connection gene in a chromosome */
void mutate_a_connection_gene(node*  chromosome)
{
	int which_connection, which_connection_gene;
	int which_node;

	which_connection_gene = get_valid_connection_gene(num_connection_genes);

	which_node = which_connection_gene/num_connection_genes_per_node;

	which_connection = which_connection_gene % num_connection_genes_per_node;

	chromosome[which_node].connection[which_connection]=get_connection_allele(which_node);

}

/* mutate a function gene in a chromosome */
void mutate_a_function_gene(node*  chromosome)
{
	int which_node;

	which_node = get_valid_function_gene(num_nodes);

	chromosome[which_node].function = get_function_allele();
}

/* mutate an argument gene in a chromosome */
void mutate_an_argument_gene(node*  chromosome)
{
	int which_argument, which_argument_gene;
	int which_node;

	which_argument_gene = get_valid_argument_gene(num_argument_genes);

	which_node = which_argument_gene/num_nodes;

	which_argument = which_argument_gene % MAX_NUM_ARGUMENTS;

	chromosome[which_node].argument[which_argument] = get_argument_allele();

}

/* carry out num_mutations mutations on the chromosome */
void mutate_a_chromosome(node* chromosome,
						 int num_connection_mutations,
						 int num_function_mutations,
						 int num_argument_mutations)
{
	int i;

	for (i = 0; i < num_connection_mutations; i++)
		mutate_a_connection_gene(chromosome);

	for (i = 0; i < num_function_mutations; i++)
		mutate_a_function_gene(chromosome);

	for (i = 0; i < num_argument_mutations; i++)
		mutate_an_argument_gene(chromosome);

}

/* (1+lambda evolutionary strategy where lamda = population size -1 */
void generate_new_pop_es(node** chromosomes,
						 node*  best_chromosome,
						 int num_connection_mutations,
						 int num_function_mutations,
						 int num_argument_mutations)
{
	int i;


     /* copy best_chromosome into last member of chromosome array */
    copy_chrom_to_pop(chromosomes, best_chromosome, population_size-1);


    /* generate new population by mutating all but last */
    for (i = 0; i < population_size-1; i++)
    {
       /* copy best chromosome to population member i */
       copy_chrom_to_pop(chromosomes, best_chromosome, i);

	   /* mutate the chromosome */
	   mutate_a_chromosome(chromosomes[i],
		                   num_connection_mutations,
			               num_function_mutations,
			               num_argument_mutations);
    }
}

/* allocate space for a 1d datatype array with size dimension1 */
data_type* create_1d_datatype_space(int dimension1)
{
	data_type* array1d = NULL;

	array1d = (data_type*) calloc(dimension1, sizeof(data_type));

    if (array1d == NULL)
    {
        printf("ERROR\n");
	    printf("In function create_1d_datatype_space\n");
		printf("Not enough memory for a 1d datatype array of length %d\n", dimension1);
        exit(0);
    }
	return array1d;
}

/* allocate space for 2 dimensional datatype array  */
data_type** create_2d_datatype_space(int dimension1, int dimension2)
{
	int			i;
	data_type**	array2d = NULL;

	/* create space for pointers to pointers */
	array2d = (data_type**) calloc(dimension2, sizeof(data_type*));
	if (array2d == NULL)
	{
       printf("ERROR\n");
	   printf("In function create_2d_datatype_space\n");
	   printf("Can not allocate space for %d many datatype pointers\n",dimension2);
       exit(0);
	}

	/* create array of pointers  */
	for (i = 0; i < dimension2; i++)
	{
	  array2d[i] = create_1d_datatype_space(dimension1);
      if (array2d[i] == NULL)
      {
		  printf("ERROR\n");
		  printf("In function create_2d_datatype_space\n");
		  printf("Not enough memory for datatype pointer arrays of length %d\n",dimension1);
		  exit(0);
      }
   }

	return array2d;
}

/* allocate space for 3 dimensional data_type array  */
data_type*** create_3d_datatype_space(int dimension1, int dimension2, int dimension3)
{
	int			i;
	data_type***	array3d = NULL;

	/* create space for pointers to int pointers */
	array3d = (data_type*** ) calloc(dimension3, sizeof(data_type**));
	if (array3d == NULL)
	{
       printf("ERROR\n");
	   printf("In function create_3d_datatype_space\n");
	   printf("Can not allocate space for %d many datatype pointers\n",dimension3);
       exit(0);
	}

	/* create array of pointers to data_type**  */
	for (i = 0;i < dimension3;i++)
	{
	  array3d[i] = create_2d_datatype_space(dimension1, dimension2);
      if (array3d[i] == NULL)
      {
		  printf("ERROR\n");
		  printf("In function create_3d_datatype_space\n");
		  printf("Not enough memory for data_type pointer arrays of length %d\n",dimension2);
		  exit(0);
      }
   }

	return array3d;
}

/* release 2d data_type memory */
void free_datatype_array2d(int dimension2, data_type** array2d)
{
	int i;

	/* free 1d array of pointers  */
	for (i = 0; i < dimension2; i++)
		free(array2d[i]);

	free(array2d);
}

/* release 3d data_type memory */
void free_datatype_array3d(int dimension2, int dimension3, data_type*** array3d)
{
	int i;

	/* free 2darray pointers  */
	for (i = 0; i < dimension3; i++)
		free_datatype_array2d(dimension2, array3d[i]);

	free(array3d);
}



/* allocate space for a 1d array of nodes with size items */
node* create_1d_array_of_nodes(int size)
{
	node* array1d = NULL;

	array1d = (node*) calloc(size,sizeof(node));

    if (array1d == NULL)
    {
		printf("ERROR.Not enough memory for a 1d array of nodes of length %d\n", size);
        exit(0);
    }

	return array1d;
}

/* allocate space for 2 dimensional array of  node structures
   e.g. a population of chromosomes */
node** create_2d_array_of_nodes(int num_horizontals, int num_verticals)
{
	int i;
	node **array2d = NULL;

	/* create space for pointers to int pointers */
	array2d = (node** ) calloc(num_verticals, sizeof(node*));
	if (array2d == NULL)
	{
       printf("ERROR.Can not allocate space for %d many node pointers\n",num_verticals);
       exit(0);
	}

	/* create array of pointers to ints  */
	for (i = 0; i < num_verticals; i++)
	{
	  array2d[i] = create_1d_array_of_nodes(num_horizontals);
      if (array2d[i] == NULL)
      {
         printf("ERROR.Not enough memory for node pointer arrays of length %d\n",num_horizontals);
         exit(0);
      }
   }

	return array2d;
}


/* release memory */
/* if this is for chromosomes then num_verticals is population_size */
void free_array2d_of_nodes(int num_verticals, node** array2d)
{
	int i;

	/* free 1darray of pointers  */
	for (i = 0; i < num_verticals; i++)
		free(array2d[i]);

	free(array2d);

}


void write_generation_to_screen(int gen_index)
{
	if (gen_index % report_interval==0)
		printf("\nGENERATION is %d",gen_index);
}

/* writes generation and fitness of best in population */
void write_progress_info_to_screen(int generation, double fit, int last_perfect_iteration)
{
	printf("\nGENERATION is %d Best fitness is now %8.5lf. Last perfect iteration is %d",generation,fit, last_perfect_iteration);
}


/* writes out chromosome to file defined by string prog */
void write_progress_info_to_file(char prog[MAX_NUM_LETTERS],
								 int gen, double best_fit,
						         node* best_chromosome)
{
	FILE* fp;

    fp=fopen(prog,"a");
    fprintf(fp,"\nGENERATION is %u     Best fitness is now %10.6lf",gen,best_fit);
    fprintf(fp,"\nThe chromosome is\n");
    fclose(fp);

	fprint_all_phenotypes(best_chromosome,prog,0, 0, 0);

}

/* checks to see if the best fitness in the population has improved.
   writes the generation, the new best fitness and the improved chromosome
   to the progress report (if progress report is set to 1)
*/
void check_if_improvement(double best_fit, double* previous_best_fit, int* best_gen, int gen,
						  char prog[MAX_NUM_LETTERS],
						  node* best_chromosome, int last_perfect_iteration)
{
	if (best_fit > *previous_best_fit) /* we have an improvement */
	{
		if (progress_report)
		{
			write_progress_info_to_screen(gen,best_fit, last_perfect_iteration);
			write_progress_info_to_file(prog, gen, best_fit, best_chromosome);
		}
		*best_gen = gen;				/* update the best generation */
		*previous_best_fit = best_fit;	/* update previous best fitness */
	}
}


/* report on results of evolutionary run in smcgp.txt */
void write_result_of_EA_to_file(int run, int bgen, double best_fit,
								node* best_chromosome)
{
	FILE* fp;

	fp=fopen("smcgp.txt","a");
	fprintf(fp,"Run %d and gen %d achieved fitness %10.6lf\n",run,bgen,best_fit);
	fprintf(fp,"Here is the chromosome\n");
	fclose(fp);

    fprint_all_phenotypes(best_chromosome, "smcgp.txt", create_dot_files, run, create_io_files);
}



/* Do a run of the EA */
double  EA(int *gen_of_best, int run,
		char prog[MAX_NUM_LETTERS])
{
	int		gen, best_gen;
	int     last_perfect_iteration = 0;
	node**	chromosomes;
	node*	best_chromosome;
	double	best_fit = 0.0, previous_best_fit = -1.0;


	#ifdef DEBUG
		printf("\nEntering EA:\n");
	#endif

	chromosomes = create_2d_array_of_nodes(num_nodes, population_size);
	best_chromosome = create_1d_array_of_nodes(num_nodes);

	initialise(chromosomes);

	for (gen = 1 ; gen <= num_generations; gen++)
	{

		write_generation_to_screen(gen);
		best_fit = get_best_chromosome(chromosomes,
									   best_chromosome,
									   previous_best_fit, gen, &last_perfect_iteration);


		check_if_improvement(best_fit, &previous_best_fit, &best_gen, gen, prog,
							 best_chromosome, last_perfect_iteration);

		/* jump out of run if maximum fitness acheived */
		if (reached_perfect_score(best_fit, perfect))
			break;
		else /* create a new population */
			generate_new_pop_es(chromosomes,best_chromosome,
			                    num_connection_mutations,
			                    num_function_mutations,
			                    num_argument_mutations);
	}

	write_result_of_EA_to_file(run, best_gen, best_fit, best_chromosome);

	*gen_of_best = best_gen;

	/* write the raw best chromosome to smcgp.chr */
	fprint_a_raw_chromosome(best_chromosome, num_nodes,"smcgp.chr",0);

	free_array2d_of_nodes(population_size, chromosomes);

	free(best_chromosome);

	#ifdef DEBUG
		printf("\nLeaving EA:\n");
	#endif

	return best_fit;
}

void setup_report_files(int run, char prog[MAX_NUM_LETTERS])
{
	char runstring[MAX_NUM_LETTERS];
	FILE* fp;

	sprintf(runstring,"%d",run); /* store run as characters */
	if (progress_report > 0)
	{
		strcpy(prog,"smcgp");
		strcat(prog,runstring);
		strcat(prog,".prg"); /* create .prg file name */
		fp=fopen(prog,"w");  /* create empty .prg file */
		fclose(fp);
	}
}



void report_final_results(double av_fitness, double st_dev,
						  double best_of_best_fit,
						  double worst_of_best_fit,
						  double av_best_gen,
						  double av_best_gen_perfect,
						  int    num_perfect,
						  double av_num_evals_perfect)
{
	FILE* best;
	best=fopen("smcgp.txt","a");
	fprintf(best,"\naverage fitness  %6.4lf\n",av_fitness);
	fprintf(best,"\nstd dev          %6.4lf\n\n",st_dev);
	fprintf(best,"\nThe average best generation is  %6.4lf\n",av_best_gen);
	fprintf(best,"\nThe best fitnes of all runs  is  %6.4lf\n",best_of_best_fit);
	fprintf(best,"\nThe worst fitness of all runs is  %6.4lf\n",worst_of_best_fit);
	if (num_perfect > 0)
	{
		fprintf(best,"\nNumber of perfect solutions is %d\n",num_perfect);
		fprintf(best,"\nOf perfect solutions found, the average number of generations is  %6.4lf\n",av_best_gen_perfect);
		fprintf(best,"\nOf perfect solutions found, the average number of evals is  %6.4lf\n",av_num_evals_perfect);
	}
	fclose(best);
}

double average(int num_items, double items[MAX_NUM_RUNS])
{
	int		i;
	double	av;

	av = 0.0;
	for (i = 0; i < num_items; i++)
		av = av + items[i];

	av = av/num_items;

	return av;
}

double get_standard_deviation(int num_items, double average,
							  double items[MAX_NUM_RUNS])
{
	int i;
	double temp, st_dev;

	st_dev = 0.0;
	for (i = 0; i < num_items;i++)
	{
		temp=(items[i]-average);
		temp=temp*temp;
		st_dev=st_dev+temp;
	}

	st_dev=st_dev/(double)num_items;
	st_dev=sqrt(st_dev);

	return st_dev;
}

/* do mutiple runs of EA and write out results */
void run_EA(int num_runs_total)
{
	int		best_gen,run, num_perfect = 0;
	double  av_num_evals_perfect;
	double	generations_for_best[MAX_NUM_RUNS];
	double  generations_for_perfect[MAX_NUM_RUNS];
	double  evals_for_perfect[MAX_NUM_RUNS];
	char	prog[MAX_NUM_LETTERS];
	double  worst_of_best_fit, best_of_best_fit = 0.0;
	double  fitness_final;
	double	st_dev, av_fitness;
	double  av_best_gen, av_best_gen_perfect;
	double	fitnesses[MAX_NUM_RUNS];

	#ifdef DEBUG
		printf("\nEntering run_EA:\n");
	#endif

	worst_of_best_fit= perfect + num_nodes;
	for (run = 0; run < num_runs_total; run++)
	{
		printf("\n\nRUN %d\n",run);

		setup_report_files(run, prog);

		fitness_final =  EA(&best_gen, run, prog);

		fitnesses[run] = fitness_final;
		generations_for_best[run] = best_gen;

        if (fitness_final < worst_of_best_fit)
			worst_of_best_fit=fitness_final;

		if (fitness_final > best_of_best_fit)
			best_of_best_fit=fitness_final;

		if (reached_perfect_score(fitness_final, perfect))
		{
			generations_for_perfect[num_perfect] = best_gen;
			evals_for_perfect[num_perfect] = 4*best_gen+1;
			num_perfect++;
		}
	}

	free_datatype_array3d(max_num_tests, num_iterations+1, data_inputs);
	free_datatype_array3d(max_num_tests, num_iterations+1, data_outputs);

	av_best_gen = average(num_runs_total, generations_for_best);
	av_num_evals_perfect = average(num_perfect, evals_for_perfect);
	av_best_gen_perfect = average(num_perfect, generations_for_perfect);
	av_fitness = average(num_runs_total, fitnesses);
	st_dev = get_standard_deviation(num_runs_total, av_fitness, fitnesses);
	report_final_results(av_fitness, st_dev,
						 best_of_best_fit, worst_of_best_fit,
						 av_best_gen, av_best_gen_perfect,
						 num_perfect, av_num_evals_perfect);

	#ifdef DEBUG
		printf("\nLeaving run_EA:\n");
	#endif
}

