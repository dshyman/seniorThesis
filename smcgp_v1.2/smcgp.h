/* smcgp.h
   Julian F. Miller (c), 2013
   version 1.2 (Aug 2013)
   Dept. of Electronics, University of York, UK
   smcgp is a generalization of classic CGP (assuming
   one row of nodes). If all self-modifying are switched off
   and num_iterations is set to zero in the .par file, then it becomes equivalent to CGP.
   However, addressing is relative rather than absolute (as in
   classic CGP) and inputs and outputs are acquired through
   functions (i.e. terminals have been abandoned).
*/

#include <stdio.h>

#ifndef __SMCGP_H__

#define __SMCGP_H__

/* invoke alternative data types by uncommenting the one you want
   and commenting out the others
*/


/* This is for compressed truth tables: 32-bit data */
/*
#define DATA_IS_UNSIGNED_INT
*/

/* At present this assumes input and output data is integer but cgp
   processes it as double */
/*
#define DATA_IS_INT
*/

/* comment this in if you want to work with real-values */
/* */
#define DATA_IS_DOUBLE
/* */

/* this if statement defines the type of data that
   is processed by cgp genotypes
   if you have a user defined data type
   that is not on this list then you need to add a case
   to the if statement
   Note this will mean that you will have to define a new function set
   this is defined in the c program file node_function.c
*/

/* this sets an upper bound on the allowed number of inputs per node */
/* defines upper bound on the number of arguments that SM functions have */
#define MAX_NUM_ARGUMENTS       3

#ifdef DATA_IS_UNSIGNED_INT
		typedef unsigned long	data_type;
		#define MAX_ARITY               3 /* includes boolean functions and MUX */
#endif

#ifdef DATA_IS_INT
		typedef double			data_type;
		#define MAX_ARITY               2 /* math functions */

#endif

#ifdef DATA_IS_DOUBLE
		typedef double			data_type;
		#define MAX_ARITY               2 /* math functions */

#endif

/* This statement defines what is meant by a node.
   In standard CGP nodes do not have arguments.
   In SMCGP, arguments are used in self-modification
   functions. You may wish to adapt the definition of node
   to be some other data structure that suits your application
*/
typedef struct
{
	int connection[MAX_ARITY];
	int argument[MAX_NUM_ARGUMENTS];
	int function;

} node;

/* This makes ToDo contain all the info needed to carry out
   sm operations */
typedef struct
{
	int position;
	int argument[MAX_NUM_ARGUMENTS];
	int function;
} ToDo;

#define PI										3.1415926535897932
#define MAX_NUM_NODES							2000
#define MAX_NUM_OUTPUTS							30
#define MAX_NUM_CONNECTION_GENES_PER_NODE		3
#define MAX_NUM_CONNECTIONS_PLUS_FUNCTION       MAX_NUM_CONNECTION_GENES_PER_NODE + 1
#define MAX_NUM_RUNS							200
#define ERROR_THRESHOLD							0.01
#define PERFECTION_THRESHOLD					0.000001
#define MAX_NUM_INPUT_FUNCTIONS                 3
#define MAX_NUM_OUTPUT_FUNCTIONS                3
#define MAX_NUM_IO_FUNCTIONS					MAX_NUM_INPUT_FUNCTIONS+MAX_NUM_OUTPUT_FUNCTIONS
#define MAX_NUM_SM_FUNCTIONS					9  /* this needs to agree with number of sm function in .par file */
#define MAX_NUM_COMP_FUNCTIONS					21 /* this needs to agree with number of comp function in .par file */
#define MAX_NUM_ITERATIONS						10
#define MAX_TODO                                10
#define CRASH_GEN                               10086


/* setting this to 1 will compile in entry and exit messages
into many functions */
/*
#define DEBUG
*/

/*
#define HITS_BASED_FITNESS
*/

#define MAX_NUM_CHROMOSOMES			1000  /* max population size */
#define MAX_NUM_LETTERS				100
#define MAX_NUM_FUNCTIONS			MAX_NUM_IO_FUNCTIONS+MAX_NUM_SM_FUNCTIONS+MAX_NUM_COMP_FUNCTIONS
/* 2^32-1 */
#define MAXNUM						4294967295


/* these are all read from the .par file and never change after that */
int				population_size;
double			per_cent_connection_mutate;
double			per_cent_function_mutate;
double			per_cent_argument_mutate;
int				num_generations;
int				num_runs_total;
int				num_nodes;
int				levels_back;
int				progress_report;
int				report_interval;
unsigned		global_seed;
int				save_best_chrom;
int				run_from_chrom;
char			targetfiles[MAX_NUM_LETTERS]; /* defines file that lists problems to be solved */
int             ToDoLength;			              /* allowed length of ToDo list */
int             num_iterations;               /* if < 0 then program allows evolution to decide */
int             create_dot_files; /* decides if .dot files should be created after each run */
int             create_io_files;  /* decides if .dat file are to be generated from each run */
int             display_junk;  /* used to choose whether to show unconnected nodes in dot graphs */

/*  global variables calculated or obtained in get_parameters.
    They never change after that
*/
int				num_functions;
int             num_output_functions;
int             num_input_functions;
int             num_connection_genes_per_node;
int             num_connection_genes;
int				num_argument_genes;
int				num_genes_per_node;
int				number[MAX_NUM_FUNCTIONS];
char			node_types[MAX_NUM_FUNCTIONS][MAX_NUM_LETTERS];
char			target_file_names[MAX_NUM_ITERATIONS][MAX_NUM_LETTERS]; /* stores the names of the target files */
int				num_connection_mutations;
int				num_function_mutations;
int				num_argument_mutations;

/* this stores the node function address in
   allowed_functions[][0] and its arity in
   allowed_functions[][1] */
int             allowed_functions[MAX_NUM_FUNCTIONS][2];
int             allowed_output_functions[MAX_NUM_OUTPUT_FUNCTIONS];


/* globals defining the computational problems */
int					num_inputs[MAX_NUM_ITERATIONS];
int					num_outputs[MAX_NUM_ITERATIONS];
int					num_tests[MAX_NUM_ITERATIONS];
int					max_num_tests, max_num_inputs, max_num_outputs;
data_type***		data_inputs;
data_type***		data_outputs;

/* calculated global constants  */
int				bit_width[MAX_NUM_ITERATIONS];
int             num_genes;
double			perfect;
double			perfect_at_it[MAX_NUM_ITERATIONS];

/* global_diagnostic gen is a global that is use
   to allow capturing diagnostics about a program
   crash that occurs at a given generation */
int global_diagnostic_gen;

/* global to phenotype. This makes it easier as no calloc or realloc
   are needed to create new phenotypes during iteration. It is
   also much faster than calling malloc, calloc or realloc
   in the phenotype iteration loop */

node phenotype[MAX_NUM_NODES];

/* macros */
#define pow2(x) (1<<x)                                   /* returns 2 to power x (x>=0) */
#define getbit(decimal,nthbit) ((decimal>>nthbit) & 01)  /* gets nth bit */



/* function prototypes */
data_type  node_type(data_type in[MAX_NUM_CONNECTION_GENES_PER_NODE],
                     int function_gene);

void validate_command_line(int argc, char* argv[], char parfile[]);

int get_num_mutant(int num_genes, double per_cent_mutate);

void get_number_of_mutations(void);

void check_parameter_error(void);

void extract_node_function_string(char string[MAX_NUM_LETTERS]);

void get_parameters(char parfile[MAX_NUM_LETTERS]);

void write_cgp_info(char command[]);

data_type myfscanf(FILE* fp);

void myfprintf(FILE* fp, data_type datum, int linefeed);

void read_data_header(char datafile[MAX_NUM_LETTERS],
					  int iteration);

int get_max_array(int array[], int num_items);

void read_target_file_info(int* max_num_inputs,
						   int* max_num_outputs,
						   int* max_num_tests);

void read_data(char datafile[MAX_NUM_LETTERS], int iteration);

double define_perfect(void);

void fprint_a_chromosome(node* chromosome,
						 int   num_nodes_chromosome,
						 char  name[],
						 int   append);

void fprint_a_phenotype(node phenotype[MAX_NUM_NODES],
						int   num_nodes_chromosome,
						char  name[],
						int   append);

void print_a_chromosome(node* chromosome, int num_nodes_chromosome);

void print_a_phenotype(node phenotype[MAX_NUM_NODES], int num_nodes_phenotype)
;


void fprint_a_raw_chromosome(node* chromosome,
							 int   num_nodes_chromosome,
							 char  name[], int append);

int get_arity(int function_gene);


int get_output_nodes_genotype(node* chromosome,
					          int   num_nodes_chromosome,
					          int   output_nodes[MAX_NUM_NODES]);

int get_output_nodes_phenotype(int   num_nodes_phenotype,
					           int   output_nodes[MAX_NUM_NODES]);

void get_node_used_in_genotype(node* chromosome,
				   int   num_nodes_chromosome,
				   int   node_used[MAX_NUM_NODES]);

void get_node_used_in_phenotype(int   num_nodes_phenotype,
				                int   node_used[MAX_NUM_NODES]);

int node_is_sm(int the_node);

int get_nodes_to_process(int   num_nodes_phenotype,
						 ToDo  ToDoList[MAX_NUM_NODES],
						 int*  num_sm_nodes_to_process,
						 int   nodes_to_process[MAX_NUM_NODES],
						 int   iteration);


void fprint_active_genes(node* chromosome,
						 int   num_nodes_chromosome,
						 char  name[MAX_NUM_LETTERS]);

void fprint_active_phenes(int   num_nodes_phenotype,
						  char  name[MAX_NUM_LETTERS],
						  int   iteration);

void fprint_dot_graph(int num_nodes_phenotype, int iteration, int run, char name[50]);

void print_evolved_input_output_table(int num_nodes_to_process,
                                      int nodes_to_process[MAX_NUM_NODES],
                                      int run,
                                      int iteration);

void fprint_all_phenotypes(node* genotype,
						   char name[MAX_NUM_LETTERS],
						   int withdot, int run, int print_io);

node copy_node(node source);

void copy_chrom(node* source,
				node* destination,
				int   number_of_nodes);

void copy_chrom_to_phenotype(node* source,
				             int   number_of_nodes);


void copy_phen(node source[MAX_NUM_NODES],
			   node destination[MAX_NUM_NODES],
			   int  number_of_nodes);


void copy_pop_to_chrom(node** chromosomes,
					   node*  chromosome,
					   int    pop_member);

void copy_chrom_to_pop(node** chromosomes,
					   node*  chromosome,
					   int    pop_member);

void read_from_chrom(node** chromosomes);

void in_nodes(int       function_type,
			  int*      input_pointer,
			  int       the_node,
			  int       nodes_to_process[MAX_NUM_NODES],
			  data_type in[MAX_NUM_CONNECTION_GENES_PER_NODE],
              int       iteration,
			  int       fitness_test);

void out_nodes(int       function_type,
			   int*      output_pointer,
			   data_type output_register[MAX_NUM_OUTPUTS],
			   int       the_node,
			   int       nodes_to_process[MAX_NUM_NODES],
			   data_type in[MAX_NUM_CONNECTION_GENES_PER_NODE],
			   data_type output[MAX_NUM_NODES],
               int       iteration);

void sm_nodes(int       the_node,
			  int       nodes_to_process[MAX_NUM_NODES],
			  data_type in[MAX_NUM_CONNECTION_GENES_PER_NODE],
			  data_type	output[MAX_NUM_NODES]);

void comp_nodes(int       function_type,
				int       the_node,
			    int       nodes_to_process[MAX_NUM_NODES],
			    data_type in[MAX_NUM_CONNECTION_GENES_PER_NODE],
				data_type output[MAX_NUM_NODES]);

int function_type_is_input(int function_gene);
int function_type_is_output(int function_gene);
int function_type_is_sm(int function_gene);
int function_type_is_comp(int function_gene);


void decode_cgp(data_type cgp_outputs[MAX_NUM_OUTPUTS],
				int num_nodes_to_process,
				int nodes_to_process[MAX_NUM_NODES],
				int fitness_test,
				int iteration);




double correctness_test(data_type data_output,
						 data_type cgp_output,
						 int iteration);

double evaluate_cgp_outputs(data_type cgp_outputs[MAX_NUM_OUTPUTS],
							int test, int iteration);

int reached_perfect_score(double fit, double perfect);

void get_bit_width(int iteration);

double fitness_at_iteration(int num_nodes_to_process,
							int nodes_to_process[MAX_NUM_NODES],
							int iteration);

void insertionsort(int* array, int num_items);

void make_arguments_and_position_valid(int  num_nodes_phenotype,
									   int* sm_position,
									   int  sm_arguments[MAX_NUM_ARGUMENTS],
			                           int  true_sm_arguments[MAX_NUM_ARGUMENTS]);

int  get_size(int   num_nodes_phenotype,
			  int   true_sm_arguments[MAX_NUM_ARGUMENTS],
			  int   sm_operation);


void get_sm_arguments(ToDo  ToDoList[MAX_TODO], int which_sm,
					  int*  sm_position,
					  int   sm_arguments[MAX_NUM_ARGUMENTS],
					  int*  sm_operation);


void  DEL(int*  num_nodes_phenotype,
		  int sm_position, int   true_sm_arguments[MAX_NUM_ARGUMENTS]);

void  ADD(int*  num_nodes_phenotype,
		  int sm_position, int   true_sm_arguments[MAX_NUM_ARGUMENTS]);

void  MOV(int*  num_nodes_phenotype,
		  int sm_position, int   true_sm_arguments[MAX_NUM_ARGUMENTS]);

void  OVR(int*  num_nodes_phenotype,
		  int sm_position, int   true_sm_arguments[MAX_NUM_ARGUMENTS]);

void  DUP(int*  num_nodes_phenotype,
		  int sm_position, int   true_sm_arguments[MAX_NUM_ARGUMENTS]);

void  CRP(int*  num_nodes_phenotype,
		  int sm_position, int   true_sm_arguments[MAX_NUM_ARGUMENTS]);

void  CHC(int*  num_nodes_phenotype,
		  int sm_position, int   true_sm_arguments[MAX_NUM_ARGUMENTS]);

void  CHF(int*  num_nodes_phenotype,
		  int sm_position, int   true_sm_arguments[MAX_NUM_ARGUMENTS]);

void  CHA(int*  num_nodes_phenotype,
		  int sm_position, int   true_sm_arguments[MAX_NUM_ARGUMENTS]);

void generate_new_phenotype(int*   num_nodes_phenotype,
                            ToDo   ToDoList[MAX_TODO],
						    int    num_sm_nodes_to_process);

void generate_new_phenotype_verbose(int*   num_nodes_phenotype,
									ToDo   ToDoList[MAX_TODO],
						            int    num_sm_nodes_to_process,
									char   name[MAX_NUM_LETTERS],
									int    iteration);

double iterate_phenotype_and_accumulate_fitness(node* genotype, int* last_perfect_iteration);

double iterate_phenotype_and_accumulate_fitness_with_comments(node* genotype);

double fitness(node* genotype, int* last_perfect_iteration);

int get_connection_allele(int which_node);

int get_argument_allele(void);

int get_function_allele(void);

int get_output_function_allele(void);

node generate_a_random_node(int position);

void generate_a_random_chromosome(node* chromosome);

void initialise(node** chromosomes);

double  get_best_chromosome(node** chromosomes,
                            node*  best_chromosome,
						    double previous_best_fitness,
						    int    gen,
						    int*   last_perfect_iteration);

int get_valid_connection_gene(int num_connection_genes);

int get_valid_function_gene(int num_functions);

int get_valid_argument_gene(int num_argument_genes);

void mutate_a_connection_gene(node*  chromosome);

void mutate_a_function_gene(node*  chromosome);

void mutate_an_argument_gene(node*  chromosome);

void mutate_a_chromosome(node* chromosome,
						 int num_connection_mutations,
						 int num_function_mutations,
						 int num_argument_mutations);

void generate_new_pop_es(node** chromosomes,
						 node*  best_chromosome,
						 int num_connection_mutations,
						 int num_function_mutations,
						 int num_argument_mutations);

data_type* create_1d_datatype_space(int dimension1);

data_type** create_2d_datatype_space(int dimension1, int dimension2);

data_type*** create_3d_datatype_space(int dimension1, int dimension2, int dimension3);

void free_datatype_array2d(int dimension2, data_type** array2d);

void free_datatype_array3d(int dimension2, int dimension3, data_type*** array3d);

node* create_1d_array_of_nodes(int size);

node** create_2d_array_of_nodes(int num_horizontals, int num_verticals);

void free_array2d_of_nodes(int num_verticals, node** array2d);

void write_generation_to_screen(int gen_index);

void write_progress_info_to_screen(int generation, double fit, int last_perfect_iteration);

void write_progress_info_to_file(char prog[MAX_NUM_LETTERS],
								 int gen, double best_fit,
						         node* best_chromosome);

void check_if_improvement(double best_fit, double* previous_best_fit, int* best_gen, int gen,
						  char prog[MAX_NUM_LETTERS],
						  node* best_chromosome, int last_perfect_iteration);

void write_result_of_EA_to_file(int run, int bgen, double best_fit,
								node* best_chromosome);

double  EA(int *gen_of_best, int run,
		char prog[MAX_NUM_LETTERS]);


void setup_report_files(int run, char prog[MAX_NUM_LETTERS]);

void report_final_results(double av_fitness, double st_dev,
						  double best_of_best_fit,
						  double worst_of_best_fit,
						  double av_best_gen,
						  double av_best_gen_perfect,
						  int    num_perfect,
						  double av_num_evals_perfect);

double average(int num_items, double items[MAX_NUM_RUNS]);

double get_standard_deviation(int num_items, double average,
							  double items[MAX_NUM_RUNS]);

void run_EA(int num_runs_total);

#endif
