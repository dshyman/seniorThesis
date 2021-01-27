/* node-function.c
   Julian F. Miller (c) 2013
   version 1.2
   Dept. of Electronics, University of York, UK

   node_type calculates the output of a node given
   the data provided in the array in
*/

#include "smcgp.h"
#include <math.h>


#ifdef DATA_IS_UNSIGNED_INT

data_type  node_type(data_type in[MAX_NUM_CONNECTION_GENES_PER_NODE],
                     int function_gene)
{
   /* assumes 32 bit operations */
   data_type result;

   switch(function_gene)
   {
		case 16:  /* 0 */
			result = 0;
			break;
		case 17:  /* 1 */
			result = MAXNUM;
			break;
		case 18: /* wire  */
			result = in[0];
			break;
		case 19: /* wire */
			result = in[1];
			break;
		case 20: /* NOT */
			result = ~in[0];
			break;
		case 21: /* NOT */
			result = ~in[1];
			break;
		case 22: /* two input gate functions */
			result = (in[0] & in[1]);   /* AND(a, b) */
			break;
		case 23:
			result = (in[0] & ~in[1]);  /* AND (a,NOT(b)) */
			break;
		case 24:
			result = (~in[0] & in[1]);  /* AND(NOT(a), b) */
			break;
		case 25:
			result = (~in[0] & ~in[1]); /* AND(NOT(a), NOT(b)) - i.e. NOR(a, b) */
			break;
		case 26:
			result = (in[0]^in[1]);     /* XOR(a, b) */
			break;
		case 27:
			result = (~in[0]^in[1]);    /* XNOR(a, b) */
			break;
		case 28:
			result = (in[0] | in[1]);   /* OR(a, b) */
			break;
		case 29:
			result = (in[0] | ~in[1]);  /* OR(a, NOT(b)) */
			break;
		case 30:
			result = (~in[0] | in[1]);  /* OR(NOT(a), b) */
			break;
		case 31:
			result = (~in[0] | ~in[1]); /* OR(NOT(a), NOT(b)) - i.e. NAND(a,b) */
			break;
		case 32:  /* mux functions - i.e IF statements */
			result = ((in[0] & ~in[2]) | (in[1] & in[2]));
			break;
		case 33:
			result = ((in[0] & ~in[2]) | (~in[1] & in[2]));
			break;
		case 34:
			result = ((~in[0] & ~in[2]) | (in[1] & in[2]));
			break;
		case 35:
			result = ((~in[0] & ~in[2]) | (~in[1] & in[2]));
		default: /* functions 0-15 */
			result = in[0];
	}
	return result;
}
#endif

#ifdef DATA_IS_INT
/* Currently this is defined to be similar to the case
   when data is double.
*/
data_type  node_type(data_type in[MAX_NUM_CONNECTION_GENES_PER_NODE],
                     int function_gene)
{
   data_type result;

   switch(function_gene)
   {
		case 16:  /*   absolute  */
			result = fabs(in[0]);
			break;
		case 17:  /*   protected sqrt  */
			result = sqrt(fabs(in[0]));
			break;
		case 18:	/* protected reciprocal */
			if (fabs(in[0]) < 0.0001)
				result = in[0];
			else
				result = 1.0/in[0];
			break;
		case 19:  /*   sin  */
			result = sin(in[0]);
			break;
		case 20:  /*   cos  */
			result = cos(in[0]);
			break;
		case 21:  /*   protected tan  */
			if (cos(in[0]) < 0.0001)
				result = in[0];
			else
				result = tan(in[0]);
			break;
		case 22:  /*   exp  */
			result = exp(in[0]);
			break;
		case 23:  /*   sinh  */
			result = sinh(in[0]);
			break;
		case 24:  /*   cosh  */
			result = cosh(in[0]);
			break;
		case 25:  /*   tanh  */
			result = tanh(in[0]);
			break;
		case 26:  /*   protected natural log  */
			if (in[0] < 0.0001)
				result = in[0];
			else
				result = log(fabs(in[0]));
			break;
		case 27:  /*   protected  log to base 10  */
			if (in[0] < 0.0001)
				result = in[0];
			else
				result = log10(fabs(in[0]));
			break;
		case 28:  /*   sin (a+b) */
			result = sin(in[0]+in[1]);
			break;
		case 29:  /*   cos(a+b)  */
			result = cos(in[0]+in[1]);
			break;
		case 30:  /*   protected remainder on division  */
			if ((int)fabs(in[1]) == 0)
				result = in[0];
			else
			{
				if (in[0] > MAXNUM/10)
					in[0] = (double) MAXNUM/10;
				if (in[1] > MAXNUM/10)
					in[1] = (double) MAXNUM/10;
				result = ((int)in[0]) % ((int) in[1]);
				break;
			}
		case 31:	/* power */
				result = pow(fabs(in[0]),in[1]);
			break;
		case 32:  /*   addition  */
			result = in[0] + in[1];
			break;
		case 33:	 /* subtraction */
			result = in[0] - in[1];
			break;
		case 34: /* multiplication */
				result = in[0]*in[1];
			break;
		case 35:	/* protected division */
			if (fabs(in[1]) < 0.0001)
				result = in[0];
			else
				result = in[0]/in[1];
			break;
		default:
			result = in[0];
   }

   return result;
}

#endif

#ifdef DATA_IS_DOUBLE

data_type  node_type(data_type in[MAX_NUM_CONNECTION_GENES_PER_NODE],
                     int function_gene)
{
   data_type result;

   switch(function_gene)
   {
		case 16:  /*   absolute  */
			result = fabs(in[0]);
			break;
		case 17:  /*   protected sqrt  */
			result = sqrt(fabs(in[0]));
			break;
		case 18:	/* protected reciprocal */
			if (fabs(in[0]) < 0.0001)
				result = in[0];
			else
				result = 1.0/in[0];
			break;
		case 19:  /*   sin  */
			result = sin(in[0]);
			break;
		case 20:  /*   cos  */
			result = cos(in[0]);
			break;
		case 21:  /*   protected tan  */
			if (cos(in[0]) < 0.0001)
				result = in[0];
			else
				result = tan(in[0]);
			break;
		case 22:  /*   exp  */
			result = exp(in[0]);
			break;
		case 23:  /*   sinh  */
			result = sinh(in[0]);
			break;
		case 24:  /*   cosh  */
			result = cosh(in[0]);
			break;
		case 25:  /*   tanh  */
			result = tanh(in[0]);
			break;
		case 26:  /*   protected natural log  */
			if (in[0] < 0.0001)
				result = in[0];
			else
				result = log(fabs(in[0]));
			break;
		case 27:  /*   protected  log to base 10  */
			if (in[0] < 0.0001)
				result = in[0];
			else
				result = log10(fabs(in[0]));
			break;
		case 28:  /*   sin (a+b) */
			result = sin(in[0]+in[1]);
			break;
		case 29:  /*   cos(a+b)  */
			result = cos(in[0]+in[1]);
			break;
		case 30:  /*   hypoteneuse  */
			result = sqrt(in[0]*in[0]+in[1]*in[1]);
			break;
		case 31:	/* power */
				result = pow(fabs(in[0]),in[1]);
			break;
		case 32:  /*   addition  */
			result = in[0] + in[1];
			break;
		case 33:	 /* subtraction */
			result = in[0] - in[1];
			break;
		case 34: /* multiplication */
				result = in[0]*in[1];
			break;
		case 35:	/* protected division */
			if (fabs(in[1]) < 0.0001)
				result = in[0];
			else
				result = in[0]/in[1];
			break;
		default:
			result = in[0];
   }

   return result;
}

#endif
