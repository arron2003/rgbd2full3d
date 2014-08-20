/* energy.h */
/* Vladimir Kolmogorov (vnk@cs.cornell.edu), 2003. 
   2007: adapted to maxflow-v3.0 */

/*
	This software minimizes certain energy functions of binary variables, as described in 

		What Energy Functions can be Minimized via Graph Cuts?
		Vladimir Kolmogorov and Ramin Zabih. 
		In IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI), February 2004. 

	It uses maxflow algorithm described in 

		An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Vision
		Yuri Boykov and Vladimir Kolmogorov.
		In IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI), September 2004.

	Note: The fact that submodular functions can be minimized via graph cuts (min-cut/max-flow algorithm) 
	is well-known in the optimization literature, see e.g.

		P.L. Hammer. "Some network flow problems solved with pseudo-Boolean programming." 
		Operations Research 13:388-399, 1965. 

	Graph construction for submodular functions with triple cliques was given in

		A. Billionnet and M. Minoux. "Maximizing a supermodular pseudo-boolean function: 
		a polynomial algorithm for supermodular cubic functions." Discrete Appl. Math. 12(1):1-11, 1985 

	For a review of pseudo-boolean optimization literature, see

		E. Boros and P.L. Hammer. "Pseudo-boolean optimization." Discrete Appl. Math. 123:155-225, 2002 
	==================================================================================================

	More specifically, the computes the global minimum of a function E of binary
	variables x_1, ..., x_n which can be written as a sum of terms involving
	at most three variables at a time:

		E(x_1, ..., x_n) = \sum_{i}     E^{i}    (x_i)
						+ \sum_{i,j}   E^{i,j}  (x_i, x_j)
						+ \sum_{i,j,k} E^{i,j,k}(x_i, x_j, x_k)

	The method works only if each term is "submodular". Definitions of submodularity
	for terms E^{i}, E^{i,j}, E^{i,j,k} are given below as comments to functions
	add_term1(), add_term2(), add_term3(). 

	In order to use it, you will also need a MAXFLOW software which can be
	obtained from http://www.adastral.ucl.ac.uk/~vnk/software.html

	NOTE: This software minimizes functions of BINARY variables only.
	However, it can also be used for minimizing certain functions of non-binary
	(multi-label) variables via a sequence of binary moves (alpha-expansion, 
	alpha-beta swap, k-jumps, etc.) as proposed in

		Efficient Approximate Energy Minimization via Graph Cuts 
		Yuri Boykov, Olga Veksler, Ramin Zabih, 
		IEEE transactions on PAMI, vol. 20, no. 12, p. 1222-1239, November 2001.

	IF YOU USE THIS SOFTWARE FOR IMPLEMENTING ALPHA-EXPANSION OR ALPHA-BETA SWAP
	ALGORITHM, YOU SHOULD CITE THIS PAPER IN ANY RESULTING PUBLICATION.

	Also note that an implementation of minimization techniques for non-binary variables
	can be downloaded from O. Veksler's homepage: http://www.csd.uwo.ca/faculty/olga/code.html .

	------------------------------------------------------------------------

	Example usage
	(Minimizes the following function of 3 binary variables:
	E(x, y, z) = x - 2*y + 3*(1-z) - 4*x*y + 5*|y-z|):

	///////////////////////////////////////////////////

	void main()
	{
		// Minimize the following function of 3 binary variables:
		// E(x, y, z) = x - 2*y + 3*(1-z) - 4*x*y + 5*|y-z|

		Energy<int,int> *e = new Energy<int,int>(3, 2); // 3 vars, 2 pairwise terms

		e -> add_variable(3);

		e -> add_term1(0,       0,  1);  // add term x 
		e -> add_term1(1,       0, -2);  // add term -2*y
		e -> add_term1(2,       3,  0);  // add term 3*(1-z)

		e -> add_term2(0, 1,    0, 0, 0, -4);  // add term -4*x*y
		e -> add_term2(1, 2,    0, 5, 5, 0);   // add term 5*|y-z|

		int Emin = e -> minimize();
		
		printf("Minimum = %d\n", Emin);
		printf("Optimal solution:\n");
		printf("x = %d\n", e->get_var(0));
		printf("y = %d\n", e->get_var(1));
		printf("z = %d\n", e->get_var(2));

		delete e;
	}

	///////////////////////////////////////////////////
*/

#ifndef __ENERGY_H__
#define __ENERGY_H__

#include <assert.h>
#include "graph.h"

/*
	Value is a type of a value in a single term.
	TotalValue is a type of a value of the total energy.
*/
template <typename Value, typename TotalValue> class Energy : Graph<Value,Value,TotalValue>
{
public:
	typedef node_id Var; // int

	/* interface functions */

	/* Constructor. 
	   For description of var_num_max and edge_num_max see the constructor of Graph in graph.h.

	   Optional argument is the pointer to the
	   function which will be called if an error occurs;
	   an error message is passed to this function. If this
	   argument is omitted, exit(1) will be called. */
	Energy(int var_num_max, int edge_num_max, void (*err_function)(char *) = NULL);

	/* Destructor */
	~Energy();

	/* Adds 'num' new binary variables. 
	   By default, one variable is added (num=1); then first call returns 0, second call returns 1, and so on. 
	   If num>1, then several variables are added, and node_id of the first one is returned. */
	Var add_variable(int num = 1);

	/* Adds a constant E to the energy function */
	void add_constant(Value E);

	/* Adds a new term E(x) of one binary variable
	   to the energy function, where
	       E(0) = E0, E(1) = E1
	   E0 and E1 can be arbitrary */
	void add_term1(Var x,
	               Value E0, Value E1);

	/* Adds a new term E(x,y) of two binary variables
	   to the energy function, where
	       E(0,0) = E00, E(0,1) = E01
	       E(1,0) = E10, E(1,1) = E11
	   The term must be submodular, i.e. E00 + E11 <= E01 + E10 */
	void add_term2(Var x, Var y,
	               Value E00, Value E01,
	               Value E10, Value E11);

	/* Adds a new term E(x,y,z) of three binary variables
	   to the energy function, where
	       E(0,0,0) = E000, E(0,0,1) = E001
	       E(0,1,0) = E010, E(0,1,1) = E011
	       E(1,0,0) = E100, E(1,0,1) = E101
	       E(1,1,0) = E110, E(1,1,1) = E111
	   The term must be submodular. It means that if one
	   of the variables is fixed (for example, y=1), then
	   the resulting function of two variables must be submodular.
	   Since there are 6 ways to fix one variable
	   (3 variables times 2 binary values - 0 and 1),
	   this is equivalent to 6 inequalities */
	void add_term3(Var x, Var y, Var z,
	               Value E000, Value E001,
	               Value E010, Value E011,
	               Value E100, Value E101,
	               Value E110, Value E111);

	/* After the energy function has been constructed,
	   call this function to minimize it.
	   Returns the minimum of the function */
	TotalValue minimize();

	/* After 'minimize' has been called, this function
	   can be used to determine the value of variable 'x'
	   in the optimal solution.
	   Returns either 0 or 1 */
	int get_var(Var x);

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

private:
	/* internal variables and functions */

	TotalValue	Econst;
	void		(*error_function)(char *);	/* this function is called if a error occurs,
											with a corresponding error message
											(or exit(1) is called if it's NULL) */
};















/***********************************************************************/
/************************  Implementation ******************************/
/***********************************************************************/

template <typename Value, typename TotalValue> 
	inline Energy<Value,TotalValue>::Energy(int var_num_max, int edge_num_max, void (*err_function)(char *)) : Graph<Value,Value,TotalValue>(var_num_max, edge_num_max, err_function)
{
	Econst = 0;
	error_function = err_function;
}

template <typename Value, typename TotalValue> 
	inline Energy<Value,TotalValue>::~Energy() {}

template <typename Value, typename TotalValue> 
	inline typename Energy<Value,TotalValue>::Var Energy<Value,TotalValue>::add_variable(int num) { return add_node(num); }

template <typename Value, typename TotalValue> 
	inline void Energy<Value,TotalValue>::add_constant(Value A) { Econst += A; }

template <typename Value, typename TotalValue> 
	inline void Energy<Value,TotalValue>::add_term1(Var x,
                              Value A, Value B)
{
	add_tweights(x, B, A);
}

template <typename Value, typename TotalValue> 
	inline void Energy<Value,TotalValue>::add_term2(Var x, Var y,
                              Value A, Value B,
                              Value C, Value D)
{
	/* 
	   E = A A  +  0   B-A
	       D D     C-D 0
	   Add edges for the first term
	*/
	add_tweights(x, D, A);
	B -= A; C -= D;

	/* now need to represent
	   0 B
	   C 0
	*/

	assert(B + C >= 0); /* check regularity */
	if (B < 0)
	{
		/* Write it as
		   B B  +  -B 0  +  0   0
		   0 0     -B 0     B+C 0
		*/
		add_tweights(x, 0, B); /* first term */
		add_tweights(y, 0, -B); /* second term */
		add_edge(x, y, 0, B+C); /* third term */
	}
	else if (C < 0)
	{
		/* Write it as
		   -C -C  +  C 0  +  0 B+C
		    0  0     C 0     0 0
		*/
		add_tweights(x, 0, -C); /* first term */
		add_tweights(y, 0, C); /* second term */
		add_edge(x, y, B+C, 0); /* third term */
	}
	else /* B >= 0, C >= 0 */
	{
		add_edge(x, y, B, C);
	}
}

template <typename Value, typename TotalValue> 
	inline void Energy<Value,TotalValue>::add_term3(Var x, Var y, Var z,
                              Value E000, Value E001,
                              Value E010, Value E011,
                              Value E100, Value E101,
                              Value E110, Value E111)
{
	register Value pi = (E000 + E011 + E101 + E110) - (E100 + E010 + E001 + E111);
	register Value delta;
	register Var u;

	if (pi >= 0)
	{
		Econst += E111 - (E011 + E101 + E110);

		add_tweights(x, E101, E001);
		add_tweights(y, E110, E100);
		add_tweights(z, E011, E010);

		delta = (E010 + E001) - (E000 + E011); /* -pi(E[x=0]) */
		assert(delta >= 0); /* check regularity */
		add_edge(y, z, delta, 0);

		delta = (E100 + E001) - (E000 + E101); /* -pi(E[y=0]) */
		assert(delta >= 0); /* check regularity */
		add_edge(z, x, delta, 0);

		delta = (E100 + E010) - (E000 + E110); /* -pi(E[z=0]) */
		assert(delta >= 0); /* check regularity */
		add_edge(x, y, delta, 0);

		if (pi > 0)
		{
			u = add_variable();
			add_edge(x, u, pi, 0);
			add_edge(y, u, pi, 0);
			add_edge(z, u, pi, 0);
			add_tweights(u, 0, pi);
		}
	}
	else
	{
		Econst += E000 - (E100 + E010 + E001);

		add_tweights(x, E110, E010);
		add_tweights(y, E011, E001);
		add_tweights(z, E101, E100);

		delta = (E110 + E101) - (E100 + E111); /* -pi(E[x=1]) */
		assert(delta >= 0); /* check regularity */
		add_edge(z, y, delta, 0);

		delta = (E110 + E011) - (E010 + E111); /* -pi(E[y=1]) */
		assert(delta >= 0); /* check regularity */
		add_edge(x, z, delta, 0);

		delta = (E101 + E011) - (E001 + E111); /* -pi(E[z=1]) */
		assert(delta >= 0); /* check regularity */
		add_edge(y, x, delta, 0);

		u = add_variable();
		add_edge(u, x, -pi, 0);
		add_edge(u, y, -pi, 0);
		add_edge(u, z, -pi, 0);
		add_tweights(u, -pi, 0);
	}
}

template <typename Value, typename TotalValue> 
	inline typename TotalValue Energy<Value,TotalValue>::minimize() { return Econst + maxflow(); }

template <typename Value, typename TotalValue> 
	inline int Energy<Value,TotalValue>::get_var(Var x) { return (int)what_segment(x); }

#endif

