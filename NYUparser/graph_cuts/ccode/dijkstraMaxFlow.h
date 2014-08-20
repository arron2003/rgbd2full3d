#ifndef __DIJKSTRA_H__
#define __DIJKSTRA_H__

#include <string.h>
#include <queue>
#include <vector>

#include "graph.h"
#include "block.h"

#include <assert.h>

template <typename edgeweighttype, typename maxflowtype> class DijkstraMaxFlow
{
public:
	typedef int node_id;

	DijkstraMaxFlow(int node_num_max, int edge_num_max, void (*err_function)(char *) = NULL);

	// Destructor
	~DijkstraMaxFlow();

	// Adds node(s) to the graph. By default, one node is added (num=1); then first call returns 0, second call returns 1, and so on. 
	// If num>1, then several nodes are added, and node_id of the first one is returned.
	node_id add_node(int num = 1);

	void add_node_weight(node_id i, edgeweighttype _weight);

	void add_edge(node_id i, node_id j, edgeweighttype weight, edgeweighttype rev_weight);
	void add_mode_edge(node_id _i, node_id _j,edgeweighttype weight, edgeweighttype rev_weight, edgeweighttype mode, edgeweighttype rev_mode);

private:
	struct node;
	struct arc;

public:
	void reset();

private:
	// internal variables and functions

	struct node
	{
		arc			*first;		// first outcoming arc
		int node_identifier;
		node  *dijkstra_parent;
		maxflowtype distance;
		unsigned int is_fixed: 1;
		int when_added;
		int position_in_tempList;
		edgeweighttype weight;
		edgeweighttype mode;
	};

	typedef node* node_pointer;

	struct arc
	{
		node		*head;		// node the arc points to
		arc			*next;		// next arc with the same originating node
		arc			*sister;	// reverse arc

		edgeweighttype weight;		// edge weight
		edgeweighttype mode;
	};

	node				*nodes, *node_last, *node_max; // node_last = nodes+node_num, node_max = nodes+node_num_max;
	arc					*arcs, *arc_last, *arc_max; // arc_last = arcs+2*edge_num, arc_max = arcs+2*edge_num_max;

	int					node_num;

	void	(*error_function)(char *);	// this function is called if a error occurs,
										// with a corresponding error message
										// (or exit(1) is called if it's NULL)

	void reallocate_nodes(int num); // num is the number of new nodes
	void reallocate_arcs();

/////////////////////////// Functions for dijkstra algorithm //////////////////////////////

public:
	typedef Graph<maxflowtype,maxflowtype,maxflowtype> GraphType;
	void add_starting_node(node_id i);
	void add_terminal_node(node_id i);
	void add_max_flow_graph(GraphType * m_f_graph);
	void add_figure_coordinates(int _x, int _y);
	void set_nrows_ncols(int _nrows,int _ncols);

	double GetWhenAdded(node_id i);
	int GetParentDijkstra(node_id i);
	
	double simple_dijkstra();

	maxflowtype dijkstra();
	double node_dijkstra(int useMax = 0);

	double mode_dijkstra();
private: 
	node_pointer first_node, last_node; 
	int nrows, ncols; //image dimensions

	//std::vector<node_pointer> terminal_nodes; 

	/*pointer to a graph where max-flow algorithm will be run
	The grahp should be constructed a-priori and the first and last node
	should be already included as hard-constraints in the graph
	
	There should also exist a correspondence between the nodes in both graphs
	*/
	GraphType *max_flow_graph;

	void put_hard_constraint(node_pointer node, int is_source, int with_parent = 1);
	void remove_hard_constraint(node_pointer node, int is_source, int with_parent = 1);

	void put_hard_constraint(node_pointer node, int is_source, node_pointer stop_node = NULL);
	void remove_hard_constraint(node_pointer node, int is_source, node_pointer stop_node = NULL);

	void ini_dijkstra();
	int in_segmentation(node_pointer node, int segment);

	class Coordinates
	{	
	private:
		int x;
		int y;
	public:
		Coordinates(int _x, int _y);
		int getX();
		int getY();
	};

	std::vector<Coordinates> figure_coordinates;


// for the heap (called tempList)
private:
	std::vector<node_pointer> tempList;
	
	void percolateDown(node_pointer node);
	void percolateUp(node_pointer node);

	void add_to_temp_list(node_pointer node);
	void remove_from_temp_list(node_pointer node);
	
	node_pointer get_next_best_node();
	
};


///////////////////////////////////////
// Implementation - inline functions //
///////////////////////////////////////



template <typename edgeweighttype, typename maxflowtype>
inline typename DijkstraMaxFlow<edgeweighttype,maxflowtype>::node_id DijkstraMaxFlow<edgeweighttype,maxflowtype>::add_node(int num)
{
	assert(num > 0);

	if (node_last + num > node_max) reallocate_nodes(num);

	if (num == 1)
	{
		node_last -> first = NULL;
		node_last ++;
		return node_num ++;
	}
	else
	{
		memset(node_last, 0, num*sizeof(node));

		node_id i = node_num;
		node_num += num;
		node_last += num;
		return i;
	}
}

template <typename edgeweighttype, typename maxflowtype>
	inline void DijkstraMaxFlow<edgeweighttype,maxflowtype>::add_edge(node_id _i, node_id _j,edgeweighttype weight, edgeweighttype rev_weight)
{
	assert(_i >= 0 && _i < node_num);
	assert(_j >= 0 && _j < node_num);
	assert(_i != _j);
	assert(weight >= 0);
	assert(rev_weight >= 0);

	if (arc_last == arc_max) reallocate_arcs();

	arc *a = arc_last ++;
	arc *a_rev = arc_last ++;

	node* i = nodes + _i;
	node* j = nodes + _j;

	a -> sister = a_rev;
	a_rev -> sister = a;
	a -> next = i -> first;
	i -> first = a;
	a_rev -> next = j -> first;
	j -> first = a_rev;
	a -> head = j;
	a_rev -> head = i;
	a -> weight = weight;
	a_rev -> weight = rev_weight;
}



template <typename edgeweighttype, typename maxflowtype>
	inline void DijkstraMaxFlow<edgeweighttype,maxflowtype>::add_mode_edge(node_id _i, node_id _j,edgeweighttype weight, edgeweighttype rev_weight, edgeweighttype mode, edgeweighttype rev_mode)
{
	assert(_i >= 0 && _i < node_num);
	assert(_j >= 0 && _j < node_num);
	assert(_i != _j);
	assert(weight >= 0);
	assert(rev_weight >= 0);

	if (arc_last == arc_max) reallocate_arcs();

	arc *a = arc_last ++;
	arc *a_rev = arc_last ++;

	node* i = nodes + _i;
	node* j = nodes + _j;

	a -> sister = a_rev;
	a_rev -> sister = a;
	a -> next = i -> first;
	i -> first = a;
	a_rev -> next = j -> first;
	j -> first = a_rev;
	a -> head = j;
	a_rev -> head = i;
	a -> weight = weight;
	a_rev -> weight = rev_weight;
	a -> mode = mode;
	a_rev -> mode = rev_mode;

}


template <typename edgeweighttype, typename maxflowtype> 
inline void DijkstraMaxFlow<edgeweighttype,maxflowtype>::add_node_weight(node_id _i, edgeweighttype _weight){

	assert(_i >= 0 && _i < node_num);
	nodes[_i].weight = _weight;
}


#endif
