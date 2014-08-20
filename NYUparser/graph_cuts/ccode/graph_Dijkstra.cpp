
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dijkstraMaxFlow.h"
#include "instances_dijkstra.inc"


template <typename edgeweighttype, typename maxflowtype> 
DijkstraMaxFlow<edgeweighttype,maxflowtype>::DijkstraMaxFlow(int node_num_max, int edge_num_max, void (*err_function)(char *))
	: node_num(0),
	  max_flow_graph(NULL),
	  first_node(NULL),
	  last_node(NULL),
	  error_function(err_function)
{
	if (node_num_max < 16) node_num_max = 16;
	if (edge_num_max < 16) edge_num_max = 16;

	nodes = (node*) malloc(node_num_max*sizeof(node));
	arcs = (arc*) malloc(2*edge_num_max*sizeof(arc));
	if (!nodes || !arcs) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }

	node_last = nodes;
	node_max = nodes + node_num_max;
	arc_last = arcs;
	arc_max = arcs + 2*edge_num_max;
}

template <typename edgeweighttype, typename maxflowtype> 
DijkstraMaxFlow<edgeweighttype,maxflowtype>::~DijkstraMaxFlow()
{
	free(nodes);
	free(arcs);
}

template <typename edgeweighttype, typename maxflowtype> 
void DijkstraMaxFlow<edgeweighttype,maxflowtype>::reset()
{
	node_last = nodes;
	arc_last = arcs;
	node_num = 0;
	
	first_node == NULL;
	last_node == NULL;

	max_flow_graph = NULL;
	tempList.clear();
}

template <typename edgeweighttype, typename maxflowtype> 
void DijkstraMaxFlow<edgeweighttype,maxflowtype>::reallocate_nodes(int num)
{
	int node_num_max = (int)(node_max - nodes);
	node* nodes_old = nodes;

	node_num_max += node_num_max / 2;
	if (node_num_max < node_num + num) node_num_max = node_num + num;
	nodes = (node*) realloc(nodes_old, node_num_max*sizeof(node));
	if (!nodes) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }

	node_last = nodes + node_num;
	node_max = nodes + node_num_max;

	if (nodes != nodes_old)
	{
		arc* a;
		for (a=arcs; a<arc_last; a++)
		{
			a->head = (node*) ((char*)a->head + (((char*) nodes) - ((char*) nodes_old)));
		}
	}
}

template <typename edgeweighttype, typename maxflowtype> 
void DijkstraMaxFlow<edgeweighttype,maxflowtype>::reallocate_arcs()
{
	int arc_num_max = (int)(arc_max - arcs);
	int arc_num = (int)(arc_last - arcs);
	arc* arcs_old = arcs;

	arc_num_max += arc_num_max / 2; if (arc_num_max & 1) arc_num_max ++;
	arcs = (arc*) realloc(arcs_old, arc_num_max*sizeof(arc));
	if (!arcs) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }

	arc_last = arcs + arc_num;
	arc_max = arcs + arc_num_max;

	if (arcs != arcs_old)
	{
		node* i;
		arc* a;
		for (i=nodes; i<node_last; i++)
		{
			if (i->first) i->first = (arc*) ((char*)i->first + (((char*) arcs) - ((char*) arcs_old)));
		}
		for (a=arcs; a<arc_last; a++)
		{
			if (a->next) a->next = (arc*) ((char*)a->next + (((char*) arcs) - ((char*) arcs_old)));
			a->sister = (arc*) ((char*)a->sister + (((char*) arcs) - ((char*) arcs_old)));
		}
	}
}
