#include <cmath>
#include <map>
#include <queue>
#include <vector>
#include <algorithm> 

#include "graph.h"
#include "dijkstraMaxFlow.h"
#include "instances_dijkstra.inc"
#include "mex.h"



#define INFINITY_D ((((unsigned)-1)/2))
#define hard_constraint 10000000


template <typename edgeweighttype, typename maxflowtype> 
void DijkstraMaxFlow<edgeweighttype,maxflowtype>::add_starting_node(node_id i)
{
	if(i < 0 || i >= node_num){
		printf("First node is out of limits \n");
	}
	else{
		first_node = &nodes[i];
}
}


template <typename edgeweighttype, typename maxflowtype> 
void DijkstraMaxFlow<edgeweighttype,maxflowtype>::add_terminal_node(node_id i)
{
	//assert(i >= 0 && i < node_num);
	if (i>=0 && i<node_num){
		last_node = &nodes[i];
	}
	else
		last_node = NULL;

}

template <typename edgeweighttype, typename maxflowtype> 
void DijkstraMaxFlow<edgeweighttype,maxflowtype>::add_max_flow_graph(GraphType * m_f_graph){
	max_flow_graph = m_f_graph;
}


template <typename edgeweighttype, typename maxflowtype> 
double DijkstraMaxFlow<edgeweighttype,maxflowtype>::GetWhenAdded(node_id i){
	return (double) nodes[i].distance;
}


template <typename edgeweighttype, typename maxflowtype> 
int DijkstraMaxFlow<edgeweighttype,maxflowtype>::GetParentDijkstra(node_id i){
	if (nodes[i].dijkstra_parent)
		return nodes[i].dijkstra_parent->node_identifier;
	else
		return -1;
		
}

template <typename edgeweighttype, typename maxflowtype> 
void DijkstraMaxFlow<edgeweighttype,maxflowtype>::put_hard_constraint(node_pointer node, int is_source, int with_parent){
	node_pointer i;	
	int node_identifier;
	int figure_node;
	typename std::vector<Coordinates>::iterator it;

	i = node;
	while( i ){
		node_identifier = i->node_identifier;
		for ( it=figure_coordinates.begin() ; it < figure_coordinates.end(); it++ ){
			figure_node = node_identifier - ((*it).getY()) + nrows*((*it).getX());
			if(figure_node >= 0 && figure_node < node_num){
				if(is_source){
					max_flow_graph -> add_tweights( figure_node, hard_constraint ,0 );
				}
				else{
					max_flow_graph -> add_tweights( figure_node, 0, hard_constraint );	
				}
				max_flow_graph -> mark_node(figure_node);
			}
		}
		i = i->dijkstra_parent;
		if (with_parent != 1) break;
	}
}

template <typename edgeweighttype, typename maxflowtype> 
void DijkstraMaxFlow<edgeweighttype,maxflowtype>::put_hard_constraint(node_pointer node, int is_source, node_pointer stop_node){
	node_pointer i;	
	int node_identifier;
	int figure_node;
	typename std::vector<Coordinates>::iterator it;

	i = node;
	while( i != stop_node && i != NULL){
		node_identifier = i->node_identifier;
		for ( it=figure_coordinates.begin() ; it < figure_coordinates.end(); it++ ){
			figure_node = node_identifier - ((*it).getY()) + nrows*((*it).getX());
			if(figure_node >= 0 && figure_node < node_num){
				if(is_source){
					max_flow_graph -> add_tweights( figure_node, hard_constraint ,0 );
				}
				else{
					max_flow_graph -> add_tweights( figure_node, 0, hard_constraint );	
				}
				max_flow_graph -> mark_node(figure_node);
			}
		}
		i = i->dijkstra_parent;
	}
}


template <typename edgeweighttype, typename maxflowtype> 
void DijkstraMaxFlow<edgeweighttype,maxflowtype>::remove_hard_constraint(node_pointer node, int is_source, int with_parent){
	node_pointer i;	
	arc* a;
	int node_identifier;
	int figure_node;
	typename std::vector<Coordinates>::iterator it;

	i = node;
	while( i ){
		node_identifier = i->node_identifier;
		for ( it=figure_coordinates.begin() ; it < figure_coordinates.end(); it++ ){
			figure_node = node_identifier - ((*it).getY()) + nrows*((*it).getX());
			if(figure_node >= 0 && figure_node < node_num){
				if(is_source){
					max_flow_graph -> add_tweights( figure_node, -hard_constraint , 0 );
				}
				else{
					max_flow_graph -> add_tweights(figure_node,0, -hard_constraint );
				}
				max_flow_graph -> mark_node(figure_node);
			}
		}
		i = i->dijkstra_parent;
		if (with_parent != 1) break;
	}
}

template <typename edgeweighttype, typename maxflowtype> 
void DijkstraMaxFlow<edgeweighttype,maxflowtype>::remove_hard_constraint(node_pointer node, int is_source, node_pointer stop_node){
	node_pointer i;	
	arc* a;
	int node_identifier;
	int figure_node;
	typename std::vector<Coordinates>::iterator it;

	i = node;
	while( i != stop_node && i != NULL){
		node_identifier = i->node_identifier;
		for ( it=figure_coordinates.begin() ; it < figure_coordinates.end(); it++ ){
			figure_node = node_identifier - ((*it).getY()) + nrows*((*it).getX());
			if(figure_node >= 0 && figure_node < node_num){
				if(is_source){
					max_flow_graph -> add_tweights( figure_node, -hard_constraint , 0 );
				}
				else{
					max_flow_graph -> add_tweights(figure_node,0, -hard_constraint );
				}
				max_flow_graph -> mark_node(figure_node);
			}
		}
		i = i->dijkstra_parent;
	}
}

template <typename edgeweighttype, typename maxflowtype> 
void DijkstraMaxFlow<edgeweighttype,maxflowtype>::ini_dijkstra(){

	node_pointer i;
	int k = 0;
	for (i=nodes; i<node_last; i++)
	{	
		i -> is_fixed = 0;
		i -> distance = (maxflowtype) INFINITY_D;
		i -> node_identifier = k;
		i -> when_added = 0;
		i -> position_in_tempList = -1;
		i -> dijkstra_parent = NULL;
		i -> mode = (maxflowtype) INFINITY_D;
		k++;
	}

	if((int)tempList.capacity() < k){
		tempList.reserve(k);
	}
}


template <typename edgeweighttype, typename maxflowtype> 
double DijkstraMaxFlow<edgeweighttype,maxflowtype>::node_dijkstra(int useMax)
{
	node_pointer current, current_neighboor;
	arc *a;
	
	double curr_distance;

	int k;

	ini_dijkstra();

	first_node -> distance = first_node -> weight;
	first_node -> dijkstra_parent = NULL;
	add_to_temp_list(first_node);
	k = 0;
	while (1){
		k++;
		current = get_next_best_node();
		if(current == NULL) break;
		current -> when_added = k;
		current -> is_fixed = 1;
		for (a = current->first; a; a=a->next)
		{
			current_neighboor = a -> head;
			if(current_neighboor -> is_fixed !=1){
				if(useMax){
					curr_distance = std::max((double) current -> distance,(double) current_neighboor -> weight);
				}
				else{
					curr_distance = current -> distance + (double) current_neighboor -> weight;
				}
				if ((current_neighboor -> distance) > curr_distance){
					current_neighboor -> dijkstra_parent = current;
					current_neighboor -> distance = curr_distance;
					add_to_temp_list(current_neighboor);
				}
			}
		}
		if (last_node != NULL && last_node->is_fixed) break;
	}
	if (last_node == NULL)
		return 0.0;
	else
		return last_node -> distance;

}



template <typename edgeweighttype, typename maxflowtype> 
double DijkstraMaxFlow<edgeweighttype,maxflowtype>::simple_dijkstra()
{

	node_pointer current, current_neighboor;
	arc *a;
	
	double curr_distance;

	int k;

	ini_dijkstra();

	first_node -> distance = 0.0;
	first_node -> dijkstra_parent = NULL;
	add_to_temp_list(first_node);
	k = 0;
	while (1){
		k++;
		current = get_next_best_node();
		if(current == NULL) break;
		current -> when_added = k;
		current -> is_fixed = 1;
		for (a = current->first; a; a=a->next)
		{
			current_neighboor = a -> head;

			if(current_neighboor -> is_fixed !=1){
				curr_distance = current -> distance + (double) a -> weight;
				if ((current_neighboor -> distance) > curr_distance){
					current_neighboor -> dijkstra_parent = current;
					current_neighboor -> distance = curr_distance;
					add_to_temp_list(current_neighboor);
				}
			}
		}
		if (last_node != NULL && last_node->is_fixed) break;
	}
	if (last_node == NULL)
		return 0.0;
	else
		return last_node -> distance;

}


template <typename edgeweighttype, typename maxflowtype> 
double DijkstraMaxFlow<edgeweighttype,maxflowtype>::mode_dijkstra()
{

	node_pointer current, current_neighboor;
	arc *a;
	
	double curr_distance;

	int k;

	ini_dijkstra();

	first_node -> distance = 0.0;
	first_node -> dijkstra_parent = NULL;
	first_node -> mode = 0.0;
	add_to_temp_list(first_node);
	k = 0;
	while (1){
		k++;
		current = get_next_best_node();
		if(current == NULL) break;
		current -> when_added = k;
		current -> is_fixed = 1;
		for (a = current->first; a; a=a->next)
		{
			current_neighboor = a -> head;

			if(current_neighboor -> is_fixed !=1){
				curr_distance = current -> distance + (double) a -> weight + current->mode;
				if ((current_neighboor -> distance) > curr_distance){
//					current_neighboor -> mode = std::min((double) current_neighboor -> mode, (double)current->mode + (double) a-> mode);
					current_neighboor -> mode = std::min((double) current_neighboor -> mode, std::max((double)current->mode,(double) a-> mode));

					current_neighboor -> dijkstra_parent = current;

					current_neighboor -> distance = curr_distance;
					
					add_to_temp_list(current_neighboor);
				}
				else{
					if(current_neighboor -> mode > current->mode + (double) a-> mode)
					{
						current_neighboor -> mode = current->mode + (double) a-> mode;
					}
				}
			}
		}
		if (last_node != NULL && last_node->is_fixed) break;
	}
	if (last_node == NULL)
		return 0.0;
	else
		return last_node -> distance;

}


template <typename edgeweighttype, typename maxflowtype> 
maxflowtype DijkstraMaxFlow<edgeweighttype,maxflowtype>::dijkstra()
{

	node_pointer current, current_neighboor, node_seg, parent_node;
	arc *a;
	maxflowtype first_node_flow, curr_flow, neigh_flow;
	int k = 0;
	int is_source;
	int segment;

	std::queue<node_pointer> nodes_segmentation;
	std::queue<node_pointer> nodes_neighbor;

	ini_dijkstra();

	first_node -> distance = (maxflowtype) 0;
	first_node -> dijkstra_parent = NULL;
	add_to_temp_list(first_node);

	first_node_flow = max_flow_graph -> maxflow();
	segment = max_flow_graph->what_segment(first_node->node_identifier);// the segment to which the first node belong
	is_source = (segment == GraphType::SOURCE ? 1:0); // is the first node in source segment

	put_hard_constraint(first_node,is_source,0);
	first_node_flow = max_flow_graph -> maxflow(1);
	remove_hard_constraint(first_node,is_source,0);
	while (1){
		k++;
		current = get_next_best_node();
		if(current == NULL) break;
		current -> is_fixed = 1;
		current -> when_added = k;
		put_hard_constraint(current,is_source,1);

		curr_flow =  (max_flow_graph -> maxflow(1)) - first_node_flow;
		if(std::abs(curr_flow - current -> distance) > 0.001){
			printf("flow changed %d \n",std::abs(curr_flow - current -> distance) );
		};	

		nodes_segmentation.push(current);
		while((int) nodes_segmentation.size()>0){
			node_seg = nodes_segmentation.front();
			nodes_segmentation.pop();
			for (a = node_seg->first; a; a=a->next)
			{
				current_neighboor = a -> head;
				if(current_neighboor -> is_fixed !=1){
					if(in_segmentation(current_neighboor,segment) == 1){
						current_neighboor -> is_fixed = 1;
						current_neighboor -> when_added = k;
						current_neighboor -> dijkstra_parent = node_seg;
						current_neighboor -> distance = curr_flow;
						nodes_segmentation.push(current_neighboor);
						remove_from_temp_list(current_neighboor);
					}
					else{
						nodes_neighbor.push(current_neighboor);
					}
				}
			}
		}
		if (last_node != NULL && last_node->is_fixed) break;
		while((int) nodes_neighbor.size()>0){
			current_neighboor = nodes_neighbor.front();
			nodes_neighbor.pop();
			parent_node = NULL;
			for (a = current_neighboor->first; a; a=a->next)
			{
				if(a -> head -> when_added == k){
					parent_node = a-> head;
					break;
				}
			}
			if(parent_node == NULL){
				printf("Error: no neighboor found");
			}

			put_hard_constraint(parent_node,is_source,current);	
			put_hard_constraint(current_neighboor,is_source,0);
			neigh_flow = max_flow_graph -> maxflow(1) - first_node_flow;
			if(neigh_flow < current_neighboor -> distance){
				current_neighboor -> distance = neigh_flow;
				current_neighboor -> dijkstra_parent = parent_node;
				add_to_temp_list(current_neighboor);
			}
			remove_hard_constraint(current_neighboor,is_source, 0);
			remove_hard_constraint(parent_node,is_source,current);	
		}
		remove_hard_constraint(current,is_source,1);
		if (last_node != NULL && last_node->is_fixed) break;
	}
	tempList.clear();
	if (last_node == NULL)
		return 0.0;
	else{
		put_hard_constraint(last_node,is_source, 1);
		curr_flow = max_flow_graph -> maxflow(1);	
		return curr_flow;
	}


}

template <typename edgeweighttype,typename maxflowtype>
inline int  DijkstraMaxFlow<edgeweighttype,maxflowtype>::in_segmentation(node_pointer node, int segment)
{
	int node_identifier;
	int figure_node;
	typename std::vector<Coordinates>::iterator it;

	node_identifier = node->node_identifier;
	for ( it=figure_coordinates.begin() ; it < figure_coordinates.end(); it++ ){
		figure_node = node_identifier - ((*it).getY()) + nrows*((*it).getX());
		if(figure_node >= 0 && figure_node < node_num){
			if( max_flow_graph -> what_segment(figure_node) != segment){
				return 0;
			}
		}
	}
	return 1;

}

template <typename edgeweighttype, typename maxflowtype> 
inline DijkstraMaxFlow<edgeweighttype,maxflowtype>::Coordinates::Coordinates(int _x,int _y):
x(_x),
y(_y)
{
}

template <typename edgeweighttype, typename maxflowtype> 
inline int DijkstraMaxFlow<edgeweighttype,maxflowtype>::Coordinates::getX()
{
	return x;
}

template <typename edgeweighttype, typename maxflowtype> 
inline int DijkstraMaxFlow<edgeweighttype,maxflowtype>::Coordinates::getY()
{
	return y;
}

template <typename edgeweighttype, typename maxflowtype> 
void DijkstraMaxFlow<edgeweighttype,maxflowtype>::add_figure_coordinates(int _x, int _y)
{
	Coordinates a(_x,_y);
	figure_coordinates.push_back(a);
}

template <typename edgeweighttype, typename maxflowtype> 
void DijkstraMaxFlow<edgeweighttype,maxflowtype>::set_nrows_ncols(int _nrows,int _ncols){
	nrows = _nrows;
	ncols = _ncols;
}
