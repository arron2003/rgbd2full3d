#include <queue>
#include <vector>
#include <iostream>

#include "dijkstraMaxFlow.h"
#include "instances_dijkstra.inc"

template <typename edgeweighttype, typename maxflowtype>
void DijkstraMaxFlow<edgeweighttype,maxflowtype>::percolateDown(node_pointer node)
{
	node_pointer child;
	int i,j,s;

	i = node -> position_in_tempList;
	s = (int) tempList.size();

	if( 2*i+1 < s){
		if(2*i+2 < s){
			child = (tempList[2*i+1] -> distance < tempList[2*i+2] -> distance)? tempList[2*i+1] : tempList[2*i+2];
		}
		else{
			child = tempList[2*i+1];
		}
		while( node -> distance > child -> distance){
			i = child -> position_in_tempList;
			j = node -> position_in_tempList;
			tempList[i] = node;
			tempList[j] = child;
			node -> position_in_tempList = i;
			child -> position_in_tempList = j;
			if( (2*i+1) < s)
				if((2*i+2) < s){
					child = (tempList[2*i+1] -> distance < tempList[2*i+2] -> distance)? tempList[2*i+1] : tempList[2*i+2];
				}
				else{
					child = tempList[2*i+1];
				}
			else 
				break;
		}
	}
}

template <typename edgeweighttype, typename maxflowtype>
void DijkstraMaxFlow<edgeweighttype,maxflowtype>::percolateUp(node_pointer node)
{
	node_pointer parent;
	int i,j;

	i = node->position_in_tempList;
	
	if((i-1)/2 >= 0){
		parent = tempList[(i-1)/2];

		while( node->distance < parent->distance){
			i = parent->position_in_tempList;
			j = node->position_in_tempList;
			tempList[i] = node;
			tempList[j] = parent;
			node->position_in_tempList = i;
			parent->position_in_tempList = j;
			if( i > 0){
				parent = tempList[(i-1)/2];
			}
			else
				break;
		}
	}

}


template <typename edgeweighttype, typename maxflowtype>
void DijkstraMaxFlow<edgeweighttype,maxflowtype>::add_to_temp_list(node_pointer node)
{
	if(node -> position_in_tempList < 0){
		tempList.push_back(node);
		node -> position_in_tempList = (int) tempList.size()-1;
		percolateUp(node);
	}
	else{
		percolateUp(node);
	}
}

template <typename edgeweighttype, typename maxflowtype>
typename DijkstraMaxFlow<edgeweighttype,maxflowtype>::node_pointer DijkstraMaxFlow<edgeweighttype,maxflowtype>::get_next_best_node()
{
	node_pointer best_node, aux;
	if( tempList.size()>0 ){
		best_node = tempList.front();
		aux = tempList.back();
		aux -> position_in_tempList = 0;
		tempList[0] = aux;
		//tempList[tempList.size()-1] = best_node;
		tempList.pop_back();
		percolateDown(aux);
		best_node -> position_in_tempList = -2;
		return best_node;
	}
	else{
		return NULL;
	}
}


template <typename edgeweighttype, typename maxflowtype>
void DijkstraMaxFlow<edgeweighttype,maxflowtype>::remove_from_temp_list(node_pointer node)
{
	int i;
	node_pointer aux;

	if(node-> position_in_tempList >= 0){
		i = node->position_in_tempList;
		aux = tempList.back();
		aux -> position_in_tempList = i;
		tempList[i] = aux;
		tempList.pop_back();
		percolateDown(aux);
		node -> position_in_tempList = -2;
	}
}
