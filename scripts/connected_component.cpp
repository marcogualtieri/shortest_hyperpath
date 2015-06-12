#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/strong_components.hpp>
#include <iostream>
#include <fstream>

using namespace std;

struct node_info {
	int id, stop_vertex;
	string lat, lon, name;  //saved as string to avoid rounding problem
};

struct way_info {
	int i, j;
	double weight;
	string name;
}; 

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
		node_info, way_info> graph;

int main(int argc, char* argv[]) {
	if(argc <= 1)
		cout<<("usage: $./connected_component filename");
	char* filename = argv[1];
	//load nodes set	
	ifstream ifs(filename);
	if(!ifs) {
    		cout<<endl<<"Failed to open file "<<filename<<endl;
    		return 1;
  	}
	int id, stop_vertex;
	string lat, lon, n_name;
	string garbage;
	int n;
	ifs >> n;
	graph g(n);
	getline(ifs, garbage);
	for(int n_counter=0;n_counter<n;n_counter++) {
		ifs >> id >> stop_vertex >> lat >> lon;
		ifs.ignore(1,' '); //elimina lo spazio prima del nome
		getline(ifs, n_name);
		g[n_counter].id = id;
		g[n_counter].stop_vertex = stop_vertex;
		g[n_counter].lat = lat; 
		g[n_counter].lon = lon;
		g[n_counter].name = n_name;
	}
	//load edges set
	int i, j;
	double weight;
	string w_name;
	int m;
	ifs >> m;
	getline(ifs, garbage);
	graph::edge_descriptor e;
	for(int w_counter=0;w_counter<m;w_counter++) {
		ifs >> i >> j >> weight;
		ifs.ignore(1,' '); //elimina lo spazio prima del nome
		getline(ifs, w_name);
		e = boost::add_edge(i,j,g).first;
		g[e].i = i;
		g[e].j = j;
		g[e].weight = weight;
		g[e].name = w_name;
	}
	//find strong connected components
	vector<int> component(num_vertices(g));
	int num = boost::strong_components(g, &component[0]);
	int component_cardinality[num];
	for(int k=0;k<num;++k) 
		component_cardinality[k] = 0;
	vector<int>::size_type index;
	for(index = 0; index != component.size(); ++index) {
		//cout << "Vertex " << index <<" is in component " << component[index] << endl;
		component_cardinality[component[index]]++;
	}
	//find largest component
	int max = 0;
	int largest_comp_index;
	for(int k=0;k<num;++k) {
		if(component_cardinality[k] > max) {
			max = component_cardinality[k];
			largest_comp_index = k;
		}
	}
	//save nodes of connected component
	char* new_filename = (char*)malloc(strlen(filename)+8+1);
	strcpy(new_filename, filename);
	strcat(new_filename, "_traffic");	
	ofstream out(new_filename);
	out << "n                              " << endl; //segnaposto
	int conn_node_counter = 0;
	for(int k=0; k<n; ++k) {
		if(component[k] == largest_comp_index) {
			g[k].id = conn_node_counter;
			out << g[k].id << " " << g[k].stop_vertex << " " << g[k].lat << " " << g[k].lon << " " << g[k].name << endl; 
			conn_node_counter++;
		}
	}
	out.seekp(0, ios_base::beg);
	out << conn_node_counter;
	out.seekp(0, ios_base::end);
	out<<endl;
	//find edges of largest connected component
	vector<graph::edge_descriptor> largest_comp_edges;	
	graph::edge_iterator e_i, e_end;
	for(boost::tie(e_i,e_end)=edges(g); e_i!=e_end; ++e_i) {
		graph::vertex_descriptor s, t;
		s = source(*e_i,g);		
		t = target(*e_i,g);
		if(component[s]==largest_comp_index && component[t]==largest_comp_index) {
			largest_comp_edges.push_back(*e_i);
		}
	}
	//save edges
	out << largest_comp_edges.size() << endl;
	graph::edge_descriptor edge;
	for(int k=0;k<largest_comp_edges.size();++k) {
		edge = largest_comp_edges[k];
		out << g[g[edge].i].id << " " << g[g[edge].j].id << " " << g[edge].weight << " " << g[edge].name << endl;
	}
	out.close();
	cout<<"Traffic graph saved in "<<new_filename<<endl;
	return 0;
}
