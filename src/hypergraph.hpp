#include <boost/graph/adjacency_list.hpp>
#include <fstream>
#include <limits>
#include <vector>
#include <deque>
#include <map>

using namespace std;

struct vertex_info {
	string name;
	bool stop_vertex; //if==1 it's a stop vertex
	string lat, lon; //string to avoid rounding problem
};

struct edge_info {
	string name;
	double weight;  //time, cost, frequency or incomfort
}; 

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
		vertex_info, edge_info> hypergraph;

typedef hypergraph::vertex_descriptor vertex_descriptor;
typedef hypergraph::edge_descriptor edge_descriptor;
typedef hypergraph::vertex_iterator vertex_iterator;
typedef hypergraph::edge_iterator edge_iterator;
typedef hypergraph::out_edge_iterator out_edge_iterator;
typedef hypergraph::in_edge_iterator in_edge_iterator;

class Hypergraph {
  
  private: 

	hypergraph h;
	static const unsigned int theta = 30;

	vector<edge_descriptor> get_bstar(vertex_descriptor v) {
		vector<edge_descriptor> bstar;
		in_edge_iterator ie_i, ie_end;
		for(boost::tie(ie_i,ie_end)=in_edges(v,h); ie_i!=ie_end; ++ie_i) {
			bstar.push_back(*ie_i);
		}
		return bstar;
	}

  	/*--------------------------------------------- INPUT ---------------------------------------------*/
	bool load_hypergraph(const char* filename) {
		ifstream ifs(filename);
		if(!ifs) {
			cout<<"Failed to open file "<<filename<<endl;
			return false;
		}
		string garbage;
		unsigned int n;
		unsigned int id;		
		bool stop_vertex;
		string lat, lon;
		string vertex_name;
		ifs >> n;
		getline(ifs, garbage);
		h = hypergraph(n);
		for(int line=0;line<n;line++) {
			ifs >> id >> stop_vertex >> lat >> lon;
			ifs.ignore(1,' '); //elimina lo spazio prima del nome
			getline(ifs, vertex_name);
			h[id].stop_vertex = stop_vertex;
			h[id].lat = lat;
			h[id].lon = lon;
			h[id].name = vertex_name;
		}
		getline(ifs, garbage);
		unsigned int m, i, j;
		double weight;
		string edge_name;
		ifs >> m;
		getline(ifs, garbage);
		edge_descriptor e;
		for(int line=0;line<m;line++) {
			ifs >> i >> j >> weight;
			ifs.ignore(1,' '); //elimina lo spazio prima del nome
			getline(ifs, edge_name);
			e = boost::add_edge(i,j,h).first;
			h[e].weight = weight;
			h[e].name = edge_name;
		}
		return true;
	}
	/*-------------------------------------------------------------------------------------------------*/

	/*--------------------------------------------- OUTPUT --------------------------------------------*/
	void export_shp(const char* filename, const char* alg, vertex_descriptor s, vertex_descriptor t, double expected_cost[], 
			map<vertex_descriptor, double> comb_freq, vector<vertex_descriptor> successors[]) {
		char* filename_shp = (char*)malloc(strlen(filename)+5+strlen(alg)+1+h[s].name.size()+1+h[t].name.size()+1);
		char* s_name = (char*)malloc(h[s].name.size());
		char* t_name = (char*)malloc(h[t].name.size());
		strcpy(filename_shp, filename);		
		strcpy(s_name, h[s].name.c_str());		
		strcpy(t_name, h[t].name.c_str());
		strcat(filename_shp, "_shp_");
		strcat(filename_shp, alg);
		strcat(filename_shp, "_");
		strcat(filename_shp, s_name);
		strcat(filename_shp, "_");
		strcat(filename_shp, t_name);
		ofstream ofs(filename_shp);
		//nodi		
		unsigned int n = num_vertices(h);
		bool in_shp[n];
		for(int i=0; i<n; i++)
			in_shp[i] = 0;	
		ofs<<"n                              "<<endl; //segnaposto
		int shp_v_counter = 0;
		print_successors(ofs, shp_v_counter, s, expected_cost, comb_freq, successors, in_shp);
		ofs.seekp(0, ios_base::beg);
		ofs<<shp_v_counter;
		ofs.seekp(0, ios_base::end);
		ofs<<endl;
		//archi
		int edges_start_pos = ofs.tellp();
		ofs<<"m                              "<<endl; //segnaposto
		int shp_e_counter = 0;		
		edge_iterator e_i, e_end;
		for(boost::tie(e_i,e_end)=edges(h); e_i!=e_end; ++e_i) {
			vertex_descriptor i, j;
			i = source(*e_i,h);		
			j = target(*e_i,h);
			if(in_shp[i] && in_shp[j]) {
				string edge_name = h[*e_i].name;				
				double edge_cost;
				double edge_prob;
				if(h[i].stop_vertex) {
					edge_cost = theta / comb_freq[i];
					edge_prob = h[*e_i].weight / comb_freq[i];
				}
				else {
					edge_cost = h[*e_i].weight;
					edge_prob = 1;
				}
				ofs<<i<<" "<<j<<" "<<edge_cost<<" "<<edge_prob<<" "<<edge_name<<endl;
				shp_e_counter++;
			}
		}
		ofs.seekp(edges_start_pos);
		ofs<<shp_e_counter;
		ofs.seekp(0, ios_base::end);
		ofs<<endl;
		ofs.close();
	}

	void print_successors(ofstream& ofs, int& shp_v_counter, vertex_descriptor v, double expected_cost[],
				map<vertex_descriptor, double> comb_freq, vector<vertex_descriptor> successors[], bool in_shp[]) {
		if(!in_shp[v]) {
			in_shp[v] = 1;
			shp_v_counter++;
			ofs<<v<<" "<<h[v].stop_vertex<<" "<<h[v].lat<<" "<<h[v].lon<<" "<<expected_cost[v]<<" "<<comb_freq[v]<<" "<<h[v].name<<endl;
			for(int i=0; i<successors[v].size(); i++) {
				vertex_descriptor suc = successors[v][i];
				print_successors(ofs, shp_v_counter, suc, expected_cost, comb_freq, successors, in_shp);
			}
		}
	}
	/*-------------------------------------------------------------------------------------------------*/

  public:

	/*------------------------------------------ PRIORITY SHT -----------------------------------------*/
	void priority_sht(vertex_descriptor s, vertex_descriptor t, const char* filename) {
		if(!load_hypergraph(filename)) return;		
		cout<<"Algoritmo SHT - priority su : "<<filename<<endl; //	
		unsigned int iteration_counter = 0; //
		unsigned int n = num_vertices(h);
		//expected_cost, successors and comb_freq as exterior property map		
		double expected_cost[n];
		for(int k=0; k<n; ++k) 
			expected_cost[k] = std::numeric_limits<double>::infinity();
		map<vertex_descriptor, double> comb_freq;
		for(int k=0; k<n; ++k) {
			if(h[k].stop_vertex)
				comb_freq[k] = 0;
		}
		vector<vertex_descriptor> successors[n];
		deque<vertex_descriptor> Q;		
		expected_cost[t] = 0;
		Q.push_back(t);
		while(!Q.empty()) {	
			vertex_descriptor v = Q.front();
			Q.pop_front();
			iteration_counter++; //
			cout << iteration_counter << ".  Q --> " << h[v].name << " , Q <-- "; //
			vector<edge_descriptor> bstar = get_bstar(v);
			for(int k=0; k<bstar.size(); ++k) {
				vertex_descriptor u = source(bstar[k], h);
				if(h[u].stop_vertex && expected_cost[u] > expected_cost[v]) {
					if(comb_freq[u]==0) {
						comb_freq[u] = h[bstar[k]].weight;
						expected_cost[u] = (theta/comb_freq[u])+expected_cost[v];
					}
					else {
						comb_freq[u] += h[bstar[k]].weight;
						expected_cost[u] -= (expected_cost[u]-expected_cost[v])*h[bstar[k]].weight/comb_freq[u];
					}
					successors[u].push_back(v);
					if(u!=s) {
						priority_insertion(u,Q,expected_cost);
						cout << h[u].name << "(" << expected_cost[u] << ")  "; //
					}
				}
				else if(!h[u].stop_vertex && expected_cost[u] > expected_cost[v] + h[bstar[k]].weight) {
					expected_cost[u] = expected_cost[v] + h[bstar[k]].weight;
					if(successors[u].empty())
						successors[u].push_back(v);
					else
						successors[u][0] = v;
					if(u!=s) { 
						priority_insertion(u,Q,expected_cost);
						cout << h[u].name << "(" << expected_cost[u] << ")  "; //
					}
				}
			}
			cout<<endl; //
		}
		cout<<"Expected cost"<<endl;
		for(int k=0;k<n;++k) {
			cout<<h[k].name<<" - "<<expected_cost[k]<<" - successore/i: ";
			for(int suc=0; suc<successors[k].size(); ++suc)
				cout<<h[successors[k][suc]].name<<" ";
			cout<<endl;
		}	
		cout<<"Numero estrazioni dalla coda: "<<iteration_counter<<endl;
		export_shp(filename, "priority", s, t, expected_cost, comb_freq, successors);
	}

	/*inserisce l'elemento in posizione corretta (se è già presente ha sicuramente un costo maggiore), quindi continuo a esplorare il resto del vettore e lo elimino se lo trovo. */
	void priority_insertion(vertex_descriptor u, deque<vertex_descriptor>& Q, double expected_cost[]) {
		int pos = 0;
		while(pos < Q.size() && expected_cost[u] > expected_cost[Q[pos]])
			pos++;				
		Q.insert(Q.begin()+pos, u);
		pos++;
		while(pos<Q.size()) {
			if(u==Q[pos]) Q.erase(Q.begin()+pos);
			pos++;
		}			
	}
	/*-------------------------------------------------------------------------------------------------*/

	/*---------------------------------------------- SHT ----------------------------------------------*/
	void sht(vertex_descriptor s, vertex_descriptor t, const char* filename) {
		if(!load_hypergraph(filename)) return;		
		cout<<"Algoritmo SHT su : "<<filename<<endl; //	
		unsigned int iteration_counter = 0; //
		unsigned int n = num_vertices(h);
		//expected_cost, successors and comb_freq as exterior property map		
		double expected_cost[n];
		for(int k=0; k<n; ++k) 
			expected_cost[k] = std::numeric_limits<double>::infinity();
		map<vertex_descriptor, double> comb_freq;
		for(int k=0; k<n; ++k) {
			if(h[k].stop_vertex)
				comb_freq[k] = 0;
		}
		vector<vertex_descriptor> successors[n];
		deque<vertex_descriptor> Q;		
		expected_cost[t] = 0;
		Q.push_back(t);
		while(!Q.empty()) {
			vertex_descriptor v = Q.front();
			Q.pop_front();
			iteration_counter++; //
			cout << iteration_counter << ".  Q --> " << h[v].name << " , Q <-- "; //
			vector<edge_descriptor> bstar = get_bstar(v);
			for(int k=0; k<bstar.size(); ++k) {
				vertex_descriptor u = source(bstar[k], h);
				if(h[u].stop_vertex && expected_cost[u]>expected_cost[v]) {
					expected_cost[u] = attractive_set(u, expected_cost, comb_freq[u], successors[u]);
					if(not_in_queue(u, Q) && u!=s) {
						Q.push_back(u);
						cout << h[u].name << "(" << expected_cost[u] << ")  "; //
					}
				}
				else if(!h[u].stop_vertex && expected_cost[u]>expected_cost[v]+h[bstar[k]].weight) {
					expected_cost[u] = expected_cost[v] + h[bstar[k]].weight;
					if(successors[u].empty())
						successors[u].push_back(v);
					else
						successors[u][0] = v;
					if(not_in_queue(u, Q) && u!=s) {
						Q.push_back(u);
						cout << h[u].name << "(" << expected_cost[u] << ")  ";
					}
				} 
			}
			cout<<endl; //
		}	
		cout<<"Expected cost"<<endl;
		for(int k=0;k<n;++k) {
			cout<<h[k].name<<" - "<<expected_cost[k]<<" - successore/i: ";
			for(int suc=0; suc<successors[k].size(); ++suc)
				cout<<h[successors[k][suc]].name<<" ";
			cout<<endl;
		}	
		cout<<"Numero estrazioni dalla coda: "<<iteration_counter<<endl;
		export_shp(filename, "no_priority", s, t, expected_cost, comb_freq, successors);
	}

	bool not_in_queue(vertex_descriptor u, deque<vertex_descriptor> Q) {
		bool not_in_queue = true;
		int pos = 0;
		while(pos<Q.size() && not_in_queue) {
			if(Q[pos]==u)
				not_in_queue = false;
			pos++;
		}
		return not_in_queue;
	}
	/*-------------------------------------------------------------------------------------------------*/

	/*----------------------------------------- ATTRACTIVE SET ----------------------------------------*/
	double attractive_set(vertex_descriptor u, double expected_cost[], double& comb_freq, vector<vertex_descriptor>& successors) {
		vector<edge_descriptor> sorted_fstar = get_sorted_fstar(u, expected_cost);
		successors.clear();
		edge_descriptor first_edge = sorted_fstar[0];
		successors.push_back(target(first_edge,h));
		comb_freq = h[first_edge].weight;
		double expected_cost_u = (theta/comb_freq) + expected_cost[target(first_edge,h)];		
		int index = 1;
		bool stop_criterion = false;
		while(index<sorted_fstar.size() && stop_criterion==false) {
			edge_descriptor next_edge = sorted_fstar[index];
			if(expected_cost[target(next_edge,h)] < expected_cost_u) {
				successors.push_back(target(next_edge,h));
				comb_freq += h[next_edge].weight;
				expected_cost_u -= ((expected_cost_u - expected_cost[target(next_edge,h)]) * h[next_edge].weight)/comb_freq;
				index++;
			}			
			else stop_criterion = true;
		}
		return expected_cost_u;
	}
		
	vector<edge_descriptor> get_sorted_fstar(vertex_descriptor i, double expected_cost[]) {
		vector<edge_descriptor> fstar;
		out_edge_iterator oe_i, oe_end;
		for(boost::tie(oe_i,oe_end)=out_edges(i,h); oe_i!=oe_end; ++oe_i) {
			int pos = 0;
			while(pos<fstar.size() && expected_cost[target(*oe_i,h)] > expected_cost[target(fstar[pos],h)])
				pos++;
			fstar.insert(fstar.begin()+pos, *oe_i);
		}
		return fstar;
	}	
	/*-------------------------------------------------------------------------------------------------*/

};
