#include <iostream>
#include <vector>
#include <map>
#include <ctime>
#include <fenv.h>
#include <time.h>
#include <fstream>
#include <chrono>
#include <pthread.h> 
#include <thread>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <math.h>
#include <stdexcept>
#include "v_repExtAMRProject.h"
#include "v_repLib.h"
#include "utils.h"

// edge returned by FindTrajectory
vector<Eigen::Vector3d> edge;

struct VERTEX qRand;

//---------------------------------------------------------------------//

int create_primitives() {
	int map_dim = ground_DIM;
	
	Eigen::Vector3d q_ini;
	q_ini << 0, 0, 0;
	
	fstream fout;
	fout.open("primitives_10x10_res05_8thetas.txt", ios::out);
	
	if (!fout) {
		cout << "Unable to write summary" << endl;
		return -1;
	}
	
	//if (fout) cout << "File created" << endl;
	
	//for (int t_i=0; t_i<discrete_thetas; t_i++) {
	for (float i=0; i<=map_dim; i+=grid_resolution) {
		for (float j=0; j<=map_dim; j+=grid_resolution) {
			for (int t_f=0; t_f<discrete_thetas; t_f++) {
				if (i==0 && j==0) continue;
				cout << "theta_fin = " << thetas[t_f] << ", i = " << i << ", j = " << j << endl;
				
				Eigen::Vector3d q_fin;
				q_fin << i, j, thetas[t_f];
				
				vector<Eigen::Vector3d> plan = steer(q_ini, q_fin);
				
				if (plan.empty()) {

					fout << q_fin(0) << " , " << q_fin(1) << " , " << q_fin(2) << endl;
					
				} else {

					int plan_size = plan.size();
					
					fout << q_fin(0) << " , " << q_fin(1) << " , " << q_fin(2) << " ; ";
					//fout << distance(q_ini, q_fin) << " ; ";
					fout << cost_plan(plan) << " ; ";
					fout << q_ini(0) << " , " << q_ini(1) << " , " << q_ini(2) << " ; ";
					Eigen::Vector3d s;
					for(int k=1; k<plan_size-1; k++){
						s = plan.at(k);
						fout << s(0) << " , " << s(1) << " , " << s(2) << " ; ";
					}
					
					fout << q_fin(0) << " , " << q_fin(1) << " , " << q_fin(2) << " ; " <<  endl;
				}
				
				q_fin << -i, -j, thetas[t_f];
				
				plan = steer(q_ini, q_fin);
				
				if (plan.empty()) {

					fout << q_fin(0) << " , " << q_fin(1) << " , " << q_fin(2) << endl;
					
				} else {

					int plan_size = plan.size();
					
					fout << q_fin(0) << " , " << q_fin(1) << " , " << q_fin(2) << " ; ";
					//fout << distance(q_ini, q_fin) << " ; ";
					fout << cost_plan(plan) << " ; ";
					fout << q_ini(0) << " , " << q_ini(1) << " , " << q_ini(2) << " ; ";
					Eigen::Vector3d s;
					for(int k=1; k<plan_size-1; k++){
						s = plan.at(k);
						fout << s(0) << " , " << s(1) << " , " << s(2) << " ; ";
					}
					
					fout << q_fin(0) << " , " << q_fin(1) << " , " << q_fin(2) << " ; " <<  endl;
				}
			}
		}
	}
	
	fout.close();			
	
	return 0;
}

//------------------------------------------------------------------------//

int is_in_tree(Eigen::Vector3d qRand) {
	int tree_size = tree.size();
	
	for (int i=0; i<tree_size; i++) {
		if (tree[i].z.head(3) == qRand) return i;
	}
	
	return -1; //index which is not in the tree
}

//-----------------------------------------------------------------------//

int FindTrajectory(Eigen::VectorXd zNear, Eigen::Vector3d zRand, float c_best) {
	
	time_t start, end;
	time(&start);
	
	float path_cost;
	int success = -1;
	
	// Translate zNear to [0 0 0] and zRand accordingly	
	Eigen::Vector3d zFin;
	zFin << zRand[0]-zNear[0], zRand[1]-zNear[1], zRand[2];
	
	float x, y, theta;
	
	x = zFin[0] * cos(-zNear[2]) - zFin[1] * sin(-zNear[2]);
	y = zFin[0] * sin(-zNear[2]) + zFin[1] * cos(-zNear[2]);
	theta = zFin[2] - zNear[2];
	if (theta < -M_PI || theta > 3*M_PI_4) {
		if (theta < 0) theta += 2*M_PI;
		else theta -= 2*M_PI;
	}
	
	zFin[0] = x;
	zFin[1] = y;
	zFin[2] = theta;
	
	
	if (((x < 0 && y >= 0) || (x >= 0 && y < 0))) { // && theta == zNear[2]) {
		zFin[1] *= -1;
		zFin[2] *= -1;
		if (theta == M_PI) zFin[2] = -M_PI;
	}
	
	
	edge.clear();	
	
	// We look for a primitive in the primitives database
	fstream fin;
	fin.open("primitives_10x10_res05_8thetas.txt", ios::in);
	
	if (!fin) {
		
		cout << "Unable to write summary" << endl;
		return success;
		
	} else {
			
		string line;
		while (getline(fin, line)) {

			string fin = line.substr(0, line.find(" ; "));
			char* str = new char [fin.length()+1];
			strcpy(str, fin.c_str());
			
			char* p = strtok(str, " , ");
			
			Eigen::Vector3d qFin;
			qFin << 0.0,0.0,0.0;
			int i = 0;
			
			while (p != NULL) {
				if (i != 2 && grid_resolution==1) qFin[i] = stoi(p);
				else qFin[i] = stof(p);
				i++;
				p = strtok(NULL, " , ");
			}
			
			
			// if the primitive doesn't bring to zFin we go to the next line
			if (abs(zFin[0] - qFin[0]) > 10e-2 || abs(zFin[1] - qFin[1]) > 10e-2 || abs(zFin[2] - qFin[2]) > 10e-2) continue;
			
			// This is the rest of the line in the file
			if (line.find(" ; ") == -1) break; //if there is no path to qFin
			
			
			string path_str = line.substr(line.find(" ; "), line.length());
			path_str = path_str.substr(path_str.find(" ; ")+2, path_str.length());
			
			// cost of the path from [0 0 0] to zFin
			string cost = path_str.substr(0, path_str.find(" ; "));
			path_cost = stof(cost) + zNear[5];
			
			if (path_cost >= c_best) break;
			
			// We are left with the path from [0 0 0] to zFin
			path_str = path_str.substr(path_str.find(" ; ")+2, path_str.length());
			string elem = path_str.substr(0, path_str.find(" ; "));
			
			while (elem != " ") {
				char* pp = new char [elem.length()+1];
				strcpy(pp, elem.c_str());
			
				char* node = strtok(pp, " , ");
				
				Eigen::Vector3d q;
				q << 0,0,0;
				i = 0;
				while(node != NULL) {
					
					try {
						q[i] = stof(node);
					} catch (const out_of_range& oor) {
						q[i] = 0.0;
					}
					
					i++;
					node = strtok(NULL, " , ");
				}
				
				float x_f, y_f, theta_f;
				
				// FLIP BACK
				if (((x < 0 && y >= 0) || (x >= 0 && y < 0))) { // && theta == zNear[2]) {
					y_f = -q[1];
					theta_f = -q[2];
				} else {
					y_f = q[1];
					theta_f = q[2];
				}					
				
				
				// ROTATE BACK
				x_f = q[0] * cos(zNear[2]) - y_f * sin(zNear[2]);
				y_f = q[0] * sin(zNear[2]) + y_f * cos(zNear[2]);
				theta_f = theta_f + zNear[2];

				//TRANSLATE BACK ALL THE CONFIGURATIONS IN THE PATH
				q[0] = x_f + zNear[0];
				q[1] = y_f + zNear[1];
				q[2] = theta_f;
				
				
				edge.push_back(q);
				
				path_str = path_str.substr(path_str.find(" ; ")+2, path_str.length());
				elem = path_str.substr(0, path_str.find(" ; "));
			}
			
			// We create a fake node in order to make the cost the last element of the edge to be returned
			Eigen::Vector3d cost_node;
			cost_node << path_cost, -1, -1;
			edge.push_back(cost_node);
						
			// we have found the primitive so we do not need to check the others
			success = 0;
			break;
		}
	}
	
	fin.close();
	
	time(&end);
	double duration_TPBVP = double(end-start);
	tpbvp_duration += duration_TPBVP;

	return success;
}

//-----------------------------------------------------------------------//

int extend(vector<Eigen::VectorXd> ZNear, Eigen::VectorXd zRand) {
	float c_best = 20000000;
	float path_cost;
	Eigen::Vector3d q_best;
	int idx_best = -1;
	vector<Eigen::Vector3d> e_best;
	int found = -1;
	
	
	for (int i=0; i<ZNear.size(); i++) {
		if (abs(ZNear[i][0]-zRand[0]) > ground_DIM || abs(ZNear[i][1]-zRand[1]) > ground_DIM) continue;
		
		int res = FindTrajectory(ZNear[i], zRand.head(3), c_best);
		
		if (res == -1) continue;
		
		path_cost = edge.back()[0];
		edge.pop_back();
		
		if (check_collisions(edge)) {
			found = 0;
			q_best = ZNear[i].head(3);
			c_best = path_cost;
			e_best = edge;
			idx_best = ZNear[i][3];
		}
	}
	
	
	Eigen::VectorXd z(6);
	z << zRand.head(3), tree.size(), idx_best, c_best;
		
	struct VERTEX nodeRand;
	nodeRand.z = z;
	nodeRand.parent_to_node = e_best;
	
	qRand = nodeRand;
	
	return found;
}

//-----------------------------------------------------------------------//

void rewire_mp(Eigen::VectorXd zRand, vector<Eigen::VectorXd> ZNear) {
	int n_neighbors = ZNear.size();
	
	for (int i=0; i<n_neighbors; i++) {
		if (abs(ZNear[i][0]-zRand[0]) > ground_DIM || abs(ZNear[i][1]-zRand[1]) > ground_DIM) continue;
		
		int res = FindTrajectory(zRand, ZNear[i].head(3), ZNear[i][5]);
		
		if (res == -1) continue;
		
		else {
			float new_cost = edge.back()[0];
			edge.pop_back();
			
			if (check_collisions(edge)) {
				cout << "4 - Rewired" << endl;
				tree[ZNear[i][3]].parent_to_node = edge;
				tree[ZNear[i][3]].z[4] = zRand[3];
				tree[ZNear[i][3]].z[5] = new_cost;
			}
		}
	}
}

//------------------------------------------------------------------------//

void RRTstar_mp(Eigen::Vector3d q_i, Eigen::Vector3d q_f){  //returns the same type we defined for the tree
	time_t start;
	time(&start);
		
	cout << "RRT* with motion primitives:" << endl;
	
	tree.clear();
	path.clear();
		
	vector<uniform_real_distribution<float_t>> distributions = initialize_generators();
	
	// Initialization of the tree : insert the initial configuration as root
	Eigen::VectorXd z_root(6);
	z_root << q_i(0), q_i(1), q_i(2), 0, -1, 0;  //-1 xk the root has no parent node
	
	struct VERTEX root;
	root.z = z_root;
	
	tree.push_back(root);
	
	double dist = 0.0;
	goal_iters = -1;
	goal_updates = 0;
	int in_tree = -1;
	
	int i;
	for(i = 0; i < maxIter; i++) {
		cout << "\nITERATION = " << i << endl;
		Eigen::VectorXd zRand = sampleRand(distributions, q_f);
		
		vector<Eigen::VectorXd> ZNear = NearByVertices(zRand.head(3));
		
		int res = extend(ZNear, zRand);
		if (res == -1) continue;
		
		in_tree = is_in_tree(zRand.head(3));
		
		if (in_tree == -1 && qRand.z[4] != -1) {
			cout << "1 - Node added" << endl;
			tree.push_back(qRand);
			rewire_mp(qRand.z, ZNear);
		} 
		
		else if (in_tree != -1) {
			cout << "2 - Node already in the tree" << endl;
			Eigen::VectorXd qPrev = tree[in_tree].z; //parent of qRand already in the tree
			Eigen::VectorXd qBest = tree[qRand.z[4]].z; //parent of qRand randomly sampled
			
			if (qPrev.head(3) != qBest.head(3) && qBest[5] < qPrev[5]) {
				cout << "3 - Node updated" << endl;
				tree[in_tree].z[4] = qRand.z[4];
				tree[in_tree].z[5] = qRand.z[5];
				tree[in_tree].parent_to_node = qRand.parent_to_node;
				
				rewire_mp(tree[in_tree].z, ZNear);
			}
		}
		
		dist = distance(qRand.z.head(3), q_f, 0.0);

		if (dist < goal_tolerance) update_goal_info_mp(qRand.z, i, dist, q_i, q_f, start);
	}
	
	planning_info(start, i, q_i, q_f);
	
}


