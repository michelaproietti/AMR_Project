#include <unistd.h>
#include "v_repExtAMRProject.h"
#include "v_repLib.h"
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
#include <random>
#include "rrt_star_primitives.h"

using namespace std;

//-----------------------------------------------------------------------------//


// given a state random_point and the tree, returns the nearet node in terms 
// of the distance fne.
Eigen::VectorXd NearestNeighbor(Eigen::VectorXd zRand) {
	int tree_size = tree.size();
	
	//initialize the nearest node with the root of the tree
	//out << "before nearest" << endl;
	Eigen::VectorXd nearest_point = tree.front().z;
	double min_distance = distance(zRand, nearest_point);  

	//loop over the tree to find the actual nearest node
	for(int i = 1; i < tree_size; i++) {
		Eigen::VectorXd point = tree[i].z;
		double new_distance = distance(zRand, point);
		if (new_distance < min_distance) {
			min_distance = new_distance;
			nearest_point = point;
		}
	}

	return nearest_point;
}

//------------------------------------------------------------------------------------

// zParent is the chosen parent for zNew
// I add a new node to the tree, with the configuration of zNew, and the index od zParent as index.
Eigen::VectorXd addNewNode(Eigen::VectorXd zParent, Eigen::Vector3d zNew, double cmin) {
	Eigen::VectorXd qNew(6);
	//double costNew = zParent[5] + distance(zParent, zNew);
	qNew << zNew[0], zNew[1], zNew[2], tree.size(), zParent[3], cmin;
	
	struct VERTEX new_node;
	new_node.z = qNew;
	new_node.parent_to_node = parent_plan;
	tree.push_back(new_node);
	
	return qNew;
}


//------------------------------------------------------------------------------------

// I take into account ZNear, that contains all the nodes N already in the tree s.t. distance(n, zNew)<radius
// I compare the costs to go from each N to zNew and return the node in ZNear that minimizes that cost.
Eigen::VectorXd ChooseParent(vector<Eigen::VectorXd> ZNear, Eigen::VectorXd zNearest, Eigen::Vector3d zNew, vector<Eigen::Vector3d> plan) {
	Eigen::VectorXd zmin = zNearest;
	double cmin = zNearest[5] + cost_plan(plan);//distance(zNearest, zNew); //subtitute distance with actual cost
	parent_plan = plan;

	for(int i = 0; i < ZNear.size(); i++) {
		
		vector<Eigen::Vector3d>  plan = steer(ZNear[i].head(3), zNew);
		
		double cost_prime = ZNear[i][5] + cost_plan(plan);//distance(ZNear[i], zNew);
		if(check_collisions(plan) && cost_prime < cmin){
			zmin = ZNear[i];
			cmin = cost_prime;
			parent_plan = plan;		//if we want to save edges between nodes
		}
	}
	
	Eigen::VectorXd qNew = addNewNode(zmin, zNew, cmin);
	
	return qNew;
}

//-------------------------------------------------------------------------------------
			

// we added a new node to the tree. Now we have to verify if an already existing node 
//can be reached from this newly added node with a smaller cost.

void rewire(vector<Eigen::VectorXd> ZNear, Eigen::VectorXd zMin, Eigen::VectorXd zNew){
	for(int i = 0; i < ZNear.size(); i++) {
		if(ZNear[i] == zMin || ZNear[i]==tree.front().z) continue;

		vector<Eigen::Vector3d> plan = steer(zNew.head(3), ZNear[i].head(3));
		
		if(check_collisions(plan) && zNew[5]+cost_plan(plan) < ZNear[i][5]) {   
			cout << "REWIRINGGGGGGGGGGGGGGGGGGGGGGGGG" << endl;
			tree[ZNear[i][3]].parent_to_node = plan;
			tree[ZNear[i][3]].z[4] = zNew[3];
			tree[ZNear[i][3]].z[5] = zNew[5]+cost_plan(plan);
		}
	}
}


//--------------------------------MAIN RRT*--------------------------------------------

void RRTstar(Eigen::Vector3d q_i, Eigen::Vector3d q_f){  //returns the same type we defined for the tree
	
	time_t start;
	time(&start);
		
	cout << "RRT*:" << endl;
	
	vector<uniform_real_distribution<float_t>> distributions = initialize_generators();
	
	tree.clear();
	path.clear();
	
	// Initialization of the tree : insert the initial configuration as root
	Eigen::VectorXd z_root(6);
	z_root << q_i(0), q_i(1), q_i(2), 0, -1, 0;  //-1 xk the root has no parent node
	
	struct VERTEX root;
	root.z = z_root;
	
	tree.push_back(root);
	
	double dist = 0.0;
	goal_iters = -1;
	goal_updates = 0;
	
	double dist_x = 0.0;
	double dist_y = 0.0;
	double dist_theta = 0.0;
	
	int i;
	for(i = 0; i < maxIter; i++) {
		cout << "ITERATION: " << i << endl;
		
		Eigen::VectorXd zRand = sampleRand(distributions, q_f);
		Eigen::VectorXd zNearest = NearestNeighbor(zRand);

		
		vector<Eigen::Vector3d> plan = steer(zNearest.head(3), zRand.head(3));
		if (plan.empty()) {
			cout << "_______________________________________________EXCEPTION AVOIDED" << endl;
			i--;
			continue;
		}
		
		Eigen::Vector3d zNew = plan.back();   //xNew is the last element of the plan computed through the TPBVP
				
		
		if (check_collisions(plan)){
			vector<Eigen::VectorXd> ZNear = NearByVertices(zNew);
			Eigen::VectorXd qNew = ChooseParent(ZNear, zNearest, zNew, plan);
			Eigen::VectorXd zMin = tree[qNew[4]].z;	//parent of qNew
			rewire(ZNear, zMin, qNew);
					
			//dist = distance(zNew, q_f);
			dist_x = abs(zNew[0] - q_f[0]);
			dist_y = abs(zNew[1] - q_f[1]);
			//dist_theta = abs(zNew[2] - q_f[2]);
			
			if (dist_x <= 0.5 && dist_y <= 0.5) update_goal_info(qNew, i, dist_x, dist_y, dist_theta, q_i, q_f, start);
		}
	}
	
	planning_info(start, i, q_i, q_f);
	
}
//-------------------------------------------------------------------------------------




