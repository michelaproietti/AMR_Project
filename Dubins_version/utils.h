#include <unistd.h>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <math.h>
#include <cmath>
#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>
#include "problem.h"
#include "dubins.h"

using namespace std;

double dt; // control timestep (can be modified directly in the V-REP scene)
int n_plan, k_plan; // number of configurations in the plan and index of current configuration

simInt hRobot, hObstacle, hSofa, hConferenceTable, hIndoorPlant, hIndoorPlant0, hIndoorPlant6, hIndoorPlant7, hObstacle1, hDoor, hDoor1, hWindow;
simInt hIndoorPlant1, hIndoorPlant2, hIndoorPlant3, hIndoorPlant4, hIndoorPlant5, hSofa0, hSofa1, hSofa2, hRack, hRack0, hCupboard, hCupboard0, hCupboard1, hCupboard2; // object handles
simInt hw1, hw2, hw3, hw4, hw5, hw6, hw7, hw8, hw9, hw10, hw11, hw12, hw13, hw14, hw15, hw16, hw17, hw18;
double zIni; // z-component of the robot position (it will be constant during the simulation)

struct VERTEX {
	Eigen::VectorXd z;	// x, y, theta, index of the node, index of the father node, cost(node, root)
	vector<Eigen::Vector3d> parent_to_node;	//path from the parent to the node
};

// Initializing
vector<int> path_config;   // contains the indices of the nodes of the optimal plan. It is filled in optimal_path()
vector<Eigen::Vector3d> path;
vector<VERTEX> tree;
vector<Eigen::Vector3d> parent_plan;
vector<Eigen::VectorXd> path_points; 

float plan_cost;

double tpbvp_duration = 0.0;
int goal_index = -1;
int goal_iters = -1;
int goal_updates = 0;

// TUNE FOR SIMULATIONS
int rrt_version = 2;		// 1 for rrt* and 2 for rrt* with motion primitives
int nScene = 4;				// scene on which we perform the simulation
int ground_DIM = 5;  		// 2.5 for scene 1 and 5 for scene 2 and 3
float goal_tolerance = 0.2; // we accept final points withing this distance from qFin
int maxIter = 10000;
float grid_resolution = 0.5;
int create_dictionary = 0;
float r = 0.2;  // for exploration-exploitation

int discrete_thetas = 8;
//vector<float> thetas = {-M_PI, -M_PI_2, 0, M_PI_2};
vector<float> thetas = {-M_PI, -3*M_PI_4, -M_PI_2, -M_PI_4, 0, M_PI_4, M_PI_2, 3*M_PI_4};


// RANDOM SAMPLES
//Uniform real distribution === Random number distribution that produces 
//floating-point values according to a uniform distribution, 
random_device rd;
mt19937 generator(rd());

//---------------------------------------------------------------------//

double distance(Eigen::VectorXd p1, Eigen::VectorXd p2, float w = 0.2) {
	double x = p1[0]-p2[0];
	double y = p1[1]-p2[1];
	double theta = p1[2]-p2[2];
	return sqrt(x*x+y*y) + w*abs(theta);
}

//---------------------------------------------------------------------//

float cost_plan(vector<Eigen::Vector3d> plan) {
	int s = plan.size();
	float dist_tot = 0.0;
	float d;
	for (int i = 0; i < s-1; i++) {
		d = distance(plan.at(i), plan.at(i+1));
		dist_tot += d;
	}
	return dist_tot;
}

//---------------------------------------------------------------------//

vector<uniform_real_distribution<float_t>> initialize_generators() {
	
	vector<uniform_real_distribution<float_t>> distributions;
	
	//std::default_random_engine generator;
	uniform_real_distribution<float_t> distribution(-ground_DIM, ground_DIM);
	distributions.push_back(distribution);
	
	if (rrt_version == 1) {
		uniform_real_distribution<float_t> orientation_dist(-M_PI, M_PI);
		uniform_real_distribution<float_t> half(0, 1);
		distributions.push_back(orientation_dist);
		distributions.push_back(half);
	}	else {
		uniform_real_distribution<float_t> orientation_dist(0, discrete_thetas-1);
		uniform_real_distribution<float_t> half(0, 1);
		distributions.push_back(orientation_dist);
		distributions.push_back(half);
	}
	
	return distributions;
}

//------------------------------------------------------------------------//


// randomly samples a state from the obstacle-free region of the state space.
Eigen::VectorXd sampleRand(vector<uniform_real_distribution<float_t>> distributions, Eigen::Vector3d q_f) {
	Eigen::VectorXd coords(6);

	// exploration - exploitation
	float p = distributions[2](generator);
	if (p < r) {
		float x = q_f[0];
		float y = q_f[1];
		float theta = q_f[2];
		coords << x, y, theta, -1, -1, 0.0;
	} else {
	
		if (rrt_version == 1) {
			float x = distributions[0](generator);
			float y = distributions[0](generator);
			float theta = distributions[1](generator);
			
			coords << x, y, theta, -1, -1, 0.0;  //-1 xk zRand will not be inserted in the tree
		} else {
			float x = round(distributions[0](generator));
			float y = round(distributions[0](generator));
					
			float half_x = round(distributions[2](generator));
			float half_y = round(distributions[2](generator));
					
			if (half_x) x += grid_resolution;
			if (half_y) y += grid_resolution;
			
			int theta_idx = distributions[1](generator);
			
			coords << x, y, thetas[theta_idx], -1, -1, 0.0;
		}
	}
	
	return coords;
}


//---------------------------------------------------------------------//

// given a state qNew, a tree and a number N, the Znear returns the vertices,
// that are the nodes in the tree, that are near qNew. 
// --> it defines q' s.t. dist(q', qNew)<radius, where radius is computed based on N (CHECK HOW)
vector<Eigen::VectorXd> NearByVertices(Eigen::Vector3d qNew) {
	vector<Eigen::VectorXd> Znear;
	for(int i = 0; i < tree.size(); i++){
		if (Znear.size() < 5)
			Znear.push_back(tree[i].z);
		else {
			double d = distance(tree[i].z, qNew);
			double max_d = 0;
			int index = -1;
			for(int j = 0; j < Znear.size(); j++){
				double dis = distance(Znear[j], qNew);
				if(dis > max_d){
					max_d = dis;
					index = j;
				}
			}
			if(max_d > d) Znear[index] = tree[i].z;
		}
	}
	return Znear;
}

//-----------------------------------------------------------------------//

void optimal_path() {
	
	//vector<int> path_config;
	path_config.clear();
	path.clear();
	
	Eigen::VectorXd root = tree.front().z;
	cout << " ROOT = " << root << endl;
	
	Eigen::VectorXd qfinal = tree[goal_index].z;
	path_config.push_back(qfinal[3]);  // index
	cout << "SUCCESS NODE = " << tree[qfinal[3]].z.head(3) << endl;
		
	Eigen::VectorXd parent = tree[qfinal[4]].z;  //node.z
	cout << "FATHER OF THE SUCCESS NODE = " << parent.head(3) << endl;
		
	while(parent[3] != 0) {	// if not root   //parent != root
		cout << " IN WHILE" << endl;
		path_config.push_back(parent[3]);   //parent[3] is the index of the node in the tree structure
		cout << "NODE = " << tree[parent[3]].z.head(3) << endl;
		parent = tree[parent[4]].z;
	}
	 
	 cout << " OUT OF WHILE" << endl;
	 path_config.push_back(parent[3]);  //aggiunta root
	 
	 cout << "number of nodes in the optimal path = " << path_config.size() << endl;
	 
	 for (int i=path_config.size()-2; i>=0; i--) {
		 cout << " IN THE FOR_1" << endl;
		 vector<Eigen::Vector3d> par_to_node = tree[path_config[i]].parent_to_node;
		 auto lenn = par_to_node.size();
		 cout << "len_path parent to node: " << lenn << endl;
		 for (int j=0; j<par_to_node.size(); j++) {
			 path.push_back(par_to_node[j]);
		 }
	 }
	 
	 auto len_path = path.size();
	 cout << " PATH len ===" <<  len_path  <<endl;
}

//------------------------------------ DUBINS VEHICLE --------------------------------------_//

double angleSignedDiff(double a, double b){
	double d = abs(a - b);
	while(d > 2.0*M_PI) d = d - 2.0*M_PI; 

	double r = 0.0;
	if(d > M_PI) r = 2.0*M_PI - d;
	else r = d;
	double sign = 0.0;
	if( (a-b>=0.0 && a-b<=M_PI) || (a-b<=-M_PI && a-b>=-2.0*M_PI) ) sign = +1.0;
	else sign = -1.0;

	r = sign * r;
	return r; 
}

double ConfigurationDistance(Eigen::Vector3d qA, Eigen::Vector3d qB){
	double alpha = 1.0;
	double posDis = sqrt( (qA(0)-qB(0))*(qA(0)-qB(0)) + (qA(1)-qB(1))*(qA(1)-qB(1)) );
	double angDis = abs(angleSignedDiff(qA(2),qB(2)));
	return posDis + alpha*angDis;
}

int fillPlan(double q[3], double x, void* user_data) {
	Eigen::VectorXd s(4);
    s << q[0], q[1], q[2], x;
	path_points.push_back(s);
    return 0;
}

vector<Eigen::Vector3d> dubins_vehicle(Eigen::Vector3d qA, Eigen::Vector3d qB){

	double vMax = 1.0; // max linear vel
	double omegaMax = 5.0; // max angular vel
	double rho = vMax/omegaMax; // turning radius
	DubinsPath path; // Dubins geometric path
	double dl = 0.01; // discretization step of the geometric path

	// compute the Dubins geometric path
    double q1[] = {qA(0), qA(1), qA(2)};
    double q2[] = {qB(0), qB(1), qB(2)};
    //cout << "c1" << endl;
    dubins_shortest_path(&path, q1, q2, rho);	// computes all the possible paths and decides which one is the one with the lowest cost
	
	// retrieve a list of points along it
	path_points.clear(); 
	//cout << "c2" << endl;
	dubins_path_sample_many(&path, dl, fillPlan, NULL);		// samples many configurations along the entire path
	//cout << "dopo c2" << endl;
	
	vector<Eigen::Vector3d> traj;
	if (path_points.size() == 0) {
		cout << "NO SOLUTION" << endl;
		traj.clear();
		return traj;
	}

	// extract the 4 via points along it
	// 0: configuration at which the robot starts turning (initial one)
	// 1: configuration at which the robot starts driving straight 
	// 2: configuration at which the robot (re)starts turning 
	// 3: configuration at which the path ends (final one) 
	vector<Eigen::VectorXd> via_points; // we use the last component of each element to keep information about the "rotation sign"	
	via_points.push_back(path_points[0]);
	bool turningPrev = true;
	bool turning = true;
	double omegaSign;		
	double omegaSignPrev;
	
	// computes via points at which the robot stops and starts turning
	//cout << "before for" << endl;
	for(int i = 1; i < path_points.size(); i++){
		double dTheta = angleSignedDiff(path_points[i](2), path_points[i-1](2));	// returns the signed difference between the two angles
		// tries to understand if it still needs to rotate and in which sense
		if(dTheta > 0.001){
			turning = true;
			omegaSign = +1.0;
		}
		else if(dTheta < -0.001){
			turning = true;
			omegaSign = -1.0;
		}
		else{
			turning = false;   	
			omegaSign = 0.0;
		}
		
		// found via point
		if(turning != turningPrev){
			via_points[via_points.size()-1](3) = omegaSignPrev; // update previous "rotation sign"
			via_points.push_back(path_points[i]); // add new via point
		}

		turningPrev = turning;
		omegaSignPrev = omegaSign;		
	}
	//cout << "after for" << endl;
	

	
	// we add the last point of the path
	//cout << "via points A  " << via_points[via_points.size()-1](3) << endl;
	via_points[via_points.size()-1](3) = omegaSignPrev;
	cout << "prima 1" << endl;
	via_points.push_back(path_points[path_points.size()-1]);
	cout << "dopo 1" << endl;
	via_points[via_points.size()-1](3) = 0.0; // not relevant
	//cout << "via points  " << via_points.size() << endl;

	// compute a (discretized) trajectory compatible with the geometric path using maximum velocities
	int N = 100;
	double dTau = dt/(double)N; // here we are using a finer discretization (time)step compared to that of the final trajectory (the one used for Coppelia play-back)
	
	/*
	if (via_points.size() < 4) {
		//cout << "NO SOLUTION" << endl;
		traj.clear();
		return traj;
	}
	* */
	
	int i_via_point = 0;	
	Eigen::Vector3d q_k = via_points[i_via_point].head(3);
	traj.push_back(q_k);
	omegaSign = via_points[i_via_point](3);
	int iters = 0;
	while(i_via_point < via_points.size()){
		//cout << "i via point = " << i_via_point << endl;
		//cout << "01" << endl;
		//cout << "via point = " << via_points[i_via_point] << endl;
		Eigen::Vector3d qViaPoint = via_points[i_via_point].head(3);
		//cout << "02" << endl;
		double dist = ConfigurationDistance(q_k, qViaPoint);			// Cartesian distance + angular difference

		//cout << "1" << endl;
		if(dist < 0.05){
			q_k = qViaPoint;
			traj.push_back(q_k);
			i_via_point++;
		}
		//cout << "2" << endl;

		// Euler integration		
		omegaSign = via_points[i_via_point-1](3);			// where to rotate
		double v_k = vMax;
		double omega_k = omegaSign * omegaMax;
		//cout << "3" << endl;
		Eigen::Vector3d q_kp1;
		q_kp1(0) = q_k(0) + dTau * v_k * cos(q_k(2));
		q_kp1(1) = q_k(1) + dTau * v_k * sin(q_k(2));
		q_kp1(2) = q_k(2) + dTau * omega_k;
		q_k = q_kp1;
		traj.push_back(q_k);
		//cout << "4" << endl;

		iters++;
		if(iters > 10000){
			traj.clear();
			return traj;
		}
	}

	// retrieve a trajectory discretized with the Coppelia (time)step
	vector<Eigen::Vector3d> subplan;
	for(int i = 0; i < traj.size(); i++){
		if(i%N == 0) subplan.push_back(traj[i]);
	}
	subplan.push_back(traj[traj.size()-1]);

    return subplan;
}

//------------------------------------------------------------------------------------

// returns a new configuration qNew at the middle of the path (computed with the TPBVP) between zNearest and zRand
vector<Eigen::Vector3d> steer(Eigen::Vector3d zNearest, Eigen::Vector3d zRand) {
	
	time_t start, end;
	time(&start);
	vector<Eigen::Vector3d> plan = dubins_vehicle(zNearest, zRand);
	time(&end);
	double duration_TPBVP = double(end-start);
	tpbvp_duration += duration_TPBVP;

	return plan;
}

//---------------------------------------------------------------------//

//has to check if the path between two points lays in the free configuration space (?).
// DA FINIRE, verificare la collisione solo sull'attuale plan.
int check_collisions(vector<Eigen::Vector3d> plan){
	
	for (int i=0; i<plan.size(); i++){
		Eigen::Vector3d s = plan.at(i);
		
		simFloat pRobot[3];
		simFloat eRobot[3];
		pRobot[0] = s(0);
		pRobot[1] = s(1);
		pRobot[2] = zIni;
		eRobot[0] = 0.0;
		eRobot[1] = 0.0;
		eRobot[2] = s(2);
		
		simSetObjectPosition(hRobot, -1, pRobot);  
		simSetObjectOrientation(hRobot, -1, eRobot);
			
		simInt check;
		if (nScene == 1){  // 5x5 simple scene
			check = simCheckCollision(hRobot, hObstacle);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hObstacle1);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hSofa);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hSofa1);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw2);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw3);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw4);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw5);
			if (check!=0) return 0;
			
		} else if (nScene == 4) {
			check = simCheckCollision(hRobot, hObstacle);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hSofa1);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hSofa2);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hRack);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hRack0);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hCupboard);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hCupboard0);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hCupboard1);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hCupboard2);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hIndoorPlant1);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hIndoorPlant2);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hIndoorPlant3);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hIndoorPlant4);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hIndoorPlant5);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hIndoorPlant6);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hIndoorPlant7);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hDoor);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hDoor1);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hWindow);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw3);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw4);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw5);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw6);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw7);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw8);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw9);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw10);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw11);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw12);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw13);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw14);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw15);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw16);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw17);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw18);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw1);
			if (check!=0) return 0;
			

		

		
		
		
 		}else{
			if(nScene == 2) {
				check = simCheckCollision(hRobot, hObstacle);
				if (check!=0) return 0;
				
				check = simCheckCollision(hRobot, hObstacle1);
			if (check!=0) return 0;
			}
			
			check = simCheckCollision(hRobot, hIndoorPlant);
			if (check!=0) return 0;
			
			simInt check = simCheckCollision(hRobot, hSofa);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw3);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw4);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw5);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw6);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw7);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw8);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw9);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw10);
			if (check!=0) return 0;
			
			check = simCheckCollision(hRobot, hw11);
			if (check!=0) return 0;
			
			
			if (nScene == 3) {
				check = simCheckCollision(hRobot, hConferenceTable);
				if (check!=0) return 0;
			
				check = simCheckCollision(hRobot, hw1);
				if (check!=0) return 0;
				
				check = simCheckCollision(hRobot, hw2);
				if (check!=0) return 0;
				
				check = simCheckCollision(hRobot, hw12);
				if (check!=0) return 0;
				
				check = simCheckCollision(hRobot, hw13);
				if (check!=0) return 0;
				
				check = simCheckCollision(hRobot, hIndoorPlant0);
				if (check!=0) return 0;
				
				check = simCheckCollision(hRobot, hIndoorPlant1);
				if (check!=0) return 0;
				
				check = simCheckCollision(hRobot, hIndoorPlant2);
				if (check!=0) return 0;
				
				check = simCheckCollision(hRobot, hIndoorPlant3);
				if (check!=0) return 0;
				
				check = simCheckCollision(hRobot, hIndoorPlant4);
				if (check!=0) return 0;
				
				check = simCheckCollision(hRobot, hIndoorPlant5);
				if (check!=0) return 0;
				
				check = simCheckCollision(hRobot, hSofa0);
				if (check!=0) return 0;
				
				check = simCheckCollision(hRobot, hSofa1);
				if (check!=0) return 0;
				
				check = simCheckCollision(hRobot, hRack);
				if (check!=0) return 0;
				
				check = simCheckCollision(hRobot, hRack0);
				if (check!=0) return 0;
				
				check = simCheckCollision(hRobot, hCupboard);
				if (check!=0) return 0;
				
				check = simCheckCollision(hRobot, hCupboard0);
				if (check!=0) return 0;
			}
		}		
	}
	return 1;
}

//----------------------------------------------------------------------//

int write_summary(double duration, int len_optimal_path, int iters) {
	fstream fout;
	fout.open("planning_history.txt", ios::out);
	
	if (!fout) {
		cout << "Unable to write summary" << endl;
		return -1;
	}
	
	if (fout) cout << "File created" << endl;
	fout << "Planning duration (seconds) = " << duration << endl;
	fout << "Total time spent computing TPBVPs = " << tpbvp_duration << endl;
		
	auto tree_size = tree.size();
	
	fout << "Total number of iterations = " << iters << endl;
	fout << "Final size of the tree = " << tree_size << endl;
	fout << "Goal found at iteration = " << goal_iters << endl;
	fout << "Goal index = " << goal_index << endl;
	
	//goal cost (in time)
	if (goal_index != -1) {
		double cost_goal = tree[goal_index].z[5];
		fout << "Goal cost (time duration of the optimal path) = " << cost_goal << endl;
		
		fout << "Goal updates = " << goal_updates << endl;

		fout << "Length of the optimal path (# traversed nodes) = " << len_optimal_path << endl;
		
		fout << "_____________________________________________________________________________________" << endl;
		fout << "Nodes in the optimal path" << endl;
		int s;
		for(int j=0; j<path_config.size(); j++){
			s = path_config[j];
			fout << tree[s].z[0] << " , " << tree[s].z[1] << " , " << tree[s].z[2] << endl;
		}
		
		fout << "_____________________________________________________________________________________" << endl;
		fout << "Configurations in the optimal path" << endl;

		Eigen::Vector3d  ss;
		for(int m=0; m<path.size(); m++){
			ss = path.at(m);
			fout << ss(0) << " , " << ss(1) << " , " << ss(2) << endl;
		}
	}
	
	fout << "_____________________________________________________________________________________" << endl;
	fout << "All the configurations in the tree" << endl;
	
	for(int i=0; i<tree.size(); i++){
		fout << tree[i].z[0] << " , " << tree[i].z[1] << " , " << tree[i].z[2] << endl;
	}

	
	
	fout.close();
	return 0;
}

//---------------------------------------------------------------------//

int write_summary_1(Eigen::Vector3d q_i, Eigen::Vector3d q_f){
	fstream fout;
	fout.open("planning_history_PATH.txt", ios::out);
	
	if (!fout) {
		cout << "Unable to write summary" << endl;
		return -1;
	}
	
	if (fout) cout << "File created" << endl;
	
	fout << q_i[0] << " , " << q_i[1] << " , " << q_i[2] << endl;
	fout << q_f[0] << " , " << q_f[1] << " , " << q_f[2] << endl;
	
	int tree_size = tree.size();
	
	for(int i=0; i<tree_size; i++){
		for(int j=0; j<tree[i].parent_to_node.size(); j++){
				fout << "[" << tree[i].parent_to_node[j][0] << "," << tree[i].parent_to_node[j][1] << "," << tree[i].parent_to_node[j][2] << "] " ;
		}
		fout << "--------------------------------------" << endl;
	}
	
	
	fout.close();
	return 0;
}


//----------------------------------------------------------------------//

int write_summary_2(Eigen::Vector3d q_i, Eigen::Vector3d q_f, int len_optimal_path_2, int i, time_t start){
	time_t end_1;	
	time(&end_1);
	double duration_goal = double(end_1-start);
	
	fstream fout;
	fout.open("planning_history_OPT_INI.txt", ios::out);
	
	if (!fout) {
		cout << "Unable to write summary" << endl;
		return -1;
	}
	
	if (fout) cout << "File created" << endl;
	
	fout << q_i[0] << " , " << q_i[1] << " , " << q_i[2]  << endl;
	fout << q_f[0] << " , " << q_f[1] << " , " << q_f[2]  << endl;
	
	fout << "Length of the optimal path (# traversed nodes), when just found = " << len_optimal_path_2 << endl;
	fout << "Goal found (the first time) at iteration = " << i << endl;
	fout << "Planning duration (seconds) = " << duration_goal << endl;
	fout << "Goal index (when first time found) = " << goal_index << endl;
	
	//goal cost (in time)
	double cost_goal = tree[goal_index].z[5];
	fout << "Goal cost (time duration of the optimal path) = " << cost_goal << endl;
	
	int path_size = path.size();
	
	Eigen::Vector3d s;
	for(int i=0; i<path_size; i++){
		s = path.at(i);
		fout << s(0) << " , " << s(1) << " , " << s(2) << endl;
	}
	
	fout.close();
	return 0;
}

//----------------------------------------------------------------------//

void update_goal_info(Eigen::VectorXd qNew, int i, double dist_x, double dist_y, double dist_theta, Eigen::Vector3d q_i, Eigen::Vector3d q_f, time_t start) {
	if (goal_index == -1) {
		goal_index = qNew[3];
		goal_iters = i;
		goal_updates = goal_updates + 1;
		
		cout << "*** GOAL FOUND ***  -->  goal_idx = " << goal_index  << endl;
		
		// I record the optimal plan at the first time I meet the goal
		optimal_path();
		int len_optimal_path_2 = path_config.size();
		
		int w2 = write_summary_2(q_i, q_f, len_optimal_path_2, i, start);
		if (w2 == -1) cout << "Unable to write summary_2" << endl;

	} else if (dist_x < abs(q_f[0]-tree[goal_index].z[0]) && dist_y < abs(q_f[1]-tree[goal_index].z[1]) && dist_theta < abs(q_f[2]-tree[goal_index].z[2]) ) {    //(qNew[5] < tree[goal_index].z[5]) {
		goal_index = qNew[3];
		goal_iters = i;
		goal_updates = goal_updates + 1;
		
		optimal_path();
		
		cout << "*** GOAL UPDATED ***  -->  goal_idx = " << goal_index << endl;
	}
}

void update_goal_info_mp(Eigen::VectorXd qNew, int i, double dist, Eigen::Vector3d q_i, Eigen::Vector3d q_f, time_t start) {
	if (goal_index == -1) {
		goal_index = qNew[3];
		goal_iters = i;
		goal_updates = goal_updates + 1;
		
		cout << "*** GOAL FOUND ***  -->  goal_idx = " << goal_index << ", distance from qFin = " << dist << endl;
		
		// I record the optimal plan at the first time I meet the goal
		optimal_path();
		int len_optimal_path_2 = path_config.size();
		
		int w2 = write_summary_2(q_i, q_f, len_optimal_path_2, i, start);
		if (w2 == -1) cout << "Unable to write summary_2" << endl;

	} else if (qNew[5] < tree[goal_index].z[5]) {
		goal_index = qNew[3];
		goal_iters = i;
		goal_updates = goal_updates + 1;
		
		optimal_path();
		
		cout << "*** GOAL UPDATED ***  -->  goal_idx = " << goal_index << ", distance from qFin = " << dist << endl;
	}
}

void planning_info(time_t start, int i, Eigen::Vector3d q_i, Eigen::Vector3d q_f) {
	time_t end;
	
	if (goal_index == -1) {
		cout << "PLANNING FAILED" << endl;
		
		time(&end);
		double duration = double(end-start);
				
		int w = write_summary(duration, -1, i);
		if (w == -1) cout << "Unable to write summary" << endl;
		
		int w1 = write_summary_1(q_i, q_f);
		if (w1 == -1) cout << "Unable to write summary_1" << endl;
				
	} else {
		
		cout << " __________________________________________________________________________ " << endl;
		cout << " ENTERING OPTIMAL PATH" << endl;
		cout << " GOAL NODE = " << q_f << endl;
		optimal_path();					
		int len_optimal_path = path_config.size();
		cout << " after optimal path " << endl;
		
		auto len_path = path.size();
		cout << " checking the PATH len ... " <<  len_path  <<endl;
		//k_path = 0;
		
		time(&end);
		double duration = double(end - start);
		
		int w = write_summary(duration, len_optimal_path, i);
		if (w == -1) cout << "Unable to write summary" << endl;
		
		int w1 = write_summary_1(q_i, q_f);
		if (w1 == -1) cout << "Unable to write summary_1" << endl;
	}
}
