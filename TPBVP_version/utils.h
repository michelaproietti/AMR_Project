#include <unistd.h>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <math.h>
#include <cmath>
#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>
#include "problem.h"

using namespace std;

double dt; // control timestep (can be modified directly in the V-REP scene)
int n_plan, k_plan; // number of configurations in the plan and index of current configuration

simInt hRobot, hObstacle, hSofa, hConferenceTable, hIndoorPlant, hIndoorPlant0, hObstacle1;
simInt hIndoorPlant1, hIndoorPlant2, hIndoorPlant3, hIndoorPlant4, hIndoorPlant5, hSofa0, hSofa1, hRack, hRack0, hCupboard, hCupboard0; // object handles
simInt hw1, hw2, hw3, hw4, hw5, hw6, hw7, hw8, hw9, hw10, hw11, hw12, hw13;
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

float plan_cost;

double tpbvp_duration = 0.0;
int goal_index = -1;
int goal_iters = -1;
int goal_updates = 0;

// TUNE FOR SIMULATIONS
int rrt_version = 1;		// 1 for rrt* and 2 for rrt* with motion primitives
int nScene = 3;				// scene on which we perform the simulation
int ground_DIM = 5;  		// 2.5 for scene 1 and 5 for scene 2 and 3
float goal_tolerance = 0.2; // we accept final points withing this distance from qFin
int maxIter = 100;
float grid_resolution = 1.0;
int create_dictionary = 0;
float r = 0.2;  // for exploration-exploitation

int discrete_thetas = 8;
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
		
	while(parent != root) {	// if not root
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

//--------------------------------------------------------------------------_//

// finds the path (that minimizes a certain cost fne) between q_ini e q_fin (zRand e zNearest)
vector<Eigen::Vector3d> TPBVP(Eigen::Vector3d q_i, Eigen::Vector3d q_f){
	
	vector<Eigen::Vector3d> plan; // the plan (configuration-space trajectory)
	
	double v; // nominal linear velocity 
	double omega; // nominal angular velocity
	double dur; // motion duration 
	double t_init = 0.0; // FIXME                
	int N = 5;
	
	// compute evolution of v, omega and motion duration solving a TPBVP   
	// 1. define the problem
	ifopt::Problem nlp;
	nlp.AddVariableSet  (make_shared<ifopt::ExVariables>());
	nlp.AddConstraintSet(make_shared<ifopt::ExConstraint>());
	nlp.AddCostSet      (make_shared<ifopt::ExCost>());   
	// 2. choose solver and options
	ifopt::IpoptSolver ipopt;
	ipopt.SetOption("linear_solver", "mumps");
	ipopt.SetOption("tol", 0.00001);
	ipopt.SetOption("jacobian_approximation", "finite-difference-values"); //"exact");    
	ipopt.SetOption("print_timing_statistics", "no"); 
   	ipopt.SetOption("print_level", 0);
   	ipopt.SetOption("gamma_theta", 1e-8);
	ipopt.SetOption("constr_viol_tol", 1e-12);   
	// 3. set initial and final configurations
	ifopt::setTPBVP(q_i, q_f);   
	// 4. solve and retrieve solution   
	ipopt.Solve(nlp);
	
	Eigen::VectorXd sol = nlp.GetOptVariables()->GetValues();       
	//cout << sol.transpose() << endl;

	Eigen::VectorXd x_list(N+1);
	Eigen::VectorXd y_list(N+1);
	Eigen::VectorXd theta_list(N+1);
	Eigen::VectorXd v_list(N);
	Eigen::VectorXd omega_list(N);
	Eigen::VectorXd dur_list(N);
	
	x_list = sol.segment(0, N+1);
	y_list = sol.segment(N+1, N+1);
	theta_list = sol.segment(2*(N+1), N+1);
	v_list = sol.segment(3*(N+1), N);
	omega_list = sol.segment(3*(N+1)+N, N);
	dur_list = sol.segment(3*(N+1)+2*N, N);	
	
	double x, y, theta, x_in, y_in;
	Eigen::Vector3d q;
	int status = ipopt.GetReturnStatus();

	if(status != 0){
		plan.clear();
		n_plan = plan.size();
		k_plan = 0;
		std::cout << "Planning Failed!" << std::endl;
		x = q_i[0];
		y = q_i[1];
		theta = q_i[2];	
	}
	else {
		Eigen::VectorXd t_list(N+1);
		t_list(0) = t_init;	
		for(int i = 1; i < N+1; i++) t_list(i) = t_list(i-1) + dur_list(i-1);

		plan.clear();
		//Eigen::VectorXd s(4); // an element of 'plan' (q,t)

		Eigen::Vector3d q_k;
		Eigen::Vector3d q_kp1;
		Eigen::Vector3d q_int;
		double t = t_list(0);
		double a = 0;
		
		int k = 0;
		float x_k = x_list(k);
		float y_k = y_list(k);
		float theta_k = theta_list(k);
		q_k << x_k, y_k, theta_k;
		q_kp1 << x_list(k+1), y_list(k+1), theta_list(k+1);
		double t_k = t_list(k);
		double t_kp1 = t_list(k+1);
		
		while(k < N){
			
			theta = theta_k + omega_list(k) * (t - t_k);
			x = x_k + v_list(k) / omega_list(k) * (sin(theta) - sin(theta_k));
			y = y_k - v_list(k) / omega_list(k) * (cos(theta) - cos(theta_k));
			
			q << x, y, theta;
			plan.push_back(q);

			x_k = x;
			y_k = y;
			theta_k = theta;
			
			t_k = t;			
			t = t + dt;	

			if(t > t_kp1){
				k++;
				a = 0;
				t = t_kp1;
				t_k = t_list(k);
				if (k < N) {
					t_kp1 = t_list(k+1);
				}
			}
		}
		plan_cost = t-t_init;
	}
	
	x_in = x;
	y_in = y;
	double theta_p, alpha, m, inter;
	
	if(abs(q_f[0] - x) >= 10e-4){
		m = (q_f[1] - y)/(q_f[0] - x);
		alpha = atan(m);		//angolo di puntamento
		inter = (q_f[0]*y -q_f[1]*x)/(q_f[0] - x);
	}
	else{
		if(q_f[1] > y) alpha = M_PI_2;	//angolo di puntamento
		else alpha = -M_PI_2;		//angolo di puntamento
	}
	theta_p = alpha - theta;	//angolo da fare
	while(abs(theta - alpha) >= 10e-4){	//prima rotazione
		//cout << theta << endl;
		theta += 0.05*theta_p;
		q << x, y, theta;
		plan.push_back(q);		
		plan_cost += dt;
	}
	//traslazione in linea retta
	if(abs(q_f[0] - x) >= 10e-4){
		while(abs(q_f[0] - x) >= 10e-4){
			x = x + 0.05*(q_f[0] - x_in);
			y = m*x + inter;
			q << x, y, theta;
			plan.push_back(q);
			plan_cost += dt;
		}
	}
	else{
		while(abs(q_f[1] - y) >= 10e-4){
			y = y + 0.05*(q_f[1] - y_in);
			q << x, y, theta;
			plan.push_back(q);
			plan_cost += dt;
		}
		
	}
	double theta_in = theta;  // theta è la mia attuale orientazione, raggiunta dopo il secondo movimento della coin
	while(abs(theta - q_f[2]) >= 10e-4){	//rotazione finale
		theta = theta + 0.05*(q_f[2] - theta_in);
		q << x, y, theta;
		plan.push_back(q);
		plan_cost += dt;
	}

	n_plan = plan.size();
	k_plan = 0;


	//cout << "Planning Completed" << std::endl;

	return plan;
}


//------------------------------------------------------------------------------------

// returns a new configuration qNew at the middle of the path (computed with the TPBVP) between zNearest and zRand
vector<Eigen::Vector3d> steer(Eigen::Vector3d zNearest, Eigen::Vector3d zRand) {
	
	time_t start, end;
	time(&start);
	vector<Eigen::Vector3d> plan = TPBVP(zNearest, zRand);
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
	fout << "Goal found at iteration = " << goal_iters << endl;
	fout << "Goal index = " << goal_index << endl;
	
	//goal cost (in time)
	double cost_goal = tree[goal_index].z[5];
	fout << "Goal cost (time duration of the optimal path) = " << cost_goal << endl;
	
	fout << "Goal updates = " << goal_updates << endl;
	
	fout << "Final size of the tree = " << tree_size << endl;
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

	} else if (dist_x <= abs(q_f[0]-tree[goal_index].z[0]) && dist_y <= abs(q_f[1]-tree[goal_index].z[1]) && dist_theta <= abs(q_f[2]-tree[goal_index].z[2]) ) {    //(qNew[5] < tree[goal_index].z[5]) {
		goal_index = qNew[3];
		goal_iters = i;
		goal_updates = goal_updates + 1;
		
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
