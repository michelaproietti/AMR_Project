// Copyright 2006-2013 Dr. Marc Andreas Freese. All rights reserved.
// marc@coppeliarobotics.com
// www.coppeliarobotics.com
//
// -------------------------------------------------------------------
// This file is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//
// You are free to use/modify/distribute this file for whatever purpose!
// -------------------------------------------------------------------
//
// This file was automatically created for V-REP release V3.0.3 on April 29th 2013

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
#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>
#include "problem.h"

#ifdef _WIN32
    #include <shlwapi.h>
    #pragma comment(lib, "Shlwapi.lib")
#endif /* _WIN32 */

#if defined (__linux) || defined (__APPLE__)
    #include <string.h>
    #include <sys/time.h>
    #define _stricmp(x,y) strcasecmp(x,y)
#endif   

#define PLUGIN_VERSION 1
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

LIBRARY vrepLib; // the V-REP library that we will dynamically load and bind

using namespace std;

simInt hRobot; // robot handle
double zIni; // z-component of the robot position (it will be constant during the simulation)
double dt; // control timestep (can be modified directly in the V-REP scene)

Eigen::Vector3d q_i, q_f; // initial and final configuration
std::vector<Eigen::Vector3d> plan; // the plan (configuration-space trajectory)
int n_plan, k_plan; // number of configurations in the plan and index of current configuration
double p = 0.23;  // raggio
double DELTA = 0.01;   // 0.01
double length_current_manuver;
double MAX_V = 1.2;
double MAX_W = 5.25;

void Initialization(){
	std::cout << "Initializing..." << std::endl;	   

	dt = simGetSimulationTimeStep();   
	
	hRobot = simGetObjectHandle("Pioneer_p3dx");	
	simFloat pRobot[3];
	simFloat eRobot[3];
	simGetObjectPosition(hRobot, -1, pRobot);
	simGetObjectOrientation(hRobot, -1, eRobot);
	zIni = pRobot[2];
	q_i << pRobot[0], pRobot[1], eRobot[2];
	
	simInt hGoal = simGetObjectHandle("GoalRegion");	
	simFloat pGoal[3];
	simFloat eGoal[3];
	simGetObjectPosition(hGoal, -1, pGoal);
	simGetObjectOrientation(hGoal, -1, eGoal);
	q_f << pGoal[0], pGoal[1], eGoal[2];

	std::cout << "Initialization Completed" << std::endl;	
}

//-----------------------------------------------------------------------------//
// DUBINS CAR


//  returns the 4 couples (x, y) of points, at the ends of the 4 tangents

vector<Eigen::VectorXd> tangent_lines(double x1, double y1, double x2, double y2) {
	double d_sq = pow(x2 - x1, 2) + pow(y2 - y1, 2);
	
	vector<Eigen::VectorXd> returnVec;
	
	double d = sqrt(d_sq);   // distanza tra i due centri dei cerchi
	double vx = (x2 - x1) / d; 
	double vy = (y2 - y1) / d;  // (vx, vy) del vettore UNITARIO che va dal centro dek cerchio 1 al centro del cerchio 2 (da la direzione)
	for (int sign1 = +1; sign1 >= -1; sign1 -= 2) {
		std::cout << "1" << std::endl;
		double c = (p - sign1 * p) / d;
		if (c*c > 1.0) continue;
		double h = sqrt(max(0.0, 1.0 - c*c));
		
		for (int sign2 = +1; sign2 >= -1; sign2 -= 2) {
			std::cout << "2" << std::endl;
			double nx = vx * c - sign2 * h * vy;
			double ny = vy * c + sign2 * h * vx;
			
			Eigen::VectorXd elem(4);
			elem << x1+p*nx, y1+p*ny, x2+sign1*p*nx, y2+sign1*p*ny;
			
			returnVec.push_back(elem);
		}
	}
	return returnVec;
}

// computes the length of the arc on which we move
// lhs è il punto iniziale della retta di tangenza. rhs è il punto finale
double arcLength(const Eigen::Vector2d center, const Eigen::Vector2d lhs, const Eigen::Vector2d rhs, bool left) {
	Eigen::Vector2d vec1, vec2;
	vec1 << lhs[0] - center[0], lhs[1] - center[1];
	vec2 << rhs[0] - center[0], rhs[1] - center[1];
	
	double theta = atan2(vec2[1], vec2[0]) - atan2(vec1[1], vec1[0]);
	if (theta < -1e-6 && left) theta += 2.0*M_PI;
	else if (theta > 1e-6 && !left) theta -= 2.0*M_PI;
	
	std::cout << "3" << std::endl;
	return fabs(theta*p);
}

// Euclidean distance function
double Norm(const Eigen::Vector2d lhs, const Eigen::Vector2d rhs) {
	return sqrt( (rhs[0] - lhs[0])*(rhs[0] - lhs[0]) + (rhs[1] - lhs[1])*(rhs[1] - lhs[1]) ); 
}


vector<Eigen::Vector2d> RSRtrajectory(vector<Eigen::VectorXd> tangents, Eigen::Vector2d right_ini, Eigen::Vector2d right_fin) {
	double arcL1, arcL2, arcL3;
	vector<Eigen::Vector2d> next;
	
	Eigen::Vector2d nextControl;  // timesteps, steering angle
	if (tangents.size() > 0) {
		
		double steeringAngle = 1.0*5.25;
		double arcL1 = arcLength(right_ini, q_i.head(2), tangents[0].head(2), false);
		double timesteps = arcL1 / DELTA;
		nextControl << timesteps, steeringAngle;
		next.push_back(nextControl);
		
		Eigen::Vector2d fin_point_tang;
		fin_point_tang << tangents[0][2], tangents[0][3];
		
		steeringAngle = 0.0;
		double arcL2 = Norm(tangents[0].head(2), fin_point_tang);
		timesteps = arcL2 / DELTA;
		nextControl << timesteps, steeringAngle;
		next.push_back(nextControl);
		
		steeringAngle = 1.0*5.25;
		double arcL3 = arcLength(right_fin, fin_point_tang, q_f.head(2), false);
		timesteps = arcL3 / DELTA;
		nextControl << timesteps, steeringAngle;
		next.push_back(nextControl);
		
		length_current_manuver = arcL1 + arcL2 + arcL3;
		
	}
	
	return next;
}



Eigen::VectorXd bestCSCTrajectory(Eigen::Vector2d ini_left, Eigen::Vector2d ini_right, Eigen::Vector2d fin_left, Eigen::Vector2d fin_right) {

	// calcoliamo le tangenti tra i cerchi. RRTangents calcola le 4 tangenti tra il cerchio di destra rispetto alla 
	// posizione iniziale e il cerchio di destra della posizione finale.
	// ini_right è il centro del cerchio di destra della posizione iniziale
	vector<Eigen::VectorXd> RRTangents = tangent_lines(ini_right[0], ini_right[1], fin_right[0], fin_right[1]);
	vector<Eigen::VectorXd> LLTangents = tangent_lines(ini_left[0], ini_left[1], fin_left[0], fin_left[1]);
	vector<Eigen::VectorXd> RLTangents = tangent_lines(ini_right[0], ini_right[1], fin_left[0], fin_left[1]);
	vector<Eigen::VectorXd> LRTangents = tangent_lines(ini_left[0], ini_left[1], fin_right[0], fin_right[1]);
	
	std::cout << "x, y ini (tangente) = " << RRTangents[0][0] << " , " << RRTangents[0][1]  << std::endl;
	std::cout << "x, y fin (tangente) = " << RRTangents[0][2] << " , " << RRTangents[0][3]  << std::endl;
	
	//vector<Eigen::Vector2d> next = RSRtrajectory(RRTangents, ini_right, fin_right);
	std::cout << "6" << std::endl;
	
	// verifica quale delle 4 manovre è la più breve...

	return RRTangents[0];
}

void TPBVP(){
	
	plan.clear();
	
	double theta = q_i[2];
	theta += M_PI_2;
	if (theta > M_PI) theta -= 2.0 * M_PI;
	
	double x_il = q_i[0] + p * cos(theta);
	double y_il = q_i[1] + p * sin(theta);
	
	Eigen::Vector2d ini_left;
	ini_left << x_il, y_il;
	
	theta = q_i[2];
	theta -= M_PI_2;
	if (theta < -M_PI) theta += 2.0 * M_PI;

	double x_ir = q_i[0] + p * cos(theta);
	double y_ir = q_i[1] + p * sin(theta);
	
	Eigen::Vector2d ini_right;
	ini_right << x_ir, y_ir;
	
	theta = q_f[2];
	theta += M_PI_2;
	if (theta > M_PI) theta -= 2.0 * M_PI;
	
	double x_fl = q_f[0] + p * cos(theta);
	double y_fl = q_f[1] + p * sin(theta);
	
	Eigen::Vector2d fin_left;
	fin_left << x_fl, y_fl;
	
	theta = q_f[2];
	theta -= M_PI_2;
	if (theta < -M_PI) theta += 2.0 * M_PI;
	
	double x_fr = q_f[0] + p * cos(theta);
	double y_fr = q_f[1] + p * sin(theta);
	
	Eigen::Vector2d fin_right;
	fin_right << x_fr, y_fr;
	/*
	std::cout << "ini_left: " << ini_left[0] << ", " << ini_left[1] << std::endl;
	std::cout << "fin_left: " << fin_left[0] << ", " << fin_left[1] << std::endl;
	std::cout << "ini_right: " << ini_right[0] << ", " << ini_right[1] << std::endl;
	std::cout << "ini_right: " << fin_right[0] << ", " << fin_right[1] << std::endl;
	*/
	Eigen::VectorXd tangent_points = bestCSCTrajectory(ini_left, ini_right, fin_left, fin_right);
	Eigen::Vector2d ini_tangent_point = tangent_points.head(2);
	Eigen::Vector2d fin_tangent_point;
	fin_tangent_point <<  tangent_points[2], tangent_points[3];

	std::cout << "8" << std::endl;
	
	double x = q_i[0];
	double y = q_i[1];
	theta = q_i[2];
	Eigen::Vector3d q;
	
	
	double theta_k;
	
	// Euler integration of the first rotation
	int j = 0;
	while( abs(x - ini_tangent_point[0])>10e-2 ||  abs(y - ini_tangent_point[1])>10e-2) {
		std::cout << "CIAO 1" << std:: endl;
		 
		//x += MAX_V*dt*cos(theta);
		//y += MAX_V*dt*sin(theta);
		//theta -= MAX_W*dt;
		
		theta_k = theta;
		theta -= MAX_W*dt;
		x += MAX_V/MAX_W*(sin(theta)-sin(theta_k));
		y -= MAX_V/MAX_W*(cos(theta)-cos(theta_k));
		
		q << x, y, theta;
		plan.push_back(q);
		j++;
	}
	std::cout << "plan size = " << plan.size() << std::endl;
	
	
	// aggiustamento angolo di puntamento
	
	double alpha, theta_p, m;
	if(abs(fin_tangent_point[0] - ini_tangent_point[1]) >= 10e-4){
		m = (fin_tangent_point[1] - ini_tangent_point[1])/(fin_tangent_point[0] - ini_tangent_point[0]);
		alpha = atan(m);		//angolo di puntamento
	}
	else{
		if(fin_tangent_point[1] > ini_tangent_point[1]) alpha = M_PI_2;	//angolo di puntamento
		else alpha = -M_PI_2;		//angolo di puntamento
	}
	theta_p = alpha - theta;	//angolo da fare
	while(abs(theta - alpha) >= 10e-4){	//prima rotazione
		//cout << theta << endl;
		theta += 0.05*theta_p;
		q << x, y, theta;
		plan.push_back(q);
	}
	
	
	
	// Euler integration of the second rotation
	j = 0;
	while( abs(x - fin_tangent_point[0])>10e-2 ||  abs(y - fin_tangent_point[1])>10e-2 ) {
		std::cout << "CIAO 2" << std:: endl;
		
		//x += MAX_V*dt*cos(theta+MAX_W*dt/2);
		//y += MAX_V*dt*sin(theta+MAX_W*dt/2);
		
		x += MAX_V*dt*cos(theta);
		y += MAX_V*dt*sin(theta);
		std::cout << x << "  " << y << "  " <<theta << std:: endl;
		std::cout << fin_tangent_point[0] << " " << fin_tangent_point[1] << std:: endl;
		q << x, y, theta;
		plan.push_back(q);
		//if (j > 50) break;
		j++;
	}
	std::cout << "plan size = " << plan.size() << std::endl;
	
	// Euler integration of the third rotation
	j = 0;
	while( abs(x - q_f[0])>10e-2  ||  abs(y - q_f[1])>10e-2 ) {
		std::cout << "CIAO 3" << std:: endl;
		
		
		theta_k = theta;
		theta -= MAX_W*dt;
		x += MAX_V/MAX_W*(sin(theta)-sin(theta_k));
		y -= MAX_V/MAX_W*(cos(theta)-cos(theta_k));
		
		//x += MAX_V*dt*cos(theta);
		//y += MAX_V*dt*sin(theta);
		//theta -= MAX_W*dt;
		q << x, y, theta;
		plan.push_back(q);
		
		std::cout << x << "  " << y << "  " <<theta << std:: endl;
		std::cout << q_f[0] << " " << q_f[1] << std:: endl;
		
		//if (j > 50) break;
		j++;
	}
	std::cout << "plan size = " << plan.size() << std::endl;
	
	
	std::cout << "7" << std::endl;
	
	n_plan = plan.size();
	std::cout << "Plan size = " << n_plan << std::endl;
	k_plan = 0;
	
	
	for (int i=0; i<n_plan; i++) {
		//std::cout << plan.at(i) << endl;
		//std::cout << "i = " << i << std::endl;
	}
	std::cout << "8" << std::endl;
}

void Planning(){
	TPBVP();
}

void Execution(){
	if(k_plan < n_plan){
		Eigen::Vector3d s = plan.at(k_plan);

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
		
		k_plan++;
	}
}

// This is the plugin start routine (called just once, just after the plugin was loaded):

VREP_DLLEXPORT unsigned char v_repStart(void* reservedPointer, int reservedInt) {
    // Dynamically load and bind V-REP functions:
    // ******************************************
    // 1. Figure out this plugin's directory:
    char curDirAndFile[1024];
#ifdef _WIN32
    GetModuleFileName(NULL, curDirAndFile, 1023);
    PathRemoveFileSpec(curDirAndFile);
#elif defined (__linux) || defined (__APPLE__)
    getcwd(curDirAndFile, sizeof (curDirAndFile));
#endif
    std::string currentDirAndPath(curDirAndFile);
    // 2. Append the V-REP library's name:
    std::string temp(currentDirAndPath);
#ifdef _WIN32
    temp += "\\v_rep.dll";
#elif defined (__linux)
    temp += "/libv_rep.so";
#elif defined (__APPLE__)
    temp += "/libv_rep.dylib";
#endif /* __linux || __APPLE__ */
    // 3. Load the V-REP library:
    vrepLib = loadVrepLibrary(temp.c_str());
    if (vrepLib == NULL) {
        std::cout << "Error, could not find or correctly load the V-REP library. Cannot start 'PluginSkeleton' plugin.\n";
        return (0); // Means error, V-REP will unload this plugin
    }
    if (getVrepProcAddresses(vrepLib) == 0) {
        std::cout << "Error, could not find all required functions in the V-REP library. Cannot start 'PluginSkeleton' plugin.\n";
        unloadVrepLibrary(vrepLib);
        return (0); // Means error, V-REP will unload this plugin
    }
    // ******************************************

    // Check the version of V-REP:
    // ******************************************
    int vrepVer;
    simGetIntegerParameter(sim_intparam_program_version, &vrepVer);
    if (vrepVer < 20604) // if V-REP version is smaller than 2.06.04
    {
        std::cout << "Sorry, your V-REP copy is somewhat old. Cannot start 'PluginSkeleton' plugin.\n";
        unloadVrepLibrary(vrepLib);
        return (0); // Means error, V-REP will unload this plugin
    }
    // ******************************************

    simLockInterface(1);

    // Here you could handle various initializations
    // Here you could also register custom Lua functions or custom Lua constants
    // etc.

    return (PLUGIN_VERSION); // initialization went fine, we return the version number of this plugin (can be queried with simGetModuleName)
}

// This is the plugin end routine (called just once, when V-REP is ending, i.e. releasing this plugin):

VREP_DLLEXPORT void v_repEnd() {
    // Here you could handle various clean-up tasks

    unloadVrepLibrary(vrepLib); // release the library
}

// This is the plugin messaging routine (i.e. V-REP calls this function very often, with various messages):

VREP_DLLEXPORT void* v_repMessage(int message, int* auxiliaryData, void* customData, int* replyData) { // This is called quite often. Just watch out for messages/events you want to handle
    // Keep following 6 lines at the beginning and unchanged:
    simLockInterface(1);
    static bool refreshDlgFlag = true;
    int errorModeSaved;
    simGetIntegerParameter(sim_intparam_error_report_mode, &errorModeSaved);
    simSetIntegerParameter(sim_intparam_error_report_mode, sim_api_errormessage_ignore);
    void* retVal = NULL;

    // Here we can intercept many messages from V-REP (actually callbacks). Only the most important messages are listed here.
    // For a complete list of messages that you can intercept/react with, search for "sim_message_eventcallback"-type constants
    // in the V-REP user manual.

    if (message == sim_message_eventcallback_refreshdialogs)
        refreshDlgFlag = true; // V-REP dialogs were refreshed. Maybe a good idea to refresh this plugin's dialog too

    if (message == sim_message_eventcallback_menuitemselected) { // A custom menu bar entry was selected..
        // here you could make a plugin's main dialog visible/invisible
    }

    if (message == sim_message_eventcallback_instancepass) { // This message is sent each time the scene was rendered (well, shortly after) (very often)
        // It is important to always correctly react to events in V-REP. This message is the most convenient way to do so:

        int flags = auxiliaryData[0];
        bool sceneContentChanged = ((flags & (1 + 2 + 4 + 8 + 16 + 32 + 64 + 256)) != 0); // object erased, created, model or scene loaded, und/redo called, instance switched, or object scaled since last sim_message_eventcallback_instancepass message
        bool instanceSwitched = ((flags & 64) != 0);

        if (instanceSwitched) {
            // React to an instance switch here!!
        }

        if (sceneContentChanged) { // we actualize plugin objects for changes in the scene

            //...

            refreshDlgFlag = true; // always a good idea to trigger a refresh of this plugin's dialog here
        }
    }

    if (message == sim_message_eventcallback_mainscriptabouttobecalled) { // The main script is about to be run (only called while a simulation is running (and not paused!))
		
    }

    if (message == sim_message_eventcallback_simulationabouttostart) { // Simulation is about to start

		// code written here is run OFF-LINE (before the simulation actually starts)
		Initialization(); 

		Planning();

    }

    if (message == sim_message_eventcallback_simulationended) { // Simulation just ended

		std::cout << "Simulation Concluded" << std::endl;
	
    }

    if (message == sim_message_eventcallback_moduleopen) { // A script called simOpenModule (by default the main script). Is only called during simulation.
        if ((customData == NULL) || (_stricmp("PluginSkeleton", (char*) customData) == 0)) // is the command also meant for this plugin?
        {
            // we arrive here only at the beginning of a simulation
        }
    }

    if (message == sim_message_eventcallback_modulehandle) { // A script called simHandleModule (by default the main script). Is only called during simulation.
	
		// code written here is run ON-LINE (the following function is invoked every dt seconds, then MUST complete in dt seconds)
		Execution();
		
        if ((customData == NULL) || (_stricmp("PluginSkeleton", (char*) customData) == 0)) // is the command also meant for this plugin?
        {
            // we arrive here only while a simulation is running
        }
    }

    if (message == sim_message_eventcallback_moduleclose) { // A script called simCloseModule (by default the main script). Is only called during simulation.
        if ((customData == NULL) || (_stricmp("PluginSkeleton", (char*) customData) == 0)) // is the command also meant for this plugin?
        {
            // we arrive here only at the end of a simulation
        }
    }

    if (message == sim_message_eventcallback_instanceswitch) { // Here the user switched the scene. React to this message in a similar way as you would react to a full
        // scene content change. In this plugin example, we react to an instance switch by reacting to the
        // sim_message_eventcallback_instancepass message and checking the bit 6 (64) of the auxiliaryData[0]
        // (see here above)

    }

    if (message == sim_message_eventcallback_broadcast) { // Here we have a plugin that is broadcasting data (the broadcaster will also receive this data!)

    }

    if (message == sim_message_eventcallback_scenesave) { // The scene is about to be saved. If required do some processing here (e.g. add custom scene data to be serialized with the scene)

    }

    // You can add many more messages to handle here

    if ((message == sim_message_eventcallback_guipass) && refreshDlgFlag) { // handle refresh of the plugin's dialogs
        // ...
        refreshDlgFlag = false;
    }
   
	// Play button is pressed
    if (simGetSimulationState()==17){  
  		
    }
	// Stop button is pressed
	else if(simGetSimulationState()==sim_simulation_stopped){    
		
	}  


    // Keep following unchanged:
    simSetIntegerParameter(sim_intparam_error_report_mode, errorModeSaved); // restore previous settings
    simLockInterface(0);
    return retVal;
}
