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
#include <math.h>
#include <random>
#include "rrt_star.h"

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

using namespace std;

LIBRARY vrepLib; // the V-REP library that we will dynamically load and bind


Eigen::Vector3d qIni, qFin; // initial and final configuration

int k_path;

 
//-------------------------------------------------------------------------//
void Initialization(){  
	cout << "Initializing..." << endl;	   
	//tree.clear();
	dt = simGetSimulationTimeStep();  
	
	hRobot = simGetObjectHandle("Pioneer_p3dx");	
	simFloat pRobot[3];
	simFloat eRobot[3]; 
	simGetObjectPosition(hRobot, -1, pRobot);
	simGetObjectOrientation(hRobot, -1, eRobot);
	zIni = pRobot[2];
	qIni << pRobot[0], pRobot[1], eRobot[2];

	         
	if (nScene == 1) {
		hObstacle = simGetObjectHandle("Cuboid");
		hObstacle1 = simGetObjectHandle("Cuboid0");
		 
		hSofa = simGetObjectHandle("sofa");
		hSofa1 = simGetObjectHandle("sofa0");
		 
		hw3 = simGetObjectHandle("240cmHighWall400cm");
		hw2 = simGetObjectHandle("240cmHighWall400cm0");
		hw4 = simGetObjectHandle("240cmHighWall400cm1");
		hw5 = simGetObjectHandle("240cmHighWall400cm2");
	} 
	
	else if (nScene == 4) {
		hSofa1 = simGetObjectHandle("sofa0");
		hSofa2 = simGetObjectHandle("sofa1");
		hRack = simGetObjectHandle("rack");
		hRack0 = simGetObjectHandle("rack0"); 
		hCupboard = simGetObjectHandle("cupboard");
		hCupboard0 = simGetObjectHandle("cupboard#0");
		hCupboard1 = simGetObjectHandle("cupboard#1");
		hCupboard2 = simGetObjectHandle("cupboard#2");
		hIndoorPlant1 = simGetObjectHandle("indoorPlant1");
		hIndoorPlant2 = simGetObjectHandle("indoorPlant2");
		hIndoorPlant3 = simGetObjectHandle("indoorPlant3");
		hIndoorPlant4 = simGetObjectHandle("indoorPlant4");
		hIndoorPlant5 = simGetObjectHandle("indoorPlant5");
		hIndoorPlant6 = simGetObjectHandle("indoorPlant6");
		hIndoorPlant7 = simGetObjectHandle("indoorPlant7");
		hDoor = simGetObjectHandle("slidingDoor#0");
		hDoor1 = simGetObjectHandle("slidingDoor");
		hWindow = simGetObjectHandle("window140cm");
		     
		hw3 = simGetObjectHandle("240cmHighWall400cm");
		hw17 = simGetObjectHandle("240cmHighWall400cm0");
		hw4 = simGetObjectHandle("240cmHighWall400cm1");
		hw5 = simGetObjectHandle("240cmHighWall400cm2");
		hw6 = simGetObjectHandle("240cmHighWall400cm3");
		hw7 = simGetObjectHandle("240cmHighWall400cm4");
		hw8 = simGetObjectHandle("240cmHighWall400cm5");
		hw9 = simGetObjectHandle("240cmHighWall400cm6");
		hw14 = simGetObjectHandle("240cmHighWall400cm7");
		hw15 = simGetObjectHandle("240cmHighWall400cm8");
		hw16 = simGetObjectHandle("240cmHighWall400cm9");
		hw18 = simGetObjectHandle("240cmHighWall400cm10");
		
		
		hw10 = simGetObjectHandle("240cmHighWall_L");
		hw11 = simGetObjectHandle("240cmHighWall_L0");
		hw12 = simGetObjectHandle("240cmHighWall_L1");
		hw13 = simGetObjectHandle("240cmHighWall_L2");
		
		hw1 = simGetObjectHandle("240cmHighPillar50cm");
		
	}
	
	else {
	
		// for scene 2 and 3
		if(nScene == 2) {
			hObstacle = simGetObjectHandle("Cuboid");
			hObstacle1 = simGetObjectHandle("Cuboid0");
		}
		
		hSofa = simGetObjectHandle("sofa");
		hIndoorPlant = simGetObjectHandle("indoorPlant");
		hw3 = simGetObjectHandle("240cmHighWall400cm");
		hw4 = simGetObjectHandle("240cmHighWall400cm1");
		hw5 = simGetObjectHandle("240cmHighWall400cm2");
		hw6 = simGetObjectHandle("240cmHighWall400cm3");
		hw7 = simGetObjectHandle("240cmHighWall400cm4");
		hw8 = simGetObjectHandle("240cmHighWall400cm5");
		hw9 = simGetObjectHandle("240cmHighWall400cm6");
		hw10 = simGetObjectHandle("240cmHighWall_L0");
		hw11 = simGetObjectHandle("240cmHighWall_L1");
	 	
		  
		if (nScene == 3) {
			// for scene 3 (office)
			hConferenceTable = simGetObjectHandle("conferenceTable");
			hw1 = simGetObjectHandle("240cmHighPillar50cm");
			hw2 = simGetObjectHandle("240cmHighWall400cm0");
			hw12 = simGetObjectHandle("240cmHighWall_L");
			hw13 = simGetObjectHandle("240cmHighWall_L2");
			hIndoorPlant0 = simGetObjectHandle("indoorPlant0");
			hIndoorPlant1 = simGetObjectHandle("indoorPlant1");
			hIndoorPlant2= simGetObjectHandle("indoorPlant2");
			hIndoorPlant3 = simGetObjectHandle("indoorPlant3");
			hIndoorPlant4 = simGetObjectHandle("indoorPlant4");
			hIndoorPlant5= simGetObjectHandle("indoorPlant5");
			hSofa0 = simGetObjectHandle("sofa0");
			hSofa1 = simGetObjectHandle("sofa1");
			hRack = simGetObjectHandle("rack");
			hRack0 = simGetObjectHandle("rack0"); 
			hCupboard = simGetObjectHandle("cupboard");
			hCupboard0 = simGetObjectHandle("cupboard#0");
			hDoor = simGetObjectHandle("slidingDoor#0");
			hDoor1 = simGetObjectHandle("slidingDoor");
			hWindow = simGetObjectHandle("window140cm");
		}
	}
	
	simInt hGoal = simGetObjectHandle("GoalRegion");	
	simFloat pGoal[3];
	simFloat eGoal[3];
	simGetObjectPosition(hGoal, -1, pGoal);
	simGetObjectOrientation(hGoal, -1, eGoal);
	qFin << pGoal[0], pGoal[1], eGoal[2];

	cout << "Initialization Completed" << endl;
}


void Planning(){
	cout << "Planning started" << endl;
	
	goal_index = -1;
	tpbvp_duration = 0.0;
	
	
	if (rrt_version == 1) RRTstar(qIni, qFin);
	else {
		if (create_dictionary) create_primitives();
		RRTstar_mp(qIni, qFin);
	}
	
	// we set back the position of the robot handle to the q_ini
	simFloat pRobot[3];
	simFloat eRobot[3];
	pRobot[0] = qIni[0];
	pRobot[1] = qIni[1];
	pRobot[2] = zIni;
	eRobot[0] = 0.0;
	eRobot[1] = 0.0;
	eRobot[2] = qIni[2];
	simSetObjectPosition(hRobot, -1, pRobot );  
	simSetObjectOrientation(hRobot, -1, eRobot );
	 
	k_path = 0;
	
	cout << "Planning Completed" << endl;
}

void Execution(){
	auto len_path = path.size();
	
	if(k_path < len_path){
		Eigen::Vector3d s = path.at(k_path);	

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
		
		k_path++;
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
    string currentDirAndPath(curDirAndFile);
    // 2. Append the V-REP library's name:
    string temp(currentDirAndPath);
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
        cout << "Error, could not find or correctly load the V-REP library. Cannot start 'PluginSkeleton' plugin.\n";
        return (0); // Means error, V-REP will unload this plugin
    }
    if (getVrepProcAddresses(vrepLib) == 0) {
        cout << "Error, could not find all required functions in the V-REP library. Cannot start 'PluginSkeleton' plugin.\n";
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
        cout << "Sorry, your V-REP copy is somewhat old. Cannot start 'PluginSkeleton' plugin.\n";
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

		cout << "Simulation Concluded" << endl;
	
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
