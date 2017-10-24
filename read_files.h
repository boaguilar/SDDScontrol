#ifndef READ_FILES_H

#define READ_FILES_H

#include <iostream>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <ctime>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <random>
#include <math.h>

using namespace std;



bool Get_all_files(string file, string& nv_file, string& varf_file,string& tt_file,string& prop_file,string& cnodes_file, string& cedges_file,

				   string& cost_file,int& NumNodes,int& p,long long int& sparse_s,int& sparse_c,int& sparse_h,int& NumSteps, bool& sparse_flag, string& Noise);

void Get_nv_and_maxinput(string nv_file, int& MaxInputs, int NumNodes, int* nv);

void Get_varf(string varf_file, int MaxInputs, int NumNodes, int** varf);

void Get_tt(string tt_file, int MaxInputStates, int NumNodes, int** tt);

void Get_prop(string prop_file, int NumNodes, float** prop);

void read_cedges(string cedges_file, int* ActionHeads, int* ActionTails, int* v_edges, float* CEdgesWeight);

void read_cnodes(string cnodes_file, int* ActionNodes, int* v_nodes, float* CNodesWeight);

void read_cost(string cost_file, int* Badstate, float* Wi, int NumNodes);

int get_lines(string file);



#endif
