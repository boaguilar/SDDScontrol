#include "read_files.h"

bool Get_all_files(string file, string& nv_file, string& varf_file,string& tt_file,string& prop_file,string& cnodes_file, string& cedges_file, string& cost_file,int& NumNodes,int& p,long long int& sparse_s,int& sparse_c,int& sparse_h, int& NumSteps, bool& sparse_flag, string& Noise, int& Lo,  int& L, int& W)
{
	ifstream input_file;
	input_file.open(file.c_str());
	if(input_file.fail()) {
		cout << "Fail to open the file: "<< file <<endl;
		return false;
	}
	string line;
	int NumLines = 0;
	while (getline(input_file, line)) { NumLines ++ ; }
	if (NumLines != 9 && NumLines != 16 && NumLines != 17 ) {
                cout << NumLines << endl;
		cout << file <<" file should contain 11 or 16 lines!"<<endl;
		return false;
	}
	input_file.clear();
	input_file.seekg(0,ios::beg);
	string symbol;
	string data;
	int i = 0;
	while (getline(input_file, line)) {
		i++;
		istringstream iss(line);
		iss >> symbol >> data;
		if (symbol == "p") { p = stoi(data); }
		else if (symbol == "n") { NumNodes = stoi(data); }
		else if (symbol == "nv") { nv_file = data; }
		else if (symbol == "varf") { varf_file = data; }
		else if (symbol == "tt") { tt_file = data; }
		else if (symbol == "cost") { cost_file = data; }
		else if (symbol == "prop") { prop_file = data; }
		else if (symbol == "cedges") { cedges_file = data; }
		else if (symbol == "cnodes") { cnodes_file = data; }
		else if (symbol == "sparse_s") { sparse_s = stoll(data); }
		else if (symbol == "sparse_c") { sparse_c = stoi(data); }
		else if (symbol == "sparse_h") { sparse_h = stoi(data); }
		else if (symbol == "NumSteps") { NumSteps = stoi(data); }
		else if (symbol == "Noise") { Noise = data; }
                else if (symbol == "Lo") { Lo = stoi(data); }
                else if (symbol == "L") { L = stoi(data); }
                else if (symbol == "W") { W = stoi(data); }
		else { cout <<"Invalid format at line "<<i<<" in the "<<file<<endl; return false; }
	}
	if ( ( p < 2) || ( NumNodes < 1) ) {
		cout <<"P shoud larger than 1 and n should larger than 0!" <<endl;
		return false;
	}
        if ( (L > W) || (W > NumNodes) || (L < 1) ) {
                cout <<"double check L and W."<<endl;
                return false;
        }
        if ( (Lo >= L)  || (Lo < 0) ) {
                cout << " double check Lo " << endl;
                return false;
        }  
	if( ( NumLines == 16 || NumLines == 17 )&& sparse_s >=0 && sparse_c >= 1 && sparse_h >= 1 && NumSteps >= 1 ) {
		sparse_flag = true;
	}
	return true;
}


void Get_nv_and_maxinput(string nv_file, int& MaxInputs, int NumNodes, int* nv)
{
	//int* nv = new int[NumNodes];
	MaxInputs = 0;
	ifstream input_file_nv;
	input_file_nv.open(nv_file);
	if(input_file_nv.fail()) {
		cout << "Fail to open: " << nv_file <<endl;
		for ( int i = 0; i < NumNodes; i ++ ) {
			nv[i] = 0;
		}
		return;
	}
	int i = -1;
	int j = -1;
	string line;
	while (getline(input_file_nv, line)) {
		istringstream iss(line);
		i ++;
		if ( i == 1 ) {
			cout << nv_file << " should contain only 1 lines."<<endl;
			break;
		}
		j = -1;
		while(!iss.eof()){
			j++;
			if ( j == NumNodes) {
				break;
			}
			iss >> nv[j];
			if(nv[j] > MaxInputs) { MaxInputs = nv[j] ; }
		}
		if( j != NumNodes - 1 ) {
			cout <<nv_file<<" should exactly contain "<<NumNodes<<" nodes in the first line."<<endl;
		}
	}
	if ( i != 0 ) {
		cout << nv_file << " is not a 1 X " <<NumNodes << " file!!!"<<endl;
	}
	input_file_nv.close();
}

void Get_prop(string prop_file, int NumNodes, float** prop)
{
	//float** prop = new float* [2];
	//for (int i = 0; i < 2; i++){
	//	prop[i] = new float[NumNodes];
	//}
	ifstream input_file_prop;
	input_file_prop.open(prop_file);
	if(input_file_prop.fail()) {
		cout << "Fail to open: " << prop_file <<endl;
		for ( int i = 0; i < 2; i ++ ) {
			for (int j = 0; j < NumNodes; j ++ ) {
				prop[i][j] = 1;
			}
		}
		return;
	}
	int i = -1;
	int j = -1;
	string line;
	while (getline(input_file_prop,line)) {
		istringstream iss(line);
		i ++;
		if ( i == 2 ) {
			cout <<prop_file << " should contain only 2 lines."<<endl;
			break;
		}
		j = -1;
		while(!iss.eof()){
			j++;
			if ( j == NumNodes) {
				break;
			}
			iss >> prop[i][j];
		}
		if ( j != NumNodes -1 ) {
			cout << prop_file <<" need to exactly contain "<<NumNodes<<" nodes in the " << i+1 << "th line"<<endl;
		}
	}
	if ( i != 1 ) {
		cout << prop_file << " does not a 2 X " <<NumNodes << " file!!!"<<endl;
	}
	input_file_prop.close();
}

void Get_tt(string tt_file, int MaxInputStates, int NumNodes, int** tt)
{
	//int max_colum = pow(p,MaxInputs);
	//cout << "maxcolum: " << max_colum <<endl;
	//int** tt = new int* [max_colum];
	for (int i = 0; i < MaxInputStates; i++){
		tt[i] = new int[NumNodes];
	}
	ifstream input_file_tt;
	input_file_tt.open(tt_file);
	if(input_file_tt.fail()) {
		cout << "Fail to open: " << tt_file <<endl;
		for ( int i = 0; i < MaxInputStates; i ++ ) {
			for (int j = 0; j < NumNodes; j ++ ) {
				tt[i][j] = 0;
			}
		}
		return;
	}
	int i = -1;
	int j = -1;
	string line;
	while (getline(input_file_tt, line)) {
		istringstream iss(line);
		i ++;
		if ( i == MaxInputStates ) {
			cout <<tt_file << " should contain only "<<MaxInputStates<< " lines."<<endl;
			break;
		}
		j = -1;
		while(!iss.eof()){
			j++;
			if ( j == NumNodes) {
				break;
			}
			iss >> tt[i][j];
		//       	cout << "tt["<<i<<"]["<<j<<"]: "<<tt[i][j]<<endl;
		}
		if ( j != NumNodes - 1 ) {
			cout << tt_file <<" need to exactly contain "<<NumNodes<<" nodes in the "<<i+1<<"th line"<<endl;
		}
	}
	if ( i != MaxInputStates - 1 ) {
		cout << tt_file << " does not a "<<MaxInputStates<<" X " <<NumNodes << " file!!!"<<endl;
	}
	input_file_tt.close();
}

void Get_varf(string varf_file, int MaxInputs, int NumNodes, int** varf)
{
	//int** varf = new int* [MaxInputs];
	//for (int i = 0; i < MaxInputs; i++){
	//	varf[i] = new int[NumNodes];
	//}
	ifstream input_file_varf;
	input_file_varf.open(varf_file);
	if(input_file_varf.fail()) {
		cout << "Fail to open: " << varf_file <<endl;
		for ( int i = 0; i < MaxInputs; i ++ ) {
			for (int j = 0; j < NumNodes; j ++ ) {
				varf[i][j] = 0;
			}
		}
		return;
	}
	int i = -1;
	int j = -1;
	string line;
	while (getline(input_file_varf, line)) {
		istringstream iss(line);
		i ++;
		if ( i == MaxInputs ) {
			cout <<varf_file << " should contain only "<<MaxInputs<< " lines."<<endl;
			break;
		}
		j = -1;
		while(!iss.eof()){
			j++;
			if ( j == NumNodes) {
				break;
			}
			iss >> varf[i][j];
//			cout <<"varf["<<i+1<<"]["<<j+1<<"]: "<<varf[i][j]<<endl;
			varf[i][j] --;
		}
		if ( j != NumNodes - 1 ) {
			cout <<varf_file <<" need to exactly contain "<<NumNodes<<" nodes in the "<<i+1<<"th line"<<endl;
		}
	}
	if ( i != MaxInputs - 1 ) {
		cout << varf_file << " does not a "<<MaxInputs<<" X " <<NumNodes << " file!!!"<<endl;
	}
	input_file_varf.close();
}

void read_cedges(string cedges_file, int* ActionHeads, int* ActionTails, int* v_edges, float* CEdgesWeight)
{
	ifstream input_file_cedges;
	input_file_cedges.open(cedges_file);
	if(input_file_cedges.fail()) {
		cout << "Fail to open: " << cedges_file <<endl;
		return;
	}
	string line;

	int i = 0;
	while (getline(input_file_cedges, line)) {
		istringstream iss(line);
		iss >> ActionHeads[i] >> ActionTails[i] >> v_edges[i] >>  CEdgesWeight[i] ;
		ActionHeads[i] --;
		ActionTails[i] --;
		i ++;
	}
	input_file_cedges.close();
}

void read_cnodes(string cnodes_file, int* ActionNodes, int* v_nodes, float* CNodesWeight)
{
	ifstream input_file_cnodes;
	input_file_cnodes.open(cnodes_file);
	if(input_file_cnodes.fail()) {
		cout << "Fail to open: " << cnodes_file <<endl;
		return;
	}
	string line;

	int i = 0;
	while (getline(input_file_cnodes, line)) {
		istringstream iss(line);
		iss >> ActionNodes[i] >> v_nodes[i] >> CNodesWeight[i];
		ActionNodes[i] --;
		i ++;
	}
	input_file_cnodes.close();
}

void read_state_cost(string cost_file, int* Badstate, float* Wi, int NumNodes)
{
	ifstream input_file_cost;
	input_file_cost.open(cost_file);
	if(input_file_cost.fail()) {
		cout << "Fail to open: " << cost_file <<endl;
		return ;
	}
	int index = 0;
	string line;
	int i = -1;
	while (getline(input_file_cost, line)) {
		istringstream iss(line);
		i++;
		if ( i == NumNodes ) {
			cout <<cost_file << " has more than "<<NumNodes <<" number of lines!!"<<endl;
			break;
		}
		iss >> index >> Badstate[i] >> Wi[i];
	}
	if ( i != NumNodes - 1 ) {
		cout <<cost_file <<" need to exactly contain "<<NumNodes<<" number of lines!"<<endl;
	}
	input_file_cost.close();
}

int get_lines(string file) 
{
	int num_lines = 0;
	ifstream input_file;
	input_file.open(file);
	if(input_file.fail()) {
		cout << "Fail to open: " << file <<endl;
		return num_lines;
	}
	string line;
	while (getline(input_file, line)) {
		num_lines ++;
	}
	input_file.close();
	return num_lines;
}
