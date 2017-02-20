#include "read_files.h"
#include "mt64.h"

//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//std::default_random_engine generator(seed);
std::random_device seed;
std::mt19937_64 generator(seed());
std::uniform_real_distribution<double> distribution(0.0,1.0);

void nextstate_ia( int *x, int *y, int NumNodes, int NumStates,  int **sdds_tt, int *sdds_nv, int **sdds_varf,  float ** sdds_prop, int *action,
		   int NumNodeEdges ,  int *ActionNodes, int NumCNodes, int *ActionHeads, int *ActionTails, int NumCEdges, int *v_nodes, int *v_edges, int p);
void sdds_nextstate ( int *x, int *y, int NumNodes, int MaxInputs, int *sdds_nv, int ** sdds_varf, int ** sdds_tt, float ** sdds_prop, int *power2 );
void dec2binary ( int *BinState, int NumNodes, int State) ;
float cost_ija (  int *BinState, int *BadState,  int NumNodes, int *actions, int NumNodeEdges, int NumCNodes, int *ActionNodes, float * Wi ); 
void prob_ia ( int *x, int NumNodes, float *Cost_ia, int NumStates,  int **sdds_tt, int *sdds_nv, int **sdds_varf,  float ** sdds_prop,
               int *action, int NumActions,  int *ActionNodes, int NumCNodes, int *ActionHeads, int *ActionTails, int NumCEdges,int * BadState,
               float * Wi, float* pia, int *v_nodes, int *v_edges, int p, float pp );
void operatorTVi ( int NumNodes, int NumStates,  int **sdds_tt, int *sdds_nv, int **sdds_varf,  float ** sdds_prop,   int NumNodeEdges, int NumActions , 
                   int *ActionNodes, int NumCNodes, int *ActionHeads, int *ActionTails, int NumCEdges, int * BadState ,float *W, float *JV, float alpha,
                   float *newJV , int * U,  int *v_nodes, int *v_edges, int p, float pp );
float RecursiveQ (int NumNodes, int NumStates,  int **sdds_tt, int *sdds_nv, int** sdds_varf,  float** sdds_prop, int NumNodeEdges, int NumActions,
                  int *ActionNodes, int NumCNodes, int *ActionHeads, int *ActionTails, int NumCEdges, int p, int * BadState, float *Wi, float alpha,
                  int c, int h, int *x, int *action, int *v_nodes, int *v_edges, float pp, int NumSteps);
float RecursiveQ_Noise(int NumNodes, int NumStates,  int **sdds_tt, int *sdds_nv, int** sdds_varf,  float** sdds_prop, int NumNodeEdges,
                 int NumActions, int *ActionNodes, int NumCNodes, int *ActionHeads, int *ActionTails, int NumCEdges, int p, int * BadState,
                 float *Wi, float alpha,  int c, int h, int *x, int *action, int *v_nodes, int *v_edges, float pp, int NumSteps);
float Sequential_TwoSteps_RecursiveQ (int NumNodes, int NumStates,  int **sdds_tt, int *sdds_nv, int** sdds_varf,  float** sdds_prop, int NumNodeEdges,
                 int NumActions, int *ActionNodes, int NumCNodes, int *ActionHeads, int *ActionTails, int NumCEdges, int p, int * BadState,
                 float *Wi, float alpha,  int c, int h, int *x, int *action_1, int *action_2, int *v_nodes, int *v_edges, float pp);
float normax ( float * array1,  float * array2 , int NumElements );
long long int multinary2dec ( int *con_number, int size, int base);
void dec2multinary ( int *con_number, int size, int base, long long int number);

int main ( int argc, char* argv[] ) {

   cout <<"***************************************************************"<<endl;	
	// Parameters of Sparse
   long long int sparse_s = -1;
   int sparse_c = -1;
   int sparse_h = -1;
   int NumSteps = 1;
   bool sparse_flag = false;
   string Noise = "No";
   // Value Iteration parameter
   float alpha = 0.9 ;  // discounting cost
   int Nmax =100;   //  maximum number of iterations
   float tol = 0.0001; //  convergence criteria
   int iter =1 ; // iteration
   float pp = 0.05;

   //Input Files 
   int NumNodes = 0; 
   int p = 2;
   string nv_file=""; 
   string varf_file="";
   string tt_file="";
   string prop_file="";
   string cnodes_file="";
   string cedges_file="";
   string cost_file="";

   if (argc != 2) {
		cout <<"Input should only be a file which contains all information of a model!!"<<endl;
		return 0;
	}
	string file = argv[1];
	bool test_get_files = Get_all_files(file,  nv_file, varf_file, tt_file,  prop_file,  cnodes_file,cedges_file,
					cost_file,NumNodes, p, sparse_s, sparse_c,sparse_h,NumSteps,sparse_flag, Noise);
	if (!test_get_files) { cout<<"Double check "<<file<<endl; return 0;}

   // Parameters of the mode
   long long int NumStates = pow(p,NumNodes); 
   int MaxInputs = 0;

   int *sdds_nv = new int[NumNodes];
   Get_nv_and_maxinput(nv_file, MaxInputs, NumNodes, sdds_nv);

   int **sdds_varf = new int* [MaxInputs]; 
   for (int i = 0; i < MaxInputs; i++){
		sdds_varf[i] = new int[NumNodes];
	}
   Get_varf(varf_file, MaxInputs, NumNodes, sdds_varf);

   int MaxInputStates = pow(p,MaxInputs);
   int** sdds_tt = new int* [MaxInputStates];
	for (int i = 0; i < MaxInputStates; i++){
		sdds_tt[i] = new int[NumNodes];
	}
   Get_tt(tt_file, MaxInputStates, NumNodes, sdds_tt);

   float** sdds_prop = new float* [2];
	for (int i = 0; i < 2; i++){
		sdds_prop[i] = new float[NumNodes];
	}
	Get_prop(prop_file, NumNodes, sdds_prop);


  // float* pia = new float[NumStates] ;
  // float* Cost_ia = new float[NumStates] ;


   // Parameters
   int * BadState = new int[NumNodes];
   float * Wi = new float[NumNodes];
   read_cost(cost_file,  BadState, Wi,  NumNodes);

   int NumCNodes = get_lines(cnodes_file); // number of control nodes
   int NumCEdges = get_lines(cedges_file);
   int NumNodeEdges = NumCNodes + NumCEdges; // Number of control  nodes plus Number of control arrows.
   int NumActions  = pow(p,NumNodeEdges);


   int *v_nodes = new int[NumCNodes];
   int *v_edges = new int[NumCEdges];
   int *ActionNodes = new int[NumCNodes];// Ids of control nodes
   int *ActionHeads = new int[NumCEdges];
   int *ActionTails = new int[NumCEdges];
   read_cnodes(cnodes_file, ActionNodes, v_nodes);
   read_cedges(cedges_file, ActionHeads,  ActionTails, v_edges);

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
cout <<"*********************************************************************"<<endl;
   cout <<"Varf_Table: " <<endl;
   for ( int i= 0 ; i < MaxInputs; i ++ ) {
           for (int j = 0; j < NumNodes; j ++ ) {
                   cout <<sdds_varf[i][j]+1 << "  ";
           }
           cout <<endl;
   }
   cout <<"nv: ";
   for ( int i = 0; i < NumNodes; i ++) {
   	   cout << sdds_nv[i] <<"  ";
   }
   cout<<endl;
   cout <<"ActionHeads: ";
   for ( int i = 0; i < NumCEdges; i ++ ) {
	cout << ActionHeads[i]+1<<" ";
   }
   cout <<endl;
   cout <<"ActionTails: ";
   for ( int i = 0; i < NumCEdges; i ++ ) {
	cout << ActionTails[i]+1<<" ";
   }
   cout <<endl;
   cout <<"v_edges:     ";
   for ( int i = 0; i < NumCEdges; i++) {
           cout <<v_edges[i] <<" ";
   }
   cout <<endl;
   cout <<"ActionNodes: ";
   for ( int i = 0; i < NumCNodes; i ++) {
	cout << ActionNodes[i]+1<<" ";
   }
   cout <<endl;
   cout<<"v_nodes:      ";
   for ( int i = 0; i < NumCNodes; i++) {
	   cout <<v_nodes[i] <<" ";
   }
   cout <<endl;
   cout << "DesiredState: ";
   for ( int i= 0 ; i < NumNodes; i ++ ) {
	   cout <<BadState[i] << "  ";
   }
   cout<<endl;
   cout << "Weight:       ";
   for ( int i= 0 ; i < NumNodes; i ++ ) {
	   cout <<Wi[i] << "  ";
   }
   cout <<endl;
   cout <<"NumNodes:          " << NumNodes <<endl;
   cout <<"NumStates:         " << NumStates <<endl;
   cout <<"NumActions:        " <<NumActions <<endl;
   cout <<"MaxIndegree:       " << MaxInputs <<endl;
   cout <<"NumCnodes:         " << NumCNodes <<endl;
   cout <<"NumCedges:         " << NumCEdges <<endl;
   cout <<"State of Interest: "<< sparse_s <<endl;
   cout <<"Num of Iterations: " << sparse_h <<endl;
   cout <<"Sample Size:       " << sparse_c <<endl;
   cout <<"NumSteps:          " << NumSteps <<endl;
   cout <<"Noise:             " << Noise <<endl;

cout <<"******************************************************************"<<endl;
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
   int * x = new int[NumNodes] ; 
   int * y = new int[NumNodes] ;
   int * z = new int[NumNodes] ;
   int * action = new int[NumNodeEdges];
/////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/*
if ( sparse_flag ) {
    ofstream output("100next.txt");
    dec2binary( action , NumNodeEdges, 3 );
    dec2multinary ( x, NumNodes, p, sparse_s);
    cout << "Action: ";
    for ( int i = 0; i < NumNodeEdges; i++ ) {
	cout << action[i];
    }
    cout << endl;
    cout <<"x: ";
    for ( int i = 0; i < NumNodes; i ++ ) {
	cout << x[i];
    }
    cout << endl;
    for ( int i = 0; i < 200; i ++ ) {

	nextstate_ia( x, y, NumNodes,  NumStates,  sdds_tt, sdds_nv, sdds_varf, sdds_prop, action,
                                NumNodeEdges, ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, v_nodes, v_edges, p);
	long long int number = multinary2dec( y, NumNodes, p);
	output << number <<" ";
    }
    output.close();
    return 0;  // program will end.
}
*/
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/*
   clock_t start_sparse = clock();

//int test_7states[7] = {42196, 34004, 54516, 56820, 32180, 15764, 13460};
//ofstream output("sparse_TLGL60_1StepH4C6_278654923981302772_pp0.05.txt");
  string file_state = to_string(sparse_s);
  string filename = "sparse_p53_3StepsH2C3_"+file_state+"_noNoise.txt";
  cout <<"Files: "<< filename <<endl;
  ofstream output(filename);
//ofstream output("sparse_TLGL6_1StepH3C6_48_pp0.05"); 
for ( int i = 0; i < 100; i ++ ) {
	//for ( int j = 0; j < NumStates; j ++ ) {
         // for ( int j = 0; j < 7; j ++ ) {
	//	dec2multinary ( x, NumNodes, p, j);
        //        dec2multinary ( x, NumNodes, p, test_7states[j]);
		dec2multinary ( x, NumNodes, p, sparse_s);
        //        dec2multinary ( x ,NumNodes, p, 278654923981302772);
		int c = sparse_c;
		int h = sparse_h ;
		int *a = new int[NumNodeEdges] ;
		int final_action = 0;

		float minQ = 10000000;
		for ( int action = 0 ; action < NumActions; action ++ ) {
			dec2binary( a , NumNodeEdges, action );
			float Q_a = RecursiveQ ( NumNodes, NumStates, sdds_tt, sdds_nv, sdds_varf, sdds_prop, NumNodeEdges,  NumActions, ActionNodes, 
			                       NumCNodes, ActionHeads, ActionTails, NumCEdges, p, BadState, Wi, alpha, c, h, x, a, v_nodes, v_edges, pp);
			if ( Q_a < minQ ) {
				minQ = Q_a ;
				final_action = action;
			}
		}
		output << final_action + 1 <<" ";
//	}
	output << "\n";
}
output.close();
	clock_t end_sparse = clock();
      double elapsed = (double) (end_sparse - start_sparse) / CLOCKS_PER_SEC;
      cout << "Running time: " << elapsed << " seconds" <<endl;

return 0;  // program will end.
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   if ( sparse_flag && Noise == "Yes" ) {
      clock_t start_sparse = clock();
      //dec2binary( x , NumNodes, sparse_s );
	  dec2multinary ( x, NumNodes, p, sparse_s);
      int c = sparse_c;
      int h = sparse_h ;
      int *a = new int[NumNodeEdges] ;


      int final_action = 0;

      float minQ = 10000000;
      for ( int action = 0 ; action < NumActions; action ++ ) {
          dec2binary( a ,NumNodeEdges, action );
          float Q_a = RecursiveQ_Noise ( NumNodes, NumStates, sdds_tt, sdds_nv, sdds_varf, sdds_prop, NumNodeEdges,  NumActions, ActionNodes, 
			           NumCNodes, ActionHeads, ActionTails, NumCEdges, p, BadState, Wi, alpha, c, h, x, a, v_nodes, v_edges,pp, NumSteps);
       
          //cout << "--- " << Q_a ;
          if ( Q_a < minQ ) {
             minQ = Q_a ;
             final_action = action;
          }
      }
      //cout << endl;
	  for ( int i=0 ; i < NumNodes; i++) {
         	 cout << x[i] << " ";
	  }
     	  cout  << " ---> " ;
          dec2multinary ( a, NumNodeEdges, p, final_action );
          for ( int u = 0 ; u <  NumNodeEdges ; u++ ){
         	 cout << " " << a[u] ;
	  }

      cout <<endl;
      clock_t end_sparse = clock();
      double elapsed = (double) (end_sparse - start_sparse) / CLOCKS_PER_SEC;
      cout << "Running time: " << elapsed << " seconds" <<endl;
      cout  <<"*********************************************************************"<<endl;
      return 0;
   }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   if ( sparse_flag && Noise != "Yes" ) {
      clock_t start_sparse = clock();
      //dec2binary( x , NumNodes, sparse_s );
	  dec2multinary ( x, NumNodes, p, sparse_s);
      int c = sparse_c;
      int h = sparse_h ;
      int *a = new int[NumNodeEdges] ;


      int final_action = 0;

      float minQ = 10000000;
      for ( int action = 0 ; action < NumActions; action ++ ) {
          dec2binary( a ,NumNodeEdges, action );
          float Q_a = RecursiveQ ( NumNodes, NumStates, sdds_tt, sdds_nv, sdds_varf, sdds_prop, NumNodeEdges,  NumActions, ActionNodes, 
			           NumCNodes, ActionHeads, ActionTails, NumCEdges, p, BadState, Wi, alpha, c, h, x, a, v_nodes, v_edges,pp, NumSteps);
       
          //cout << "--- " << Q_a ;
          if ( Q_a < minQ ) {
             minQ = Q_a ;
             final_action = action;
          }
      }
      //cout << endl;
	  for ( int i=0 ; i < NumNodes; i++) {
         	 cout << x[i] << " ";
	  }
     	  cout  << " ---> " ;
          dec2multinary ( a, NumNodeEdges, p, final_action );
          for ( int u = 0 ; u <  NumNodeEdges ; u++ ){
         	 cout << " " << a[u] ;
	  }

      cout <<endl;
      clock_t end_sparse = clock();
      double elapsed = (double) (end_sparse - start_sparse) / CLOCKS_PER_SEC;
      cout << "Running time: " << elapsed << " seconds" <<endl;
      cout  <<"*********************************************************************"<<endl;
      return 0;
   }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
   if ( sparse_flag ) {
      clock_t start_sparse = clock();
      //dec2binary( x , NumNodes, sparse_s );
      dec2multinary ( x, NumNodes, p, sparse_s);
      int c = sparse_c;
      int h = sparse_h ;
      int *a = new int[NumNodeEdges];
      int *b = new int[NumNodeEdges];


      int final_action_1 = 0;
      int final_action_2 = 0;

      float minQ = 10000000;
      int action_i = 1;
      int action_j = 1;
      for ( int i = 0 ; i < NumNodeEdges; i ++ ) {
         dec2binary( a , NumNodeEdges, action_i );
	 action_j = 1;
	 for ( int j = 0; j < NumNodeEdges; j ++ ) {
	 	dec2binary( b , NumNodeEdges, action_j );
         	float Q_a =  Sequential_TwoSteps_RecursiveQ ( NumNodes, NumStates, sdds_tt, sdds_nv, sdds_varf, sdds_prop, NumNodeEdges,  NumActions,
					ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, p, BadState, Wi, alpha, c, h, x, a, b, v_nodes, v_edges, pp);
		cout <<"action_1: "<<action_i<<" action_2: "<<action_j<< " cost: "<<Q_a<<endl;
         	if ( Q_a < minQ ) {
			minQ = Q_a ;
			final_action_1 = action_i;
			final_action_2 = action_j;
		}
		action_j = action_j * 2;
	}
	action_i = action_i * 2;
      }
       
          //cout << "--- " << Q_a ;

          
      //cout << endl;
		for ( int i=0 ; i < NumNodes; i++) {
			cout << x[i] << " ";
		}
		cout  << " ---> " ;
		dec2multinary ( a, NumNodeEdges, p, final_action_1 );
		for ( int u = 0 ; u <  NumNodeEdges ; u++ ){
			cout << " " << a[u] ;
		}
		cout <<" and";
		dec2multinary ( b, NumNodeEdges, p, final_action_2 );
		for ( int u = 0 ; u <  NumNodeEdges ; u++ ){
			cout << " " << b[u] ;
		}

      cout <<endl;
      clock_t end_sparse = clock();
      double elapsed = (double) (end_sparse - start_sparse) / CLOCKS_PER_SEC;
      cout << "Running time: " << elapsed << " seconds" <<endl;
      cout  <<"*********************************************************************"<<endl;
      return 0;
   }
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
   if ( sparse_flag ) {
	   ofstream output_nodes("p53_nodes.csv");
	   ofstream output_edges("p53_edges.csv");
	   output_nodes << "Id;Label;Color";
	   output_nodes << "\n";
	   output_edges << "Source;Target;Type;Weight;Color";
	   output_edges << "\n";



      clock_t start_sparse = clock();
      //dec2binary( x , NumNodes, sparse_s );
	  dec2multinary ( x, NumNodes, p, sparse_s);
      int c = sparse_c;
      int h = sparse_h ;
      int *a = new int[NumNodeEdges] ;


      int final_action = 0;

      float minQ = 10000000;
   for ( int f = 0 ; f < 10; f ++ ) {
      for ( int action = 0 ; action < NumActions; action ++ ) {
          dec2binary( a , NumNodeEdges, action );
          float Q_a = RecursiveQ ( NumNodes, NumStates, sdds_tt, sdds_nv, sdds_varf, sdds_prop, NumNodeEdges,  NumActions, ActionNodes, 
			           NumCNodes, ActionHeads, ActionTails, NumCEdges, p, BadState, Wi, alpha, c, h, x, a, v_nodes, v_edges,pp);
       
          //cout << "--- " << Q_a ;
          if ( Q_a < minQ ) {
             minQ = Q_a ;
             final_action = action;
          }
       }
      //cout << endl;

	    output_nodes << f <<";";
	    for ( int i=0 ; i < NumNodes; i++) {
         	 cout << x[i] << " ";
			 output_nodes << x[i];
	    }
		output_nodes << ";" << "\n";


     	  cout  << " ---> " ;
          dec2multinary ( a, NumNodeEdges, p, final_action );
          for ( int u = 0 ; u <  NumNodeEdges ; u++ ){
         	 cout << " " << a[u] ;
	      }
      cout <<endl;
	  nextstate_ia( x, z, NumNodes,  NumStates,  sdds_tt, sdds_nv, sdds_varf, sdds_prop, a, 
				NumNodeEdges, ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, v_nodes, v_edges, p);

	  output_edges << f <<";" << f+1 <<";" << "Directed;" <<"1;Black"<<"\n";
	  cout << "State("<<f+1<<"): ";
	  for ( int i = 0; i < NumNodes; i++ ) {
		 cout << z[i] <<" ";
		 x[i] = z[i];
	  }
      cout << endl;
   }

     output_nodes << 10 <<";";
	    for ( int i=0 ; i < NumNodes; i++) {
			 output_nodes << z[i];
	    }
		output_nodes << ";";


      clock_t end_sparse = clock();
      double elapsed = (double) (end_sparse - start_sparse) / CLOCKS_PER_SEC;
      cout << "Running time: " << elapsed << " seconds" <<endl;
      cout  <<"*********************************************************************"<<endl;
      return 0;
   }
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
   if ( sparse_flag ) {

      clock_t start_sparse = clock();
      //dec2binary( x , NumNodes, sparse_s );
	  dec2multinary ( x, NumNodes, p, sparse_s);
      int c = sparse_c;
      int h = sparse_h ;
      int *a = new int[NumNodeEdges] ;


      int final_action = 0;

//      float minQ = 10000000;
   for ( int f = 0 ; f < 15; f ++ ) {
	float minQ = 10000000;
      for ( int action = 0 ; action < NumActions; action ++ ) {
          dec2binary( a , NumNodeEdges, action );
          float Q_a = RecursiveQ ( NumNodes, NumStates, sdds_tt, sdds_nv, sdds_varf, sdds_prop, NumNodeEdges,  NumActions, ActionNodes, 
			           NumCNodes, ActionHeads, ActionTails, NumCEdges, p, BadState, Wi, alpha, c, h, x, a, v_nodes, v_edges,pp);
       
          //cout << "--- " << Q_a ;
//	cout << "action: " << action <<"  cost: " << Q_a <<endl;
          if ( Q_a < minQ ) {
             minQ = Q_a ;
             final_action = action;
          }
       }
      //cout << endl;



	    for ( int i=0 ; i < NumNodes; i++) {
         	 cout << x[i] << " ";
	    }

     	  cout  << " ---> " ;
        //  dec2multinary ( a, NumNodeEdges, p, final_action );
            dec2binary ( a, NumNodeEdges, final_action);
          for ( int u = 0 ; u <  NumNodeEdges ; u++ ){
         	 cout << " " << a[u] ;
	      }
          cout <<endl;


	       nextstate_ia( x, z, NumNodes,  NumStates,  sdds_tt, sdds_nv, sdds_varf, sdds_prop, a, 
					NumNodeEdges, ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, v_nodes, v_edges, p);

		   

			cout << "State("<<f<<"): ";
		    for ( int i = 0; i < NumNodes; i++ ) {
				cout << z[i] <<" ";
				x[i] = z[i];
			}
			cout << endl;



			nextstate_ia( x, z, NumNodes,  NumStates,  sdds_tt, sdds_nv, sdds_varf, sdds_prop, a, 
					NumNodeEdges, ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, v_nodes, v_edges, p);

			cout << "State("<<f<<"): ";
		    for ( int i = 0; i < NumNodes; i++ ) {
				cout << z[i] <<" ";
				x[i] = z[i];
			}
			cout << endl;

			nextstate_ia( x, z, NumNodes,  NumStates,  sdds_tt, sdds_nv, sdds_varf, sdds_prop, a, 
					NumNodeEdges, ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, v_nodes, v_edges, p);

			cout << "State("<<f<<"): ";
		    for ( int i = 0; i < NumNodes; i++ ) {
				cout << z[i] <<" ";
				x[i] = z[i];
			}
			cout << endl;

   }

      clock_t end_sparse = clock();
      double elapsed = (double) (end_sparse - start_sparse) / CLOCKS_PER_SEC;
      cout << "Running time: " << elapsed << " seconds" <<endl;
      cout  <<"*********************************************************************"<<endl;
      return 0;
   }
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
   if ( sparse_flag ) {
      clock_t start_sparse = clock();

   int * action  = new int [ NumNodeEdges ] ;
   int * x = new int [ NumNodes ] ;
   float * pia = new float [ NumStates ];
   float * Cost_ia = new float [ NumStates ] ;
dec2binary( action , NumNodeEdges, 15 );
dec2multinary ( x, NumNodes, p, sparse_s);

clock_t end_sparse = clock();

for ( int j = 0; j < 15; j ++ ) {
	prob_ia (  x, NumNodes,  Cost_ia, NumStates,  sdds_tt, sdds_nv, sdds_varf,
               sdds_prop,  action, NumNodeEdges , ActionNodes, NumCNodes, ActionHeads,
               ActionTails, NumCEdges,  BadState, Wi,  pia, v_nodes, v_edges,  p,  pp);
    
      
     float max = 0.0;
     int next = 0;
     for ( int i = 0; i < NumStates; i ++ ) {
	if( pia[i] > max ) {
		max = pia[i];
		next = i;
	}
    }	
    int count = 0;
    for ( int i = 0; i < NumStates; i ++) {
	if ( pia[i] == max ) {
		count ++;
	}
    }
   // cout << "MaxProp: " << max << endl;
   // cout << "# of maxprop: " <<count <<  endl;
    cout << "NextState["<<j <<"]: "<<next << " ";
    dec2multinary ( x, NumNodes, p, next);
   // cout << " Binary: ";
    for ( int i = 0; i < NumNodes; i ++ ) {
	cout << x[i];
    }
    cout << endl;
}
 
      double elapsed = (double) (end_sparse - start_sparse) / CLOCKS_PER_SEC;

      cout << "Running time: " << elapsed << " seconds" <<endl;
      cout  <<"*********************************************************************"<<endl;
      return 0;
   }

*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ofstream output_policy("full_policy.txt");
   
   clock_t start_policy = clock();
   
   float maxdiff = 0.0;
   float * JV = new float[NumStates];
   float * newJV = new float[NumStates];
   int * U = new int[NumStates];

   for ( int i = 0; i < NumStates; i++ ) {
      JV[i] = 4.0;
      newJV[i] = JV[i]; 
   }

   operatorTVi( NumNodes, NumStates, sdds_tt, sdds_nv, sdds_varf, sdds_prop, NumNodeEdges, NumActions,  ActionNodes, NumCNodes,
                ActionHeads, ActionTails, NumCEdges, BadState, Wi,  JV, alpha, newJV ,  U,  v_nodes, v_edges, p, pp); 
   
   while( ( normax(JV, newJV, NumStates) > tol ) && ( iter < Nmax ) ) {

      for ( int i=0 ; i < NumStates; i++ ) 
          JV[i] = newJV[i] ;

      operatorTVi( NumNodes, NumStates, sdds_tt, sdds_nv, sdds_varf, sdds_prop, NumNodeEdges, NumActions,  ActionNodes, NumCNodes,
                   ActionHeads, ActionTails, NumCEdges,BadState, Wi, JV, alpha, newJV ,  U ,  v_nodes, v_edges, p, pp);
      iter++; 
   }

   for ( int i=0 ; i < NumStates; i++ )
          JV[i] = newJV[i] ;

   operatorTVi( NumNodes, NumStates, sdds_tt, sdds_nv, sdds_varf, sdds_prop, NumNodeEdges, NumActions,  ActionNodes, NumCNodes,
                ActionHeads, ActionTails, NumCEdges, BadState, Wi,  JV, alpha, newJV ,  U,  v_nodes, v_edges, p, pp);

   clock_t end_policy = clock();

   for ( int State = 0; State < NumStates  ; State++ ) {
     // cout <<  "  " << U[State] << endl;
      cout << State  << "  --->  " ;
      output_policy <<U[State]+1 << "\n";
      dec2multinary( action, NumNodeEdges, p, U[State] );
      for ( int a = 0 ; a <  NumNodeEdges ; a++ ) {  
          cout << " " << action[a] ;
      }

      cout  << endl;
   }
   double run_time = (double)(end_policy - start_policy) / CLOCKS_PER_SEC;
   cout << "Running time: " << run_time <<" seconds" <<endl;
output_policy.close();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


   //delete [] BinState;
   delete [] JV ;
   delete [] newJV ;
   delete [] U ;
   delete [] sdds_nv ;
   delete [] x;
   delete [] y;
   //system("pause");
   return 0;
}

//////////////////////////////////////////////////////////////
// RecursiveQ  Function of Sparse Sample Algorithm
// Inputs
//  SDDS     : The Stochastic Discrete Dynamicl system
//  Actions  : Action class
//  c        : parameter of SSA algorithm
//  h        : parameter of SSA algorithm
//  state    : initial state of the systems

float RecursiveQ(int NumNodes, int NumStates,  int **sdds_tt, int *sdds_nv, int** sdds_varf,  float** sdds_prop, int NumNodeEdges,
                 int NumActions, int *ActionNodes, int NumCNodes, int *ActionHeads, int *ActionTails, int NumCEdges, int p, int * BadState,
                 float *Wi, float alpha,  int c, int h, int *x, int *action, int *v_nodes, int *v_edges, float pp, int NumSteps)
{

   if (h == 0) {
      return 0.0 ;
   }

   float Q = 0.0 ;
   float R = 0.0 ;

   int *y = new int[NumNodes] ;
   int *z = new int [NumNodes];

   int *a = new int[NumNodeEdges];

//   init_genrand64(UINT64_C(0x12345));
//   std::uniform_int_distribution<long long int> dist_states(0,NumStates-1);
   for (int i = 0 ; i < c ; i++ ) {
//	 double r = distribution(generator);
//	 if ( r < pp ) {
//	       long long int rand_state = dist_states(generator);
//                cout <<"rand_state: "<< rand_state <<endl;
//		dec2multinary ( y, NumNodes, p, rand_state);
//	 }
	 for ( int i = 0; i < NumNodes; i ++ ) {
       		 y[i] = x[i];
  	 }

	 for (int v = 0; v < NumSteps ; v++ ) {

		
		nextstate_ia( y, z, NumNodes,  NumStates,  sdds_tt, sdds_nv, sdds_varf, sdds_prop, action, 
				NumNodeEdges, ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, v_nodes, v_edges, p);
		for ( int i = 0; i < NumNodes; i ++ ) {
			y[i] = z[i];
		}
	 }

/*
	for ( int i = 0; i < NumNodes; i ++ ) {
		z[i] = y[i];
	}
	long long int number = multinary2dec( y, NumNodes, p); 
	cout << "dec: " << number;
*/
/*
//	cout << "Action_1: ";
//	for ( int i = 0; i < NumNodeEdges; i ++ ) {
//		cout << action_1[i] << " ";
//	}
//	cout <<"State(x): ";
//	for ( int i = 0; i < NumNodes; i ++ ) {
//		cout << x[i] <<" ";
//	}
	cout << " State: ";
	for ( int i = 0; i < NumNodes; i ++ ) {
		cout << y[i] << " ";
	}
	cout << "Step: "<< i+1 << endl;
*/


	 float minQ = 10000000;
         for ( int action = 0 ; action < NumActions; action ++ ) {
         	 dec2binary( a , NumNodeEdges, action );
         	 float Q_a =  RecursiveQ ( NumNodes, NumStates, sdds_tt, sdds_nv, sdds_varf, sdds_prop, NumNodeEdges,  NumActions,
                       ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, p, BadState, Wi, alpha, c, h-1, z, a, v_nodes,v_edges,pp,NumSteps);
         	 if ( Q_a < minQ ) {
                	 minQ = Q_a ;
                 }  
         }
      
         Q = Q + minQ ;
         R = R + cost_ija(z, BadState, NumNodes, action, NumNodeEdges, NumCNodes, ActionNodes, Wi);
	}
   delete [] y;
   delete [] z;
   delete [] a;
//   return R;
   return ( R + alpha*Q ) / c ; 
}

float RecursiveQ_Noise(int NumNodes, int NumStates,  int **sdds_tt, int *sdds_nv, int** sdds_varf,  float** sdds_prop, int NumNodeEdges,
                 int NumActions, int *ActionNodes, int NumCNodes, int *ActionHeads, int *ActionTails, int NumCEdges, int p, int * BadState,
                 float *Wi, float alpha,  int c, int h, int *x, int *action, int *v_nodes, int *v_edges, float pp, int NumSteps)
{

   if (h == 0) {
      return 0.0 ;
   }

   float Q = 0.0 ;
   float R = 0.0 ;

   int *y = new int[NumNodes] ;
   int *z = new int [NumNodes];

   int *a = new int[NumNodeEdges];

//   init_genrand64(UINT64_C(0x12345));
//   std::uniform_int_distribution<long long int> dist_states(0,NumStates-1);
   for (int i = 0 ; i < c ; i++ ) {
//	 double r = distribution(generator);
//	 if ( r < pp ) {
//	       long long int rand_state = dist_states(generator);
//                cout <<"rand_state: "<< rand_state <<endl;
//		dec2multinary ( y, NumNodes, p, rand_state);
//	 }
		 for ( int i = 0; i < NumNodes; i ++ ) {
       			 y[i] = x[i];
  		 }	


		for (int v = 0; v < NumSteps ; v++ ) {
			double r = distribution(generator);
			if ( r < pp) {
				long long int rand_state = genrand64_intNumNodes(NumNodes);
				dec2multinary ( z, NumNodes, p, rand_state);
			}
			else {
				nextstate_ia( y, z, NumNodes,  NumStates,  sdds_tt, sdds_nv, sdds_varf, sdds_prop, action, 
						NumNodeEdges, ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, v_nodes, v_edges, p);
			}
			for ( int i = 0; i < NumNodes; i ++ ) {
				y[i] = z[i];
			}
		}
/*
	for ( int i = 0; i < NumNodes; i ++ ) {
		z[i] = y[i];
	}
	long long int number = multinary2dec( y, NumNodes, p); 
	cout << "dec: " << number;
*/
/*
//	cout << "Action_1: ";
//	for ( int i = 0; i < NumNodeEdges; i ++ ) {
//		cout << action_1[i] << " ";
//	}
//	cout <<"State(x): ";
//	for ( int i = 0; i < NumNodes; i ++ ) {
//		cout << x[i] <<" ";
//	}
	cout << " State: ";
	for ( int i = 0; i < NumNodes; i ++ ) {
		cout << y[i] << " ";
	}
	cout << "Step: "<< i+1 << endl;
*/


/*
	r = distribution (generator);
	if ( r < pp ) {
		int rand_state = dist_states(generator);
		dec2multinary ( q, NumNodes, p, rand_state);
	}
        else {
		nextstate_ia( z, q, NumNodes, NumStates, sdds_tt, sdds_nv, sdds_varf, sdds_prop, action,
				NumNodeEdges, ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, v_nodes, v_edges, p);
	}
*/

	 float minQ = 10000000;
         for ( int action = 0 ; action < NumActions; action ++ ) {
         	 dec2binary( a , NumNodeEdges, action );
         	 float Q_a =  RecursiveQ ( NumNodes, NumStates, sdds_tt, sdds_nv, sdds_varf, sdds_prop, NumNodeEdges,  NumActions,
                       ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, p, BadState, Wi, alpha, c, h-1, z, a, v_nodes, v_edges, pp,NumSteps);
         	 if ( Q_a < minQ ) {
                	 minQ = Q_a ;
                 }  
         }
      
         Q = Q + minQ ;
         R = R + cost_ija(z, BadState, NumNodes, action, NumNodeEdges, NumCNodes, ActionNodes, Wi);
	}
   delete [] y;
   delete [] z;
   delete [] a;
//   return R;
   return ( R + alpha*Q ) / c ; 
}

float Sequential_TwoSteps_RecursiveQ(int NumNodes, int NumStates,  int **sdds_tt, int *sdds_nv, int** sdds_varf,  float** sdds_prop, int NumNodeEdges,
                 int NumActions, int *ActionNodes, int NumCNodes, int *ActionHeads, int *ActionTails, int NumCEdges, int p, int * BadState,
                 float *Wi, float alpha,  int c, int h, int *x, int *action_1, int *action_2, int *v_nodes, int *v_edges, float pp)
{

   if (h == 0) {
      return 0.0 ;
   }

   float Q = 0.0 ;
   float R = 0.0 ;

   int *y = new int[NumNodes] ;
   int *z = new int [NumNodes];
   int *q = new int [NumNodes];
   for ( int i = 0; i < NumNodes; i ++ ) {
	q[i] = x[i];
   }
//   int *q = new int [NumNodes];
//   int *m = new int [NumNodes];
   int *a = new int[NumNodeEdges];
   int *b = new int[NumNodeEdges];
   
//   std::uniform_int_distribution<long long int> dist_states(0,NumStates-1);
   for (int i = 0 ; i < c ; i++ ) {
//	 double r = distribution(generator);
//	 if ( r < pp ) {
//	       long long int rand_state = dist_states(generator);
//                cout <<"rand_state: "<< rand_state <<endl;
//		dec2multinary ( y, NumNodes, p, rand_state);
//	 }
//	 else {

		nextstate_ia( q, y, NumNodes,  NumStates,  sdds_tt, sdds_nv, sdds_varf, sdds_prop, action_1, 
				NumNodeEdges, ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, v_nodes, v_edges, p);
//	 }
/*
	for ( int i = 0; i < NumNodes; i ++ ) {
		z[i] = y[i];
	}
	long long int number = multinary2dec( y, NumNodes, p); 
	cout << "dec: " << number;
*/
	cout << "Action_1: ";
	for ( int i = 0; i < NumNodeEdges; i ++ ) {
		cout << action_1[i] << " ";
	}
//	cout <<"State(x): ";
//	for ( int i = 0; i < NumNodes; i ++ ) {
//		cout << x[i] <<" ";
//	}
	cout << " State: ";
	for ( int i = 0; i < NumNodes; i ++ ) {
		cout << y[i] << " ";
	}
	cout << "Step: "<< i+1 << endl;


//         r = distribution (generator);
//         if ( r < pp ) {
//		long long int rand_state = dist_states (generator);
//                dec2multinary (z, NumNodes, p, rand_state);
//       	}
//	else {
		nextstate_ia( y, z, NumNodes, NumStates, sdds_tt, sdds_nv, sdds_varf, sdds_prop, action_2,
				NumNodeEdges, ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, v_nodes, v_edges, p);
//	}

		for ( int i = 0; i < NumNodes; i ++) {
			q[i] = z[i];
		}
		cout <<"Action_2: ";
		for ( int i = 0; i < NumNodeEdges; i ++ ) {
			cout << action_2[i]<<" ";
		}
		cout << " State: ";
		for (int i = 0; i < NumNodes; i ++ ) {
			cout << z[i] << " ";
		}
		cout << "Step: "<<i+1 <<endl;
/*
	r = distribution (generator);
	if ( r < pp ) {
		int rand_state = dist_states(generator);
		dec2multinary ( q, NumNodes, p, rand_state);
	}
        else {
		nextstate_ia( z, q, NumNodes, NumStates, sdds_tt, sdds_nv, sdds_varf, sdds_prop, action,
				NumNodeEdges, ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, v_nodes, v_edges, p);
	}

	if ( r < pp ) {
		int rand_state = dist_states (generator);
		dec2multinary ( m, NumNodes, p, rand_state);
	}
	else {
		nextstate_ia ( q, m, NumNodes, NumStates, sdds_tt, sdds_nv, sdds_varf, sdds_prop, action,
				NumNodeEdges, ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, v_nodes, v_edges, p);
	}
*/
/*
         float minQ = 10000000;
	 int action_i = 1;
	 int action_j = 1;
         for ( int i = 0 ; i < NumNodeEdges; i ++ ) {
         	 dec2binary( a , NumNodeEdges, action_i );
		 action_i = action_i*2;
		 action_j = 1;
		 for ( int j = 0; j < NumNodeEdges; j ++ ) {
		 	dec2binary( b , NumNodeEdges, action_j );
		 	action_j = action_j*2;
         	 	float Q_a =  Sequential_TwoSteps_RecursiveQ ( NumNodes, NumStates, sdds_tt, sdds_nv, sdds_varf, sdds_prop, NumNodeEdges,  NumActions,
					ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, p, BadState, Wi, alpha, c, h-1, z, a, b, v_nodes, v_edges, pp);
         		if ( Q_a < minQ ) {
                		minQ = Q_a ;
			}  
		 }
	 }
         Q = Q + minQ ;
*/}
         R = R + cost_ija(z, BadState, NumNodes, action_2, NumNodeEdges, NumCNodes, ActionNodes, Wi);
//    }

   delete [] y;
   delete [] z;
   delete [] q;
//   delete [] m;
   delete [] a;
   delete [] b;
   return R;
//   return ( R + alpha*Q ) / c ; 
}

//////////////////////////////////////////////////////////////
///  Define Operator T
///   x current state ( array of 1 or 0, with a size of NumNodes)
///   pia       : array of 2^n douable elements, this will hold the probabilitites
///   NumNodes  : number of nodes of SDDS
///   MaxInputs : Maximun indegree
void operatorTVi( int NumNodes, int NumStates,  int **sdds_tt, int *sdds_nv, int **sdds_varf,  float ** sdds_prop,   int NumNodeEdges,
                  int NumActions ,  int *ActionNodes, int NumCNodes, int *ActionHeads, int *ActionTails, int NumCEdges, int * BadState,
                  float *Wi, float *JV, float alpha, float *newJV , int * U, int *v_nodes, int *v_edges, int p, float pp  )
{
  
   int * action  = new int [ NumNodeEdges ] ;
   int * x = new int [ NumNodes ] ;
   float * pia = new float [ NumStates ]; 
   float * Cost_ia = new float [ NumStates ] ; 
   float * temp_J = new float [ NumActions  ] ;



   for (int i = 0 ; i < NumStates; i++  ) {
      for (int u = 0 ; u < NumActions; u++ ) {

         dec2binary( action , NumNodeEdges, u );
         //dec2binary( x, NumNodes, i ) ;
		 dec2multinary ( x, NumNodes, p, i);


         float avg_J = 0; 
         float avg_g = 0; 

         prob_ia ( x, NumNodes, Cost_ia, NumStates,  sdds_tt, sdds_nv, sdds_varf, sdds_prop, action, NumNodeEdges, ActionNodes, 
                   NumCNodes, ActionHeads, ActionTails, NumCEdges, BadState, Wi, pia, v_nodes, v_edges, p, pp ) ;


         for (int j =0 ; j < NumStates; j++ ){
            avg_J = avg_J +  pia[j] * JV[j] ; 
            avg_g = avg_g +  pia[j] * Cost_ia[j]  ;
         }

         temp_J[u] = avg_g +  alpha*avg_J ;

      }

      newJV[i] = 1000000.0 ;
      U[i]  = 0 ;
      for (int u = 0 ; u < NumActions; u++ ) {
          if ( temp_J[u] < newJV[i]  ) {
              newJV[i] = temp_J[u] ;
              U[i] =  u ;
          }
      }


   }
   
   delete [] action ;  
   delete [] pia ;
   delete [] Cost_ia ;
   delete [] x ; 
   delete [] temp_J ;
}





//////////////////////////////////////////////////////////////
///  This function compute the probabilities of going from state x to the other states  
///   x current state ( array of 1 or 0, with a size of NumNodes)
///   pia       : array of 2^n douable elements, this will hold the probabilitites
///   NumNodes  : number of nodes of SDDS
///   MaxInputs : Maximun indegree
///   sdds_nv   : array of number of variables per node
///   sdds_varf : array with the nodes ids for each nodes ( size MaxInputs x NumNodes)  
///   sdds_tt   : array with transition table of the SDDS ( size MaxInputStates x NumNodes )
/// Output
///   y   new state

void prob_ia ( int *x, int NumNodes, float *Cost_ia, int NumStates,  int **sdds_tt, int *sdds_nv, int **sdds_varf,
               float ** sdds_prop,  int *action, int NumNodeEdges ,  int *ActionNodes, int NumCNodes, int *ActionHeads,
               int *ActionTails, int NumCEdges, int * BadState, float *Wi, float* pia, int *v_nodes, int *v_edges, int p, float pp)
{
   
   // internal variables
   int * z  = new int [NumNodes] ;
   int * y  = new int [NumNodes] ;
   int * FreeNodes  = new int [NumNodes] ;
   int * newx = new int [NumNodes] ;
   int * copy_newx = new int [NumNodes] ;

  // float* pia = new float[NumStates];
   // initialize the probabilies and cost
   for (int i = 0 ; i < NumStates ; i ++ ) { 
       pia[i] = 0.0 ;
       Cost_ia[i] = 1000;
   }
   
   // initialize the z 
   for (int i = 0 ; i < NumNodes ; i ++ ) {
     //  z[i] = 0 ;
       FreeNodes[i] = 1; 
       newx[i] = x[i] ; 
       //prob_up[i] = sdds_prop[0][i]; //  this probably is not needed
       //prob_dn[i] = sdds_prop[1][i]; //
   }

   // Set control nodes to be zero if control is 1
   for ( int i = 0 ; i < NumCEdges; i++ ) {
       if ( action[i + NumCNodes ] == 1 ) {
       		FreeNodes[ ActionHeads[i] ]  = -1; // Node is marked for control edge
       }
   }

   // Set control nodes to be zero if control is 1
   for (int i = 0; i  <  NumCNodes ; i ++ ) {
       if ( action[i] == 1) {
           newx[ ActionNodes[i] ] = v_nodes[i]; // Node is disable 
           FreeNodes[ ActionNodes[i] ] = 0;//Node is marked as a control node
		   sdds_prop[0][ActionNodes[i]] = 1;
		   sdds_prop[0][ActionNodes[i]] = 1;
       }
   }
   // set up z with initial values
   for (int i = 0 ; i < NumNodes ; i ++ ) {
	   z[i] = newx[i];
   }


/*

       cout << " action: ";
	  for (int i = 0; i < NumNodeEdges; i ++) {
	   cout << action[i] <<" ";
   }
	  cout<<endl;
      cout << " x: ";
   for (int i = 0; i < NumNodes; i ++) {
	   cout << x[i] <<" ";
   }
   cout <<endl;
   cout << " newx: ";
   for (int i = 0; i < NumNodes; i ++) {
	   cout << newx[i] <<" ";
   }
   cout <<endl;

*/



   // compute z, the deterministi next state 
   for (int i = 0; i < NumNodes; i++ ) {

       if ( FreeNodes[i] == 1) {
           int nv = sdds_nv[i];
           int* new_function = new int[nv];
		   for (int h = 0; h < nv; h ++) {
			  // cout<< "sdds_varf[h][i]: "<< sdds_varf[h][i] - 1<< "  h: " << h << "  i:" << i<<endl;
			   new_function[h] = newx[ sdds_varf[h][i] ];
		   }
		    int new_index = multinary2dec ( new_function, nv, p);
           z[i] = sdds_tt[new_index][i];
		   delete [] new_function;
       }
       else if ( FreeNodes[i] == -1 ) {
		   for (int j = 0; j < NumNodes ; j++) {
			   copy_newx[j] = newx[j];
		   }
		   int nv = sdds_nv[i];
		   int* new_function = new int[nv];
		   for (int g = 0; g < NumCEdges; g ++) {
			   if (ActionHeads[g] == i) {
				   copy_newx[ActionTails[g]] = v_edges[g];
			   }
		   }
		   for (int h = 0; h < nv; h++) {
			    new_function[h] = copy_newx[ sdds_varf[h][i] ];
		   }
		   int new_index = multinary2dec ( new_function, nv, p);
           z[i] = sdds_tt[new_index][i];
		   delete [] new_function;
	   }
   }

/* 
   cout <<" Action: ";
   for (int i = 0; i < NumNodeEdges; i++){
	   cout << action[i] <<" ";
   }
   cout<<endl;
   cout <<" X: ";
   for (int i = 0; i < NumNodes; i++){
	   cout << x[i] <<" ";
   }
   cout<<endl;
   cout << " z: ";
  for (int i = 0; i < NumNodes; i ++) {
	   cout << z[i] <<" ";
   }
  cout <<endl;
   cout<<"****************************************"<<endl;
*/

   for ( int state  = 0 ; state < NumStates; state++ ) {
       //dec2binary( y ,  NumNodes, state);
	   dec2multinary ( y, NumNodes, p, state);
	   pia[state] = 1.0;
           Cost_ia[state]= cost_ija( y , BadState,  NumNodes, action, NumNodeEdges , NumCNodes, ActionNodes, Wi );


       for ( int k=0; k<NumNodes; k++ ){
           float c = 0.0;
           
           if (z[k] > newx[k]) {
              if ( y[k] == z[k] ) 
                 c = c + sdds_prop[0][k]; 
 
              if ( y[k] == newx[k] ) 
                 c = c + ( 1.0 - sdds_prop[0][k]) ;             
           }
           else if ( z[k] < newx[k] ) {
              if ( y[k] == z[k] ) 
                 c = c + sdds_prop[1][k];

              if ( y[k] == newx[k] )
                 c = c + ( 1.0 - sdds_prop[1][k]);              
           }
           else {
              if ( y[k] == newx[k]  ) 
                 c = 1.0;
           }
           pia[state] = pia[state] * c ;
       }
   }

   //delete [] Cost_ia ;   
   delete [] FreeNodes ;
   delete [] z ;   
   delete [] y ;   
   delete [] newx ;
   delete [] copy_newx ;

   float p_to_n = 1/(float)NumStates;
   for (int j = 0; j < NumStates ; j ++ ) {
	pia[j] = ( 1 - pp ) * pia[j] + pp*p_to_n;
   }
}

void nextstate_ia (int *x, int *y, int NumNodes, int NumStates,  int **sdds_tt, int *sdds_nv, int **sdds_varf,  float ** sdds_prop, int *action, 
	           int NumNodeEdges, int *ActionNodes, int NumCNodes, int *ActionHeads, int *ActionTails, int NumCEdges, int *v_nodes, int *v_edges, int p)
{

   int * z  = new int [NumNodes] ;
   int * FreeNodes  = new int [NumNodes] ;
   int * newx = new int [NumNodes] ;
   int * copy_newx = new int [NumNodes] ;

   // initialize the z 
   for (int i = 0 ; i < NumNodes ; i ++ ) {
     //  z[i] = 0 ;
       FreeNodes[i] = 1; 
       newx[i] = x[i] ; 
   }

   // Set control nodes to be zero if control is 1
   for ( int i = 0 ; i < NumCEdges; i++ ) {
       if ( action[i + NumCNodes ] == 1 ) {
       		FreeNodes[ ActionHeads[i] ]  = -1; // Node is marked for control edge
       }
   }

   // Set control nodes to be zero if control is 1
   for (int i = 0; i  <  NumCNodes ; i ++ ) {
       if ( action[i] == 1) {
           newx[ ActionNodes[i] ] = v_nodes[i]; // Node is disable 
           FreeNodes[ ActionNodes[i] ] = 0;//Node is marked as a control node
		   sdds_prop[0][ActionNodes[i]] = 1;
		   sdds_prop[0][ActionNodes[i]] = 1;
       }
   }
   // set up z with initial values
   for (int i = 0 ; i < NumNodes ; i ++ ) {
	   z[i] = newx[i];
   }

/*
            cout << " action: ";
	  for (int i = 0; i < NumNodeEdges; i ++) {
	   cout << action[i] <<" ";
   }
	  cout<<endl;
      cout << " x: ";
   for (int i = 0; i < NumNodes; i ++) {
	   cout << x[i] <<" ";
   }
   cout <<endl;
   cout <<" FreeNodes: ";
   for ( int i = 0; i < NumNodes; i ++ ) {
	cout << FreeNodes[i] << " ";
   }
   cout << endl;
   cout << " newx: ";
   for (int i = 0; i < NumNodes; i ++) {
	   cout << newx[i] <<" ";
   }
   cout <<endl;
   cout<<"****************************************"<<endl;
*/
 
   // compute z, the deterministi next state 
   for (int i = 0; i < NumNodes; i++ ) {

       if ( FreeNodes[i] == 1) {
           int nv = sdds_nv[i];
           int* new_function = new int[nv];
		   for (int h = 0; h < nv; h ++) {
			  // cout<< "sdds_varf[h][i]: "<< sdds_varf[h][i] - 1<< "  h: " << h << "  i:" << i<<endl;
			   new_function[h] = newx[ sdds_varf[h][i] ];
		   }
      	   long long int new_index = multinary2dec ( new_function, nv, p);
           z[i] = sdds_tt[new_index][i];
           delete [] new_function;
       }
       else if ( FreeNodes[i] == -1 ) {
		   for (int j = 0; j < NumNodes ; j++) {
			   copy_newx[j] = newx[j];
		   }
		   int nv = sdds_nv[i];
		   int* new_function = new int[nv];
		   for (int g = 0; g < NumCEdges; g ++) {
			   if (ActionHeads[g] == i) {
				   copy_newx[ActionTails[g]] = v_edges[g];
			   }
		   }
		   for (int h = 0; h < nv; h++) {
			    new_function[h] = copy_newx[ sdds_varf[h][i] ];
		   }
		   long long int new_index = multinary2dec ( new_function, nv, p);
                   z[i] = sdds_tt[new_index][i];
		   delete [] new_function;
	   }
   }


/*

   cout <<" Action: ";
   for (int i = 0; i < NumNodeEdges; i++){
	   cout << action[i] <<" ";
   }
   cout<<endl;
   cout <<" X: ";
   for (int i = 0; i < NumNodes; i++){
	   cout << x[i] <<" ";
   }
   cout<<endl;
   cout << " z: ";
   for (int i = 0; i < NumNodes; i ++) {
	   cout << z[i] <<" ";
   }
   cout <<endl;
   cout<<"****************************************"<<endl;

*/

   for ( int k=0; k<NumNodes; k++ ) {
       double r  = distribution(generator); // random generator

       if ( newx[k] > z[k] ) {
           if ( r < sdds_prop[1][k] )
              y[k] = z[k];
           else
              y[k] = newx[k];
       }
       else if ( newx[k] < z[k]  ) {
           if ( r < sdds_prop[0][k] )
              y[k] = z[k];
           else
              y[k] = newx[k];          
       }
       else {
           y[k] = newx[k] ; 
       }
   }

   delete [] FreeNodes ;
   delete [] z ;
   delete [] newx ;
   delete [] copy_newx ;
}




//////////////////////////////////////////////////////////////
/// This function compute the new state after  one step of a SDDS
///  The current state is x  (binary array of n elements, 
///  where  n is # of genes).
/// Inputs 
///   x   current state ( array of 1 or 0, with a size of NumNodes)
///   NumNodes  : number of nodes of SDDS
///   MaxInputs : Maximun indegree
///   sdds_nv   : array of number of variables per node
///   sdds_varf : array with the nodes ids for each nodes ( size MaxInputs x NumNodes)  
///   sdds_tt   : array with transition table of the SDDS ( size MaxInputStates x NumNodes )
/// Output
///   y   new state 
void sdds_nextstate(int *x, int *y, int NumNodes, int MaxInputs, int *sdds_nv, int ** sdds_varf, int ** sdds_tt, float ** sdds_prop, int *power2 ) {

   int *z = new int [NumNodes] ; 
   for (int i=0 ; i < NumNodes ; i++ ) {

       int nv = sdds_nv[i];
       int k = 0;
       for (int j = 0 ; j < MaxInputs ; j ++ ) {
           int Indx =  MaxInputs - nv  + j;
           k = k + power2[Indx]*x[ sdds_varf[j][i] ] ;
       }
       z[i] = sdds_tt[k][i];
   }
   

   for (int i=0 ; i < NumNodes ; i++ ) {
       double r  = distribution(generator); // random generator
       if ( x[i] > z[i] ) {
           if ( r < sdds_prop[1][i] )
              y[i] = z[i];
           else
              y[i] = x[i];
       }
       else if (  x[i] < z[i] ) {
           if ( r < sdds_prop[0][i] )
              y[i] = z[i];
           else
              y[i] = x[i];
       }
       else
           y[i] = x[i] ;
   }
}

//////////////////////////////////////////////////////
/// Comppute the cost( or reward) when the systems
/// transition from a current state "i" to an state "j"
/// ehwn the action "a" was taken.
/// Inputs 
/// BinState   : Represent the sate of "j", it is an array that store the binary number
///              most significant digit is first (index=0)
/// NumNodes   : Number of nodes (size of BinState)
/// actions    : Vector that stores the states of nodes and edges of
///              of an specific action
/// NumActions : Number of control nodes + N of control edlges
/// NumCNodes  : Number of control nodes
float cost_ija(int *BinState, int *BadState,  int NumNodes, int *actions, int NumNodeEdges, int NumCNodes, int *ActionNodes, float * Wi ) { 
   ///// Parameters to compute Cost
   float  Cost = 0.0 ;
   float  Wn = 1.0 ;  // Weigth for control nodes
   float  We = 1.0 ; // Weigth for control edges

   //////////////////////////////////////////////// 

   // cost of control nodes
   for ( int i = 0 ; i < NumCNodes; i++ ) 
      Cost = Cost +  Wn * actions[i] ;

   // cost of control edges
   for ( int i = NumCNodes ; i < NumNodeEdges  ; i++ ) 
      Cost = Cost +  We  * actions[i] ;

   // cost of bad state
   for ( int  i = 0; i < NumNodes; i++){ 
      if ( BinState[i] != BadState[i]) {
         Cost = Cost + Wi[i] ;
      }
   }

   //delete [] Wi ; 
   return Cost ;
   
}



//////////////////////////////////////////////////////
///Convert decimal number to binary   
///Inputs  
///BinState : Vector that store the binary number
///           most significant digit is first (index=0)
///NumNodes : Number of bits(Nodes) of the final binary number
///State    : Decimal number to be converted to binary

///output 
///Binstate
/////////////////////////////////////////////////////
void dec2binary( int *BinState, int NumNodes, int State) {

   int bit = NumNodes -1;
   int number = State;
   do{
      if ( (number & 1) == 0 )
         BinState[bit] = 0 ;
      else
         BinState[bit] = 1 ;

      number >>= 1;
      bit--;
   } while ( number );

   for (int i = bit ; i >= 0 ; i--)
       BinState[i] = 0 ;
} 

/////////////////////////////////////////////////////
/// Function tht compute the maximum diference between
/// two arrays
float normax ( float * array1,  float * array2 , int NumElements ) {

   float maximum = -1.0;
   for (int i = 0 ; i <  NumElements; i++) {
       float diff = array1[i] - array2[i] ;

       if ( diff < 0.0 ) 
           diff = -1.0 * diff ;
              
       if ( diff > maximum )
           maximum = diff ;
   }
   return maximum ;
}

long long int multinary2dec( int *con_number, int size, int base) {

	long long int deci_number = 0;
	long long int tempt = 0;
	deci_number = deci_number + con_number[size - 1];
	if (size == 1) {
		return deci_number;
	}
	int j = size - 1;
	for ( int i = 0; i < size - 1; i ++ ){
		tempt = pow(base,j)*con_number[i];
		//cout << con_number[i] << " " <<base << " " << j << " "<< tempt <<endl;
		deci_number = deci_number + tempt;
		j--;
	}
	return deci_number;
}

void dec2multinary ( int *con_number, int size, int base, long long int number)
{
	long long int divid_number = number;
	for ( int i = size - 1; i >= 0; i -- ) {
		con_number[i] = divid_number % base;
		divid_number = divid_number / base;
	}
}
