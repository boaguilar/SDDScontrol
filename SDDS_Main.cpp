#include "read_files.h"
#include "mt64.h"

//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//std::default_random_engine generator(seed);
std::random_device seed;
std::mt19937_64 generator(seed());
std::uniform_real_distribution<double> distribution(0.0,1.0);

void random_nextstate(int *z, int NumNodes);
void nextstate_ia( int *x, int *y, int NumNodes, int NumStates,  int **sdds_tt, int *sdds_nv, int **sdds_varf,  float ** sdds_prop, int *action,
		   int NumNodeEdges ,  int *ActionNodes, int NumCNodes, int *ActionHeads, int *ActionTails, int NumCEdges, int *v_nodes, int *v_edges, int p);
void sdds_nextstate ( int *x, int *y, int NumNodes, int MaxInputs, int *sdds_nv, int ** sdds_varf, int ** sdds_tt, float ** sdds_prop, int *power2 );
void dec2binary ( int *BinState, int NumNodes, int State) ;
float cost_ija_LW(int *BinState, int *BadState,  int NumNodes, int *action1, int *action2, int NumNodeEdges, int NumCNodes, float * Wi, float *CNodesWeight, float *CEdgesWeight );
float cost_ija (  int *BinState, int *BadState,  int NumNodes, int *actions, int NumNodeEdges, int NumCNodes, int *ActionNodes, float * Wi, float *CNodesWeight, float *CEdgesWeight  ); 
void prob_ia ( int *x, int NumNodes, float *Cost_ia, int NumStates,  int **sdds_tt, int *sdds_nv, int **sdds_varf,  float ** sdds_prop,
               int *action, int NumActions,  int *ActionNodes, int NumCNodes, int *ActionHeads, int *ActionTails, int NumCEdges,int * BadState,
               float * Wi, float* pia, int *v_nodes, int *v_edges, int p, float pp, float *CNodesWeight, float *CEdgesWeight );
void operatorTVi ( int NumNodes, int NumStates,  int **sdds_tt, int *sdds_nv, int **sdds_varf,  float ** sdds_prop,   int NumNodeEdges, int NumActions , 
                   int *ActionNodes, int NumCNodes, int *ActionHeads, int *ActionTails, int NumCEdges, int * BadState ,float *W, float *JV, float alpha,
                   float *newJV , int * U,  int *v_nodes, int *v_edges, int p, float pp, float *CNodesWeight, float *CEdgesWeight );

float RecursiveQ_Noise_LW(int NumNodes, int NumStates, int  MaxInputs,  int **sdds_tt, int *sdds_nv, int** sdds_varf,  float** sdds_prop, int* power2, int NumNodeEdges, int *ActionNodes, int NumCNodes, int *ActionHeads, int *ActionTails, int NumCEdges, int p, int * BadState, float *Wi, float alpha_over_m,  int c, int h, int *x, int ActionIdx, int NumActions,  int *v_nodes, int *v_edges, float pp, int Lo, int L, int W, float *CNodesWeight, float *CEdgesWeight);



float normax ( float * array1,  float * array2 , int NumElements );
long long int multinary2dec ( int *con_number, int size, int base);
void dec2multinary ( int *con_number, int size, int base, long long int number);
void ActionHash ( int *action1, int *action2, int NumNodeEdges, int ActionIndex) ;

int main ( int argc, char* argv[] ) {

   cout <<"***************************************************************"<<endl;	
   // Parameters of Sparse
   long long int sparse_s = -1;
   int sparse_c = -1;
   int sparse_h = -1;
   int Lo = 0 ;
   int L = 2 ;  
   int W = 2 ;
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
   bool test_get_files = Get_all_files(file,  nv_file, varf_file, tt_file,  prop_file,  cnodes_file,cedges_file, cost_file,NumNodes, p, sparse_s, sparse_c,sparse_h,NumSteps,sparse_flag, Noise,Lo, L, W);
  
   if (!test_get_files) { cout<<"Double check "<<file<<endl; return 0;}

   // Parameters of the sdd model
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


   // Parameters accions and cost
   int * BadState = new int[NumNodes];
   float * Wi = new float[NumNodes];
   read_state_cost(cost_file,  BadState, Wi,  NumNodes);

   int NumCNodes = get_lines(cnodes_file); // number of control nodes
   int NumCEdges = get_lines(cedges_file);
   int NumNodeEdges = NumCNodes + NumCEdges; // Number of control  nodes plus Number of control arrows.
   int NumActions = 0; 
   if ( sparse_flag && (Lo > 0) ) {
        NumActions  =  ( NumNodeEdges + 1 )*( NumNodeEdges + 1 ) ;
   }
   else {
        NumActions  = pow(p,NumNodeEdges) ;
   }

   int *v_nodes = new int[NumCNodes];
   int *v_edges = new int[NumCEdges];
   int *ActionNodes = new int[NumCNodes];// Ids of control nodes
   int *ActionHeads = new int[NumCEdges];
   int *ActionTails = new int[NumCEdges];
   float *CNodesWeight = new float[NumCNodes];
   float *CEdgesWeight = new float[NumCEdges];
   
   read_cnodes(cnodes_file, ActionNodes, v_nodes, CNodesWeight);
   read_cedges(cedges_file, ActionHeads,  ActionTails, v_edges, CEdgesWeight);

   ///////////////////////////////////////////////////////////////////
   //// This part creates an array like [2^(n-1) 2^(n-2) ... 2^1  2^0 ] 
   //// can be put in an small function
   int *power2 ;  // matrix used to transform from binary to decimal
   int prod2 ; // helper variable to fill power2
   prod2 = 1 ;
   power2 = new int [MaxInputs] ;
   for ( int i = 0 ; i < MaxInputs ; i++) {  // MaxInputsStates = 2^(MaxInputs)
       power2[ MaxInputs - i - 1  ] = prod2 ;
       prod2   =  prod2   <<  1 ;
   }
   ///////////////////////////////////////////////////////////////

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
   cout<<"CEdgesWeight:      ";
   for ( int i = 0; i < NumCEdges; i++) {
       cout <<CEdgesWeight[i] <<" ";
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
   cout<<"CNodesWeight:      ";
   for ( int i = 0; i < NumCNodes; i++) {
       cout <<CNodesWeight[i] <<" ";
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
   int * action1 = new int[NumNodeEdges];
   int * action2 = new int[NumNodeEdges];


   if ( sparse_flag  ) {
     
      if ( Noise != "Yes" ) {
          pp = 0.0;
      }

      dec2multinary ( x, NumNodes, p, sparse_s);
      clock_t start_sparse = clock();
      int c = sparse_c;
      int h = sparse_h;
      int BestActionIdx = 0;
      float alpha_over_m = alpha  ;

      for ( int i = 1 ; i < W ; i ++ ) {
         alpha_over_m = alpha * alpha_over_m ;         
      }


      float minQ = 10000000;
      for ( int a_idx  = 0; a_idx < NumActions;  a_idx ++ ) {
          
           float Q_a = RecursiveQ_Noise_LW( NumNodes, NumStates, MaxInputs, sdds_tt, sdds_nv, sdds_varf, sdds_prop, power2, NumNodeEdges, ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, p, BadState, Wi, alpha_over_m,  c, h, x, a_idx, NumActions, v_nodes, v_edges, pp, Lo, L, W, CNodesWeight, CEdgesWeight) ;

           if ( Q_a < minQ ) {
              minQ = Q_a;
              BestActionIdx  = a_idx ;
           }
      }


      // displaying
      if ( Lo > 0 ) { 
          ActionHash(action1,action2,NumNodeEdges,BestActionIdx);
          cout << "state :  " ;
          for ( int i = 0; i < NumNodes; i++ ) {
              cout << x[i] << " ";
          } cout << endl ;
          cout << "action1 : "; 
          for ( int i = 0; i < NumNodes; i++ ) {
              cout << action1[i] << " ";
          } cout << endl ;
          cout << "action2 : ";
          for ( int i = 0; i < NumNodes; i++ ) {
              cout << action2[i] << " ";
          } cout << endl ;
      } 
      else {
          dec2binary( action1,  NumNodeEdges, BestActionIdx ) ;
          for ( int i = 0; i < NumNodes; i++ ) {
              cout << x[i] << " ";
          } cout << endl ;
          cout << " -->  ";
          for ( int i = 0; i < NumNodes; i++ ) {
              cout << action1[i] << " ";
          } cout << endl ; 
      }
      
   }
   else { // exact method 
       ofstream output_policy("full_policy.txt");
       clock_t start_policy = clock();

       float maxdiff = 0.0;
       float * JV = new float[NumStates];
       float * newJV = new float[NumStates];
       int * U = new int[NumStates];       
 
       for ( int i = 0; i < NumStates; i++ ) {
          JV[i] = 0.0;   // What this value ?? shoulnt it be 0?
          newJV[i] = JV[i];
       }
   
       operatorTVi( NumNodes, NumStates, sdds_tt, sdds_nv, sdds_varf, sdds_prop, NumNodeEdges, NumActions,  ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, BadState, Wi,  JV, alpha, newJV ,  U,  v_nodes, v_edges, p, pp, CNodesWeight, CEdgesWeight );

       while( ( normax(JV, newJV, NumStates) > tol ) && ( iter < Nmax ) ) {

           for ( int i=0 ; i < NumStates; i++ )
               JV[i] = newJV[i] ;

           operatorTVi( NumNodes, NumStates, sdds_tt, sdds_nv, sdds_varf, sdds_prop, NumNodeEdges, NumActions,  ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges,BadState, Wi, JV, alpha, newJV ,  U ,  v_nodes, v_edges, p, pp, CNodesWeight, CEdgesWeight );

           iter++;
       }

       for ( int i=0 ; i < NumStates; i++ )
          JV[i] = newJV[i] ;
      
       operatorTVi( NumNodes, NumStates, sdds_tt, sdds_nv, sdds_varf, sdds_prop, NumNodeEdges, NumActions,  ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, BadState, Wi,  JV, alpha, newJV ,  U,  v_nodes, v_edges, p, pp, CNodesWeight, CEdgesWeight );

       clock_t end_policy = clock();

       for ( int State = 0; State < NumStates  ; State++ ) {
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

       delete [] JV ;
       delete [] newJV ;
       delete [] U ;

   }


   delete [] x ;
   delete [] y ;
   delete [] z ;
   delete [] action ;
   delete [] action1 ;
   delete [] action2 ;


   return 0;
}



float RecursiveQ_Noise_LW(int NumNodes, int NumStates, int  MaxInputs,  int **sdds_tt, int *sdds_nv, int** sdds_varf,  float** sdds_prop, int* power2, int NumNodeEdges, int *ActionNodes, int NumCNodes, int *ActionHeads, int *ActionTails, int NumCEdges, int p, int * BadState, float *Wi, float alpha_over_m,  int c, int h, int *x, int ActionIdx, int NumActions,  int *v_nodes, int *v_edges, float pp, int Lo, int L, int W, float *CNodesWeight, float *CEdgesWeight) 
{

    if (h == 0) {
      return 0.0 ;
    }

    float Q = 0.0 ;
    float R = 0.0 ;

    int *y = new int[NumNodes] ;
    int *z = new int [NumNodes];

    int *action1 = new int[NumNodeEdges];
    int *action2 = new int[NumNodeEdges];
   
    if ( Lo > 0 ) { 
         ActionHash( action1, action2, NumNodeEdges, ActionIdx) ;
    }
    else { // se podria reducir el codigo haciendo que ActionHash lea Lo
         dec2binary( action2,  NumNodeEdges, ActionIdx ) ;  
         for (int i = 0; i < NumNodeEdges;  i++ ) { action1[i] = 0 ; }  
    }      

    for (int i = 0 ; i < c ; i++ ) { 
        for ( int i = 0; i < NumNodes; i ++ ) {
            y[i] = x[i];
        }

        for (int v = 0; v < Lo ; v++ ) {
            if ( (pp > 0 ) && ( pp >= distribution(generator) ) )  {
                random_nextstate(z, NumNodes);
            }
            else { nextstate_ia( y, z, NumNodes,  NumStates,  sdds_tt, sdds_nv, sdds_varf, sdds_prop, action1, NumNodeEdges, ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, v_nodes, v_edges, p);
            }

            for ( int i = 0; i < NumNodes; i ++ ) {
                y[i] = z[i];
            }
        }

        for (int v = 0; v < ( L - Lo) ; v++ ) {
           if ( ( pp > 0 ) && ( pp >= distribution(generator) ) ) {
               random_nextstate(z, NumNodes);
           }
           else { nextstate_ia( y, z, NumNodes,  NumStates,  sdds_tt, sdds_nv, sdds_varf, sdds_prop, action2, NumNodeEdges, ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, v_nodes, v_edges, p);
           }

           for ( int i = 0; i < NumNodes; i ++ ) {
              y[i] = z[i];
           }
        }

        for (int v = 0; v < ( W - L) ; v++ ) {
           if ( ( pp > 0 ) && ( pp >= distribution(generator) ) )  {
               random_nextstate(z, NumNodes);
           }
           else { sdds_nextstate( y, z, NumNodes, MaxInputs, sdds_nv, sdds_varf, sdds_tt,  sdds_prop, power2 ) ;
           }

           for ( int i = 0; i < NumNodes; i ++ ) {
               y[i] = z[i];
           }
        }

        float minQ = 10000000;
        for ( int a_idx  = 0; a_idx < NumActions;  a_idx ++ ) {
            float Q_a = RecursiveQ_Noise_LW(NumNodes, NumStates, MaxInputs, sdds_tt, sdds_nv, sdds_varf, sdds_prop, power2, NumNodeEdges, ActionNodes, NumCNodes, ActionHeads, ActionTails, NumCEdges, p, BadState, Wi, alpha_over_m, c, h-1, z, a_idx, NumActions,  v_nodes, v_edges,  pp, Lo,  L,  W, CNodesWeight, CEdgesWeight) ;
            if ( Q_a < minQ ) {
                minQ = Q_a ;
            }
        }
        Q = Q + minQ ;
        R = R + cost_ija_LW(z, BadState, NumNodes, action1, action2,NumNodeEdges, NumCNodes, Wi, CNodesWeight, CEdgesWeight );
       
    }

    delete [] y;
    delete [] z;
    delete [] action1;
    delete [] action2;

    
    return ( R + alpha_over_m*Q ) / c ;
}


//////////////////////////////////////////////////////
/// Comppute the cost( or reward) when the systems
/// transition from a current state "i" to an state "j" 
/// after W steps with L steps drug duration
/// ehwn the action "a" was taken.
/// Inputs 
/// BinState   : Represent the sate of "j", it is an array that store the binary number
///              most significant digit is first (index=0)
/// NumNodes   : Number of nodes (size of BinState)
/// actions    : Vector that stores the states of nodes and edges of
///              of an specific action
/// NumActions : Number of control nodes + N of control edlges
/// NumCNodes  : Number of control nodes
float cost_ija_LW(int *BinState, int *BadState,  int NumNodes, int *action1, int *action2, int NumNodeEdges, int NumCNodes, float * Wi, float *CNodesWeight, float *CEdgesWeight ) {

   ///// Parameters to compute Cost
   float  Cost = 0.0 ;

   //////////////////////////////////////////////// 

   // cost of control nodes
   for ( int i = 0 ; i < NumCNodes; i++ ) {
      Cost = Cost +  CNodesWeight[i] * action1[i] ;
      Cost = Cost +  CNodesWeight[i] * action2[i] ;
   }

   // cost of control edges
   for ( int i = NumCNodes ; i < NumNodeEdges  ; i++ ) {
      Cost = Cost +  CEdgesWeight[ i - NumCNodes] * action1[i] ;
      Cost = Cost +  CEdgesWeight[ i - NumCNodes] * action2[i] ;
   }

   // cost of bad state
   for ( int  i = 0; i < NumNodes; i++){
      if ( BinState[i] != BadState[i]) {
         Cost = Cost + Wi[i] ;
      }
   }   
   
   return Cost ;
}

//////////////////////////////////////////////////////////////
///  Define Operator T
///   x current state ( array of 1 or 0, with a size of NumNodes)
///   pia       : array of 2^n douable elements, this will hold the probabilitites
///   NumNodes  : number of nodes of SDDS
///   MaxInputs : Maximun indegree
void operatorTVi( int NumNodes, int NumStates,  int **sdds_tt, int *sdds_nv, int **sdds_varf,  float ** sdds_prop,   int NumNodeEdges,
                  int NumActions ,  int *ActionNodes, int NumCNodes, int *ActionHeads, int *ActionTails, int NumCEdges, int * BadState,
                  float *Wi, float *JV, float alpha, float *newJV , int * U, int *v_nodes, int *v_edges, int p, float pp, float *CNodesWeight, float *CEdgesWeight )
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
                   NumCNodes, ActionHeads, ActionTails, NumCEdges, BadState, Wi, pia, v_nodes, v_edges, p, pp, CNodesWeight, CEdgesWeight ) ;


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
///   x         : current state ( array of 1 or 0, with a size of NumNodes)
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
               int *ActionTails, int NumCEdges, int * BadState, float *Wi, float* pia, int *v_nodes, int *v_edges, int p, float pp, float *CNodesWeight, float *CEdgesWeight  )
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


   for ( int state  = 0 ; state < NumStates; state++ ) {

       dec2multinary ( y, NumNodes, p, state);
       pia[state] = 1.0;
       Cost_ia[state]= cost_ija( y , BadState,  NumNodes, action, NumNodeEdges , NumCNodes, ActionNodes, Wi, CNodesWeight, CEdgesWeight );


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
/// This function generate a random state randomly.
/// Inputs 
///   z         :  output array 
///   NumNodes  : number of nodes of SDDS
/// Output
///   void  
void random_nextstate(int *z, int NumNodes) {
   for (int i=0 ; i < NumNodes ; i++ ) {
       double r  = distribution(generator); // random generator
       if (  r < 0.5 ) 
           z[i] = 0;    
       else
           z[i] = 1;
   }
   return;
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
/// when the action "a" was taken.
/// Inputs 
/// BinState   : Represent the sate of "j", it is an array that store the binary number
///              most significant digit is first (index=0)
/// NumNodes   : Number of nodes (size of BinState)
/// actions    : Vector that stores the states of nodes and edges of
///              of an specific action
/// NumActions : Number of control nodes + N of control edlges
/// NumCNodes  : Number of control nodes
float cost_ija(int *BinState, int *BadState,  int NumNodes, int *actions, int NumNodeEdges, int NumCNodes, int *ActionNodes, float * Wi, float *CNodesWeight, float *CEdgesWeight ) { 

   float  Cost = 0.0 ;

   // cost of control nodes
   for ( int i = 0 ; i < NumCNodes; i++ ) 
      Cost = Cost + CNodesWeight[i] * actions[i] ;

   // cost of control edges
   for ( int i = NumCNodes ; i < NumNodeEdges  ; i++ ) 
      Cost = Cost + CEdgesWeight[ i - NumCNodes]  * actions[i] ;

   // cost of bad state
   for ( int  i = 0; i < NumNodes; i++){ 
      if ( BinState[i] != BadState[i]) {
         Cost = Cost + Wi[i] ;
      }
   }

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

//////////////////////////////////////////////////////
///Convert decimal index to two vector arrays 
///Inputs  
///action1      : Binary vector that defines active nodes or edges initially
///action2      : Binary vector that defines active nodes or edges after L/2
///NumNodeEdges : Number of bits (actions), size of actions1 and 2
///ActionIndex    : Decimal number to be converted to binary

///output 
///Binstate of action1 and action 2
/////////////////////////////////////////////////////
void ActionHash ( int *action1, int *action2, int NumNodeEdges, int ActionIndex)
{
   int Idx_a1 = ActionIndex / ( NumNodeEdges + 1 );
   int Idx_a2 = ActionIndex % ( NumNodeEdges + 1 );

   for (int i = 0; i < NumNodeEdges; i++) {
      action1[i] = 0;
      action2[i] = 0;
   }		
  
   if ( Idx_a1 > 0 ){
      action1[ Idx_a1 - 1] = 1;
   }
   
   if ( Idx_a2 > 0 ){
      action2[ Idx_a2 - 1] = 1;
   }
}

