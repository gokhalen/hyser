#ifndef HAVE_MESH_HPP
#define HAVE_MESH_HPP
#endif

#ifndef HAVE_UTIL_HPP
#include "../utillib/utils.hpp"
#endif

#ifndef HAVE_ELEMLIB_HPP
#include "../elemlib/elemlib.hpp"
#endif

#ifndef HAVE_FUNCDECLS_H
#include "../../headers/funcdecls.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>


class mesh{

public:
  // long int is too big.
  // int on zodiac does not seem to be restricted to 32728

  // a flag which decides what type of mesh it is
  // set to true, it is a inverse mesh
  // should be set from outside by the user
  bool inverseflag;


  // this flag is specifically built in for linear elasticity
  // in this case lambda and \mu are linked through 
  // a spatially varying poisson's ratio
  bool linklambdamu; 
  

  int **fillin;
  int *fillinno;
  int *row_ptr;  int *row_ptr2; int *row_ptr_mass;
  int *col_ptr;  int *col_ptr2; int *col_ptr_mass; // the column indices for the matrix
  
  int nnz;        // number of nonzeros

  long int elem_in_dir[3];   // elements in direction regular meshes
  long int node_in_dir[3];   // nodes in direction
  
  int nfields; // number of boundary condition sets.
  int which_field; // to store which field is being operated on.
  long int num_elem,num_nodes,num_prop_nodes; 

  // eigen & dispersion stuff
  double fmin, fmax;  // min, max frequency (Hz)
  double rayalpha,raybeta; // rayleigh damping alpha and beta (needed for w2k)
  // num_eigen = number of requested 
  long int num_eigen,num_conv_eigen;  
  //double dspkmin,dspkdel,dspkmax,dspktheta;
  int ndispersionseg;
  int maxdispersionpts;
  int *ndispersionpts;
  double ***ibzcoord;          // (nseg,2[start,end],ndime)
  double **w2karray;          // (nseg,5)

  // arrays to hold eigenvalues at each k point along the dispersion curve
  double ***eigdisperreal;   //  (num_eigen,ndispersionptsm,ndispersionseg )
  double ***eigdisperimag;   //  ( ----------------" --------------------  )
  double ***eigdispererror;       //(num_eigen,ndispersionptsm,ndispersionseg)
  double ***disperkvecs;      //  ( (ndime,maxdispersionpts,ndispersionseg)the k vectors along the boundary    ) 
  double **disperwvecs;      //  ( (maxdispersionpts,ndispersionseg)the k vectors along the boundary    ) 

  int which_disperseg;    // which dispersion segment is current
  int which_disperpt;    // which point in the segment is current
  
  // total number of dirichlet and traction dofs specified
  // dirdofn is the maximum number of dirichlet degrees of freedom over all fields (used for memory allocation only)

  long int dirdofn,tracdofn,perdofn,pdofn; 
  // an array containing the number of NODES on which dirichlet data is specified for this field
  long int *ndirich, *ntrac, *nperiodic, *npforc;
  long int *ndirdofn, *ntracdofn, *nperdofn, *npforcdofn;

  // max number of properties, max dofn per node, number of dimensions
  // nregparm == max number of regularization parameters per field
  int mnode,mdofn,ndime,nprop,nregprop,mregparm,loadits,newtonits; 
  
  // store the current load and newton iteration
  int iloadits;
  int inewtonits;
  int intelligentnewton;
  
  int band;               //bandwidth
  int max_nnz_per_row;    
  int conseqnno;          // constitutive equation number (used by cons. eqn. lib)
  long int *num_equations; //num_equations, is calculated in read_mesh
  long int  num_equations_max;  // max_number of equations
  
  double   **coord;     // coordinate arrays
  long int **conn;      // connectivity array
  long int **nodeid;    // node number to node id 
  long int *elemsattached; // nodes attached to a particular node
  double   **prop;      // material property array (numnodes,nprop)
  double   ***propcons; // bound constraints for nodal property values (numnodes,nprop,2)
  double    **density;   // density
  int      **lmref;     // element type and number of nodes per elem array
  int        ninte;     // number of integration points in the element

 
  
  // dirich and trac hold the dirichlet boundary conditions 
  // the boundary conditions are a function of the load iteration number
  
  // we don't have a problem with the dirichlet boundary conditions 
  // being a function of the load iteration iteration number
  // however, we probably have a problem with the traction being a function 
  // of the load iteration number. Check this. This is interesting, 
  // this needs to verified: should we be taking into account the non-dead loading

  int   *initstrn; // type of initial strain forcing -- used for homogenization
  double ***dirich; // dirich(dirichletdofs,loadits,nfields)
  double ***trac;   // trac(tracdofs,loadits,nfields)
  double ***pforc;  // point force (pdofs,loadits,nfields)
  double ***perbc ; // periodic boundary conditions (pdofs, loadits, nfields)
   
  // finite element arrays, see TJR Hughes(linear, pg 94, dover)
  // Feb 23 2009, added ID tag

  long int ***id;    //INPUT:  global node num, dof number. 
                    //OUTPUT: global eqn number, zero if dirichlet

  long int **ien;   //INPUT: element number, local node number. 
                    //OUTPUT: global node number 

  //do not construct this array
  //long int **lm;    //LM(dofnumber,local_node_num,elem_number)=global eqn_number


  // Feb 23 2009: made isbc a function of the number of displacement fields
  long int ***isbc;  
  long int ***ispointbc;  // if a point force exists at this node
  long int ***isperbc;

  // if a (node,dofn) is a dirichlet boundary condition:
  // stores location in dirich (positive number) 
  // if a (node,dofn) is a traction  boundary condition:
  // stores neg of location in trac
  // else zero

  // A list of traction elements and dirichlet boundary condition elements
  long int *tracelems;
  long int *direlems;

  

  
  //  long int **matcoord;
  //  double *stiffmat;
  
  
  double **rhsvec; //(num_equations,nfields)
  double **solrhsreal; //(num_equations,nfields)
  double **solrhsimag; //(num_equations,nfields)
  
  double ***solvec;       //(num_nodes,mdofn,nfields)
  double **eigenvalreal;  //(num_equations,nfields) - more space then requested
  double **eigenvalimag;  // -----------" " ----------------------------
  double **eigenerror;    // -----------" "--------------------------------

  // eigenvectors are not stored permenently
  // overallocated
  double ***eigenvecreal; //(num_nodes,mdofn,nfields); 
  double ***eigenvecimag; //(num_nodes,mdofn,nfields); 
 
  // array to store the last update to the solution \Delta\b{U} frm newton it.
  double ***lastupdate; //(num_nodes,mdofn,nfields)

  // the enhancements
  // this is a rough estimation of size --works for the pilter taylor element
  // a more general framework to accomodate enhancements is needed
  double ***enhvec; //( num_elem,4,numfields)
  
  double rhsnorm;
  double fdstep;
  
  int fdgradflag;  // flag to test gradient with finite differences
  int taoiterates;      // number of tao (bfgs) iterations
  int ctangflag;   // check consistent tangent is evaluated at every newton it.
  int primalguessflg; // check if load increments used for inverse its > 1
  int propdiscontinflag; // check if the properties are discontinuous
  int firstguessflag;
  int linearflag; // 1 for linear, 0 for nonlinear
  int solvetype;  // 1 for straight-boundary valu problem solve, 2 for eigen, 3 for dispersion curves,
  int which_mass;   // consistent mass or lumped
  int printflag;  // 
  int writeeigenvectorflag; // flag to write eigenvectors 
  
  bool localizelast;
  bool init_primal_goodguess_flag; //this variable decides whether to use good guess for the primal or to use the bad guess. this is set by calc_forward

  // physics 
  int physics;
  
  // movie options
  int movieflag;
  
  // Inversion Stuff
  int numoptvars;
  int inverseiteration;
  int regfunctionaltype;    // chooses type of regfunctional
  int functionaltype;       // chooses type of functional
  int opttype;              // optimization type (0=all params free, 1=lambda is opt var, 2=mu is optvar,
  double poissons;          // possions=poissons ratio to use) to fix the optimizaton parameter
  double lamtomu,mutolam;
  double totalvolume;
  double totalmass;
  double ***measdata;   //(num_nodes,mdofn,nfield)
  double ***fdgradgrad; //(num_nodes,nprop,nfield) 
  double ***dualfield;  //(num_nodes,mdofn,nfield)
  double ***gradient;   //(num_nodes,nprop,nfield)
  

  // a guess of the displacement field. used to speed up the inverse prob
  // (num_nodes,mdofn,nfield)
  double ***gooddispguess; 

  double  **combgrad;   // (num_nodes,nprop) the combined gradient. 
                        // obtained by adding all gradients, from each of the files.
  double  **combfdgrad; // combined finite difference gradient.
  double  **functens;   // (ndime,ndime) functens is the tensor T in |T(u^{m} - u)|

  double *functarray;    // (nfield) to store functional for each field.
  double *regfunctarray; // (nfield) to store the regularization functional 
  //  the regularization functional is independent of the deformation field
  //  because the material properties are the same for every deformation field (unlike disp fields)
  //  however, since the regularization functional is calculated for every deformation field
  //  this is made into an array here
  //  it should be calculated only once....make this change later

  double *normmeasdata;  //norm of the measured data |T(u^{m}|  
  

  // regularization parameter
  double **regparm;    // (nprop) contains the reg param for each prop
  double *fieldfuncweight; // the weight to use for this field
  
  double functional;
  double regfunctional;

  // homogenization stuff
  double ****EHIJKL;    // the homogenized properties
  double ****EHAVIJKL; // average properties
  double ****PMIJKL;    // the target properties
  double ****GIJKL;     // the difference 

  double  **KHIJ;       // the dielectric tensor
  double  ***PZHKIJ;     // the piezoelectric tensor
  double  **KHAVIJ;
  double  ***PZHAVKIJ;  

  double volume;         
  double unitcellvolume; //  volume of the unit cell
  double *coordmin;
  double *coordmax;
  double *meshdim;

  // stiffness matrix
  MATCSR *stiff; 
  MATCSR *stiffimag;
  MATCSR *mass;
  MATCSR *massimag;

  // quadratic eigenproblem
  MATCSR *quadmr;
  MATCSR *quadcreal;
  MATCSR *quadcimag;
  MATCSR *quadkr1;
  MATCSR *quadkr2;
  MATCSR *quadki;

  // a one element stiffness matrix --hack for armero validation purposes
  //double **stiffoneelem;

public:
  mesh();         // default constructor

public:  
  ~mesh();        // destructor deallocates memory

//==============================================================================
  
public:  
  // alphabetically arranged (approximately)
  void accumulatenormmeas(double);   // calculate the norm of the measured data 
  void allocmisc();                  // allocate miscellaeneous memory -- called in readmesh
  void allocate_arrays();    // allocate memory: see also "deallocate_arrays"
  void applybcstosol(double **);      // applies dirichlet bcs to the solution from the solver

  void calculate_optvars();  // calculate number of optimization variables
  void combine_gradients();  // combine gradient vectors into one useful gradient vector
  void combine_fdgradients();// combine finite difference gradient vectors into one vector
  void computeaveprops();
  void computevolumemass();  
  void computenumeqns();     // compute number of equations
  // delete the stiffness, mass and damping matrix
  void delete_stiff();       
  void delete_mass();
  void delete_damping();

  void deallocate_arrays();  // deallocate memory see also "allocate_arrays"

  void getbcs(long int ***,double ***,long int*, int,FILE *);  // get bcs

  void init_primal();           // initalizes primal guess consistent with dirich BC
  bool isenhanced(int);            // returns true if the element has enhancement; false otherwise
 

  void make_ien();          // ien(element num, local node num) = glb node number
  //void make_id();           // id(global node num, dofn) = glb eqn num, 0 if dirichlet
  void make_eqn_num();      // makes equation num, starting at zero
  void make_fillin();       // makes the fillin for the CSR

  bool isdirichlet(long int ,int, int); // true if dirichlet bound else false
  bool istraction(long int ,int, int);  // true if traction  bound else false
  bool isslave(long int ,int, int);  // true if slave node - e.g. in periodic master-slave else false
  
  double getdirichlet(long int, int, int, int); //returns the dirichlet value
  double gettraction(long int, int, int, int); //returns the traction
  long int getmasterperiodic(long int, int, int); // get the master periodic node

  void printrhsnorm();    // prints norm of the rhs

  void read_mesh();       // read mesh and allocate arrays
  void readmeasdata();    // reads measured data from meas.in

  // stores the dual field in the array this->dualfield. (make more general by accepting pointer argument)
  void writeonefieldfunctional();
  void resetprops(); // reset the props (make them consistent) if opttype == 1,or 2
  void storedual();      
  bool issolvedispersion();
  bool issolveeigen();

  void updategoodguess();
  void view_and_store_sol();

  void writecoord();
  void writedualrhs();
  void writedual();
  void writeeigenvalues();
  void writefdgradgrad(); 
  void writecombfdgrad();
  void writefunctional();
  void writeonefieldgradient();
  void writegradient();
  void writemeasnorm();
  void writeproperties();
  void writeehijkl();
  void writeehijkl3d();
  void writeprimalrhs();
  void writeenhancedparm();
  void zerofunctional(); // zeroes "functional" and "functarray"
  void zeroenhanced();   // zeroes the enhanced modes
  void zerogradient();
  void zerorhs();

  // vtk functions
  void initvtk(char*);        // initializes the vtk file with common information
  void writevtk();       // writes stored the mesh output in a VTK format at the end of forward.cpp
  void writeeigenvectorsvtk(int);  
  
  // dispersion stuff
  void computedisperkvecs();  // individual values of k
  void computedisperwvecs();  // individual values of w

  void seteigendisper(int,double,double,double);
  void writedispersion();
  
}; 
/*=========================================================================*/

