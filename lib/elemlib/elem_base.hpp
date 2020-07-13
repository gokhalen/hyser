/*========================================================================
  Nachiket Gokhale gokhalen@bu.edu: Feb 2005 elem_base.cpp 
  Base class for finite elements. Keep the common data and functions here.
  This defines virtual functions that other elements must implement
==========================================================================*/

#include <math.h>

#ifndef HAVE_ELEM_BASE_HPP
#define HAVE_ELEM_BASE_HPP
#endif

#define L2NORM          1
#define H1SEMINORM      2
#define TOTALVARIATION  3
#define PAULNORM        4
#define PENTANORM       5
#define PENTANORM2      6
#define PENTANORM3      7
#define PENTAREG1       8

// common to elem_base and globals - must be consistent with globals
#define BUILD_REAL_STIFFNESS        1 
#define BUILD_COMPLEX_STIFFNESS_K2W 2
#define BUILD_COMPLEX_STIFFNESS_W2K 3

// common to elem_base and globals - must be consistent
#define BUILD_CONSISTENT_MASS 1
#define BUILD_LUMPED_MASS 2 


#ifndef HAVE_SHAPESTRUCT_HPP
#include "../datastruct/shapestruct.hpp"
#endif


class elem_base {

public:
  //ELEMENT ID

  int nprop;          // number of properties at each node
  int nregprop;       // regularization props;
  int mregparm;       // max number of regularization parameters
  int mdofn;          // max number of degrees of freedom at each node
  int *dofnatnode;    // (local node number, number of degrees of freedom at that node) 

  int npropgrad;      // number of properties to take the gradient with respect to
  int *wch_grad_prps; // A 1d array of properties to take the gradient wrt.

  // COORDINATE DATA
 

  //LOCAL FINITE ELEMENT ARRAYS
  long int  **ideqn;       // INPUT: (local node number, dofn) OUTPUT: (global eqn number)
  long int *idnode;       // INPUT: (local node number)       OUTPUT: (global node number)
  long int   **isbc;       // INPUT: (local node number, dofn) OUTPUT: (1==dirichlet,0==none,-1
  
  
  //CALCULATED QUANTITIES
  double     **ecmass;                    // element consistent mass matrix
  double     **estiff;                    // elemental stiffness matrix
  double     **estiffimag;                // imaginary part of the stiffness matrix
  double     **K11, **K12,**K22, **K22INV;// element matrices
  double      *erhs;                      // elemental right hand side
  double      *erhsenh;                   // erhs enhanced
  double      *erhsnc;                    // not condensed
  double       functional;                // functional 
  double       regfunctional;             // regularization functional 
  double     **egrad;                     // gradient vector
  double     **stiffkr1,**stiffkr2,**stiffki; // stiffness matrices for w2k formulation
  
  //DATA
  double    **primal;      // primal field -- the main qty we solve for
  double    **lastupdate;  // last update to the primal field
  double    *enh;          // enhancements;
  double    *enhinc;       // increment in the enhancement
  double    **measdata;    // measured data
  double    **prop;        // elemental properties (\lambda and \mu or more)
  double    **dualfield;   // dualfield
  double    **gradient; 
  double     *density;
  double    *propvec;     // (nprop) a vector containes properties at each point - deprecated
  double    *kvec;        // (ndime) wave vector 
  double    wdisp;        // the circular frequency at which at w-k calculation has to be carried out
    

  double    ***primalwhole;  // primal data corresponding to all displacement fields - used in homogenization
  double    ***dualwhole;    // dual   data  --------"---------  dual         ----------"--------------------

  // ELASTICITY QUANTITIES
  double **CIJ;
  double **krchstrs;
  double ****sptan;
  

  // INITIAL STRAIN BODY FORCE FOR HOMOGENIZATION
  int initstrntype;
  double *initstrn;
  double ****BASE_EHIJKL;
  double ****BASE_EHAVIJKL;
  double ****BASE_PMIJKL;
  double ****BASE_GIJKL;
  double ****BASE_GHATIJKL;
  double ****BASE_CIJKL;

  double ***BASE_PZHKIJ;
  double ***BASE_PZHAVKIJ;

  double **BASE_KHIJ;
  double **BASE_KHAVIJ;
  
  // 
  double elemvolume;
  double elemmass;
 

  // USEFUL GLOBAL DATA
  // this is put in purely for convenience
  // it shouldn't be, because it doesn't vary at the element level
  // this is used to evaluate the homogenized Eijkl tensor a la Hassani & Hinton 1998 Computers and Structures Vol 69 (1998).
  double unitcellvolume;    
  int nfields;
  

  //BOUDARY DATA

  double **dirdata;        // dirichlet data (numnodes,mdofn)
  double **tracdata;       // traction data  (numnodes,mdofn)

  // constitutive equation number
  int conseqnno;

  // element option number
  int elemopt;
  
  // inverse problem parameters
  int opttype;

  // the tensor which modifies the functional 
  double **functens;
  
  // the transformed primal and measured data
  double **Tprim;  // (numnodes, mdofn)
  double **Tmeas;  // (numnodes, mdofn)

  // (nprop) allocation done by this class, should match size
  double **regparm;
  
  double funcweight; // weight for this field

  // boolean 

  // is it nonlinear or not: probably unused. .
  // should be used to gain speed in the global scheme of things
  bool nonlinflag;
  bool condense_matrix; 
  int propdiscontinflag;


  // miscellaneous
  double rhsnorm;  
  double normmeasdata;

  // inverse problem options
  int which_functional;      // set the functional to minimize
  int which_regularization;  // set the regularization to use
 
  shapestruct *shpstrct;   
  


  // FINITE ELEMENT STUFF
  
   void setmatprops();            // set material properties
   void setfemarrays();           // set femarrays: local coordinates

  
   
  // warning: virtual functions must be implemented
  //This file is the header and there are 
  virtual void elemstiff(int);      // calculate elemental stiffness 
  virtual void elemcmass(int); // compute a consistent mass matrix 
  virtual void elemrhs();        // calculate right hand side
  virtual void calcenh();        // calculate the enhanced parameters 
  virtual void eleminitstrn();   // initial strain body forcing

  
   

  // inverse problem functions

  virtual void addreggradient();     // adds the contrib of reg to grad (common)
  virtual void calcdual();      // calculate forcing for dual
  virtual void calcfunc();      // calculate functional
  virtual void calcgradient();  // calculate gradient
  virtual void calcnormmeas();  // calculate norm of measured field

  // calculates the gradient of the primal field
  // a variational problem is solved to find the nodal strains
  //virtual void calcstrain();    
  


  //INTEGRATION POINTS          // gets integration points 
  virtual void getinteg();
  
  //SHAPE FUNCTIONS
  virtual void getshape();

  // memory functions
  virtual void memalloc();
  virtual void memdealloc();

  // utility functions
  virtual void calculate_Tmeas_Tprim();
   
  elem_base();
  virtual ~elem_base();   // the virtual destructor is extended not overridden
   
  // DATA SETTING: 
  virtual void setcoord(double, int, int);
  virtual void setprop (double, int, int);
  virtual void setdensity(double,int);
  virtual void setideqn(int,int,long int);
  virtual void setidnode(int,long int);
  virtual void setisbc(int,int,int);
  virtual void setdirdata (double, int, int);
  virtual void settracdata(double, int, int);  
  virtual void setmeasdata(double, int,int);
  virtual void setkvec(double,int);
  virtual void setwdisp(double);
  
  void setdualfield(double,int,int);
  void setprimal(double vval, int inode, int idofn);
  void setlastupdate(double, int, int);

  virtual void homogenizehassanihinton();
  virtual void computevolumemass();
  virtual void computeaveprops();
  
  // miscellaneous utility functions 
  virtual void viewid();

};
