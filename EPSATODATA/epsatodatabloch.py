# This script makes an element.tec file from tape5.elem
# 
import sys;
import math;
import numpy;

if len ( sys.argv ) < 4 : 
   print "Incorrect Usage. Three arguments required: Inputfile Outputfile unitcelltype \n"
   print "Exiting ...\n"
   sys.exit(1);

infile       =str(sys.argv[1]);
outfile      =str(sys.argv[2]);
inputfile    =open(infile,'r');
outputfile   =open(outfile,'w');
unitcelltype =str(sys.argv[3]);

if ( unitcelltype == "hexagon" ):
   print " Hexagon Detected\n";
   nibzseg = 5;
   if ( len (sys.argv) != 7 ):
      print "Hexagon requires: LL, HH, theta (degrees)\n"
      print "Usage: InputFile OutputFile hexagon LL HH theta(degrees)"
      print "..exiting..\n";
      sys.exit(1);
   else:
      hexll = float(sys.argv[4]);
      hexhh = float(sys.argv[5]);
      hexthetarad = (math.pi/180.0)*float(sys.argv[6]);
      print "Hexagon unit cell dimensions: LL= " + str(hexll) + " HH= " + str(hexhh) + " hextheta= " + str((180.0/math.pi)*hexthetarad) + str("\n");

elif ( unitcelltype == "rectangle"):
   nibzseg = 4;
   print " Rectangle Detected\n ";
   if ( len(sys.argv) != 6):
      print "Rectangle requires: LL, HH\n";
      print "Usage: InputFile OutputFile rectangle LL HH ";
      print "..exiting..\n";
      sys.exit(1);
   else:
      rectll=float(sys.argv[4]);
      recthh=float(sys.argv[5]);
      print "Rectangle unit cell dimensions: LL= " + str(rectll) + " HH= " + str(recthh) + str("\n");

elif ( unitcelltype == "ribs"):
   nibzseg = 1;
   print "Ribs Detected\n";
   if ( len(sys.argv) !=6 ): 
      print "Ribs Require: kmin,kmax\n";
      print "Usage: Inputfile Outputfile ribs LL HH";
      print "....exiting ...\n";
   else:
      ribkmin = float(sys.argv[4]);
      ribkmax = float(sys.argv[5]);
      print " kmin = " + str(ribkmin) +  " kmax = " + str(ribkmax) + "\n";

else:
   print " Unknown unit cell type ... exiting ...\n"
   sys.exit(1);
   
listfile   = inputfile.readlines();
# set up dictionaries for python and 
elems = {};
nodes = {};
prop  = {};  
# prop is indexed with elem, not with nodes
# XX means normal to the X-axis; YY means normal to the Y-axis, ZZ means norma
suportside   = {};
suporttop    = {};
suportcorner = {};
matepsa = {};
perbc  = {};

# the first field occuring in the list is the separator
separator = ':';
for ii in range(14):
   # throw away the first 14 lines separator
   tempstring = listfile.pop(0);

# throw away the first separator
# initialize FE data
nfields    =  1;
ndime      =  2;
mnode      =  4;
mdofn      =  2;
nprop      =  3;
mregparm   =  2;
num_nodes  =  int(0);
num_elem   =  int(0);
ndispersionpts = int(32);
ninteriorseg   = int(8);
#
if (( unitcelltype == "rectangle" ) or ( unitcelltype == "hexagon")):
   ndispersionseg = ninteriorseg*(nibzseg-2) + nibzseg;
elif (unitcelltype == "ribs"):
   ndispersionseg = nibzseg;

ndirich    =  int(0);
dirdofn    =  int(0);
ntrac      =  int(0);
tracdofn   =  int(0);
perdofn    =  int(0);
nperiodic  =  int(0);
pdofn      =  int(0);
npforc     =  int(0);
npforcdofn =  int(0);
band       =  int(40);
loadits    =  int(1);
newtonits  =  int(1);
intelligentnewton = int(1);
conseqno   = int(4);

# phase-I parse input 
# parse input file; compute quantitites of interest - ndirich etc; write output file
while ( len(listfile) > 0 ):
	listline  = listfile.pop(0);
	listline  = listline.lstrip();
        splitline = listline.split();
        if ( len(splitline) > 0 ) : 
           if splitline[0] == separator:
              temp=int(0)  ;# do nothing
           elif splitline[0] == 'matlq':
              #print splitline
              if ( int(splitline[2]) != 1 ):
                 print "Unsupported Material Type Found"
              else:
                 # EE, nu, rho
                 matepsa[int(splitline[1])] = [ float(splitline[3]), float(splitline[4]), float(splitline[5]) ];
                 #print matepsa[int(splitline[1])];

           elif splitline[0] == 'node':
              inode = int(splitline[1]);
              num_nodes=num_nodes+1;
              nodes[inode] = [ float(splitline[2]), float(splitline[3]) ]
           elif splitline[0] == 'quad4':
              num_elem = num_elem+1;
              ielem        = int(splitline[1]);
              prop[ielem]  = int(splitline[2]);
              elems[ielem] = [ int(splitline[3]),int(splitline[4]),int(splitline[5]),int(splitline[6])];
           elif splitline[0] == 'suport':
              nn    = len(splitline);
              inode = int(splitline[2]);
              if splitline[nn-1] == 'Sides':
                 dirdofn = dirdofn+1;   ndirich=ndirich+1; suportside[inode]=int(1);
              elif splitline[nn-1] == 'Top':
                 dirdofn = dirdofn+1;   ndirich=ndirich+1; suporttop[inode]=int(1);
              elif splitline[nn-1] == 'Corner':
                 dirdofn = dirdofn+2;   ndirich=ndirich+1; suportcorner[inode]=int(1);
              else:
                 print "Unrecognized label in suport -> should be Side Top "
           elif splitline[0] == 'endsht':
              print "endsht"+"\n";
           elif  ( (splitline[0] == 'momcon') or (splitline[0] == 'pincon' ) ):
              masternode = int(splitline[2]);  slavenode = int(splitline[7]);
              perbc[masternode] = slavenode;
              perdofn = perdofn + mdofn;
              nperiodic = nperiodic + 1;
           else:
              print "Unrecognized symbol: " + splitline[0] + " In input file"
              
print " ""Nodes: " + str(num_nodes) + " Elems " + str(num_elem) + " suport " + str(ndirich)

# end parse input
# begin computing quantities of interest
# end computing quantities of interest


outputfile.write("TITLE: Homogenization Data File written by EPSATODATA\n");
outstring = str(nfields) + " " + str(num_elem) + " " + str(num_nodes) + " " + str(ndime) + " " + str(mnode) + " " + str(mdofn) + " " + str(nprop) + " " + str(mregparm) + " " + str(ndispersionseg) + " " +  str(dirdofn) + " " + str(tracdofn) + " " + str(perdofn) + " " + str(pdofn) + " " + str(band) + " " + str(loadits) + " " + str(newtonits) + " " + str(intelligentnewton) + " " + str(conseqno)  + "\n";
#
outputfile.write(outstring + "\n");
outstring =  str(ndirich) + " " + str(dirdofn) + " " +  str(ntrac) + " " +str(tracdofn) + " " + str(nperiodic) + " " + str(perdofn) + " " + str(npforc) + " " + str(npforcdofn) + "\n";
#
outputfile.write("ndirich\n");
for ifield in range(nfields):
   outputfile.write(outstring);

outstring = "\n options:fggradflag fdstep fevals ctangflg resusesol firstgss \n 0 0 10 1 0 0 \n\n linearflag printflag propdiscontin  \n 1 1 1"
outstring = outstring + "\n\n solvetype \n 5 1 100 \n"
outputfile.write(outstring);

if ( unitcelltype == "hexagon" ):
   print "Hexagon in solvetype"
   ct = math.cos(hexthetarad);
   st = math.sin(hexthetarad);
   ll = hexll;
   hh = hexhh;
   pi = math.pi;

   lphct = ll + (hh*ct);
   denom = 2*hh*st*lphct;

   print "ll = " + str(ll) + " hh = " + str(hh) + " st = " + str(st) + " ct= " + str(ct) + " lphct= " + str(lphct);

   BB1 = (2*pi/denom)*numpy.array([ -lphct,hh*st]);
   BB2 = (2*pi/denom)*numpy.array([ lphct,hh*st]);
   BB1DOTBB2 = numpy.dot(BB1,BB2);
   BB1DOTBB1 = numpy.dot(BB1,BB1);
   BB2DOTBB2 = numpy.dot(BB2,BB2);

   OO = numpy.array([0.0, 0.0]);

   if ( BB1DOTBB2 <= 0.0 ):
      print "Negative BB1DOTBB2"
      AA = numpy.array([BB1DOTBB1*hh*st/(2*pi),0]);
      BB = numpy.array([-BB1[0]/2.0,BB1[1]/2.0]);
      CC = numpy.array([-BB1DOTBB2*hh*st/(2*pi),BB1[1]]);
      DD = numpy.array([0.0,BB1[1]]);
   else:
      print "Positive BB1DOTBB2"
      AA = numpy.array([BB2[0],0.0]);
      BB = numpy.array([BB2[0], (lphct/(2*pi))*BB1DOTBB2]);
      CC = numpy.array([BB2[0]/2.0, BB2[1]/2.0]);
      DD = numpy.array([0.0,(lphct/(2*pi))*BB2DOTBB2]);

   print "BB1 = " + str(BB1) + " BB2= " + str(BB2) + "\n"
   print "BB1DOTBB1 = " + str(BB1DOTBB1) + " BB1DOTBB2 = " + str(BB1DOTBB2)
   
   outputfile.write(str(OO[0]) + " " + str(OO[1]) + " " + str(AA[0]) + " " + str(AA[1]) + " " +str(ndispersionpts)+"\n");
   outputfile.write(str(AA[0]) + " " + str(AA[1]) + " " + str(BB[0]) + " " + str(BB[1]) + " " +str(ndispersionpts)+"\n");
   outputfile.write(str(BB[0]) + " " + str(BB[1]) + " " + str(CC[0]) + " " + str(CC[1]) + " " +str(ndispersionpts)+"\n");
   outputfile.write(str(CC[0]) + " " + str(CC[1]) + " " + str(DD[0]) + " " + str(DD[1]) + " " +str(ndispersionpts)+"\n");
   outputfile.write(str(DD[0]) + " " + str(DD[1]) + " " + str(OO[0]) + " " + str(OO[1]) + " " +str(ndispersionpts)+"\n");

   DAB = (BB - AA)/(ninteriorseg+1);
   DBC = (CC - BB)/(ninteriorseg+1);
   DCD = (DD - CC)/(ninteriorseg+1);

   PPAB = AA; PPBC = BB; PPCD = CC;

   for ipt in range(ninteriorseg):
      PPAB = PPAB + DAB;
      PPBC = PPBC + DBC;
      PPCD = PPCD + DCD;
      outputfile.write(str(OO[0]) + " " + str(OO[1]) + " " + str(PPAB[0]) + " " + str(PPAB[1]) + " " +str(ndispersionpts)+"\n");
      outputfile.write(str(OO[0]) + " " + str(OO[1]) + " " + str(PPBC[0]) + " " + str(PPBC[1]) + " " +str(ndispersionpts)+"\n");
      outputfile.write(str(OO[0]) + " " + str(OO[1]) + " " + str(PPCD[0]) + " " + str(PPCD[1]) + " " +str(ndispersionpts)+"\n");


elif ( unitcelltype == "rectangle"):
   print "Rectangle in solvetype"
   # O = (0,0), A = (0,2pi/hh),  B = (2pi/ll,2pi/hh),  C=(2pi/ll,0)

   OO  = numpy.array([0.0, 0.0]); 
   AA  = numpy.array([0.0, (2.0*math.pi)/recthh]); 
   BB  = numpy.array([(2.0*math.pi)/rectll,(2.0*math.pi)/recthh]); 
   CC  = numpy.array([(2.0*math.pi)/rectll,0.0]);

   outputfile.write(str(OO[0]) + " " + str(OO[1]) + " " + str(AA[0]) + " " + str(AA[1]) + " " +str(ndispersionpts)+"\n");
   outputfile.write(str(AA[0]) + " " + str(AA[1]) + " " + str(BB[0]) + " " + str(BB[1]) + " " +str(ndispersionpts)+"\n");
   outputfile.write(str(BB[0]) + " " + str(BB[1]) + " " + str(CC[0]) + " " + str(CC[1]) + " " +str(ndispersionpts)+"\n");
   outputfile.write(str(CC[0]) + " " + str(CC[1]) + " " + str(OO[0]) + " " + str(OO[1]) + " " +str(ndispersionpts)+"\n");
   
   DAB = (BB - AA)/(ninteriorseg+1);
   DBC = (CC - BB)/(ninteriorseg+1);

   for ipt in range(ninteriorseg):
      PPAB = AA + DAB;
      PPBC = BB + DBC;
      outputfile.write(str(OO[0]) + " " + str(OO[1]) + " " + str(PPAB[0]) + " " + str(PPAB[1]) + " " +str(ndispersionpts)+"\n");
      outputfile.write(str(OO[0]) + " " + str(OO[1]) + " " + str(PPBC[0]) + " " + str(PPBC[1]) + " " +str(ndispersionpts)+"\n");
elif ( unitcelltype == "ribs" ):
   print "Ribs in Solvetype";
#   print "kmin= "  + str(ribkmin) + " 0.0 " +  "kmax= " + str(ribkmax) + " 0.0 " + "\n";
   outputfile.write(str(ribkmin) + " 0.0 " + str(ribkmax) + " 0.0 " + " " + str(ndispersionpts) + "\n");

outstring = "\n" + "opttype\n 2 0.34\n\n movie\n 1\n\n elem \n"
outputfile.write(outstring);
outstring = "1 " + str(num_elem) + " 202 4 1 0\n\n boundingvolume \n 1.0 \n\n coord\n"
outputfile.write(outstring);
for icoord in range(num_nodes):
   coord = nodes[icoord+1];
   outputfile.write(str(icoord+1) + " " + str(coord[0]) + " " + str(coord[1]) + "\n" );


outputfile.write("\n conn\n");
for ielem in range(num_elem):
    conn = elems[ielem+1];
    outputfile.write(str(ielem+1) + " " + str(conn[0]) + " " + str(conn[1]) + " " + str(conn[2]) + " " + str(conn[3]) + " 0\n" );

outputfile.write("\n prop \n");
for ielem in range(num_elem):
   #lam  = 0.16568; mu=0.000188017;
   #if ( prop[ielem+1] == 2 ):
   # Aluminum 
   # EE = 69; nu = 0.33; rho = 2700
   #lam = 50.2810; mu = 25.902;
   # Brass
   # EE = 100; nu = 0.4  
   # Silicon 
   # EE = 112.8; nu = 0.28  
   # Steel
   # EE = 210; nu = 0.30
   # Lead
   #EE = 16; nu = 0.44; rho = 11340 
   #- matweb
   #lam = 55.8807; mu = 43.9063
   if ( prop[ielem+1] in matepsa ) :
      temp  = matepsa[prop[ielem+1]];
      EE = temp[0];  nu=temp[1]; rho=temp[2];
      outputfile.write( str(ielem+1) + " 2 " + str(EE) + " " + str(nu) + " " + str(rho) + "\n");
   else:
      print "Material " + str(prop[ielem]) + " used in quad not defined"
   
    
outputfile.write("\n propcons \n");
for ielem in range(num_elem):
   outputfile.write( str(ielem+1) + "\n 0.001 1729.0 \n 0.001 1000\n" );

#outputfile.write("\n prop \n");
#for icoord in range(num_nodes):
#    outputfile.write( str(icoord+1) + " 3.14152 3.14152 \n");     

#outputfile.write("\n propcons \n");
#for icoord in range(num_nodes):
#    outputfile.write( str(icoord+1) + "\n 0.001 1729.0 \n 0.001 1000\n");     
 
outputfile.write("\n dirich\n");
for ifield in range(nfields):
    idx=1; idy=0;
    if ( ifield == 2 ) :
       idx=0; idy=1;

    if ( len(suportside) > 0):
       for iside in suportside:
          outputfile.write(str(iside) + "\n" + str(idx) + " 0 \n" + str(idy) + " 0\n"); 
          
    if ( len(suporttop) > 0):
       for itop in suporttop:  
# note flip in idx and idy 
          outputfile.write(str(itop) + "\n" + str(idy) + " 0 \n" + str(idx) + " 0\n"); 

    if ( len(suportcorner) > 0):
       for icorner in suportcorner:
          outputfile.write(str(icorner) + "\n"  + " 1 0 \n" + " 1 0\n"); 
        
    outputfile.write("0\n");



# write initstrn
outputfile.write("\n initstrn \n");
for ifield in range(nfields):
    outputfile.write(str(ifield+1) + "\n");

outputfile.write("0\n\n trac\n");
for ifield in range(nfields):
    outputfile.write(str(0) + "\n");


# write periodic boundary conditions
outputfile.write("\nperiodic \n");

for ifield in range(nfields):
   for masternode in perbc:
      outputfile.write(str(masternode) + "\n" + str(perbc[masternode]) + " 1 1 " +  "\n");
   outputfile.write("0\n");

outputfile.write("\n");
outputfile.write("pdofn\n");
for ifield in range(nfields):
   outputfile.write("0\n");
outputfile.write("\n ninte\n 3\n\n goodguess\n");

for ifield in range(nfields):
   for inode in range(num_nodes):
      outputfile.write("0 0\n");
   outputfile.write("0\n");

outputfile.write("\n functionaltype \n 5 2 \n\n regparm \n 0 0 \n 0  0\n\n functens\n 1 0 0 1\n\n");
outputfile.write("funcweight\n");
for ifield in range(nfields):
   outputfile.write("1 ");

outputfile.write("\n\n Pentatarget - D11 D12 D22 D66\n 1 1 1 1 \n\n end");
# close outputfile
outputfile.close();
