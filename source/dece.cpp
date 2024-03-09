/******************************************************************************/
/**                                                                          **/
/**     DeCE : The Descriptive Correction of ENDF-6 Format Code              **/
/**                                          1.2.6 (Garnet)   Oct.  2022     **/
/**                                                            T. Kawano     **/
/**                                       Los Alamos National Laboratory     **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cstring>
#include <unistd.h>

using namespace std;

#include "dece.h"
#include "command.h"
#include "global.h"
#include "terminate.h"

static string version  = "1.2.6 Garnet (Oct. 2022)";
static bool   verbflag = false;
static bool   justquit = false;
static bool   filescan = false;
static int    newsec   = 0;
static ENDF   *lib[MAX_SECTION];

static void DeceMain          (string, string, ENDFDict *);
static void DeceStoreData     (ENDFDict *, ifstream *);
static void DeceReadMonitor   (const int, const int, const int, const int, const int, const int, const int);
static void DeceHelp          (void);
static void DeceHelpENDF      (void);
static void DeceFreeMemory    (void);
static void DeceBanner        (void);

#undef PeekObject // for debugging

/**********************************************************/
/*      Global Parameters                                 */
/**********************************************************/
#define DECE_TOPLEVEL
GlobalOption  opt;
string        tempfile = "DECETempfile.dat";
ostringstream message;

/*** defined in command.cpp */
extern CLine cmd;


/**********************************************************/
/*      DeCE Main                                         */
/**********************************************************/
int main(int argc, char *argv[])
{
  cout.setf(ios::scientific, ios::floatfield);
  cerr.setf(ios::scientific, ios::floatfield);

  ENDFDict dict;
  int      p, mfin = 0, mtin = 0;
  double   ein = 0.0;
  string   libname_in = "",  libname_out = "";
  bool     reconr = false;

  /*** command line options */
  while((p = getopt(argc,argv,"o:f:t:e:rqsnvhm")) != -1){
    switch(p){
    case 'o':  libname_out = optarg;   break;
    case 'f':  mfin = atoi(optarg);    break;
    case 't':  mtin = atoi(optarg);    break;
    case 'e':  ein  = atof(optarg);    break;
    case 'r':  reconr   = true;        break;
    case 'q':  justquit = true;        break;
    case 's':  filescan = true;
               justquit = true;        break;
    case 'n':  opt.LineNumber = true;
               ENDFPrintLineNumber(opt.LineNumber); break;
    case 'v':  verbflag = true;        break;
    case 'h':  DeceHelp();             break;
    case 'm':  DeceHelpENDF();         break;
    default:                           break;
    }
  }

  /*** check if ENDF file is given */
  if(optind < argc) libname_in = argv[optind];
  if(libname_in == ""){ message << "ENDF-6 formattted file not given"; TerminateCode("main"); }
  if(libname_in == libname_out){ message << "same in/out file names"; TerminateCode("main"); }

  /*** allow to omit when MF=2 (resonance) */
  if((mfin == 2) || (mfin == 32)) mtin = 151;

  /*** tabulate data from one section specified by -f and -t command line options */
  if((mfin > 0) && (mtin > 0)){
    ifstream  fpin;
    fpin.open(&libname_in[0]); if(!fpin){ message << "ENDF file " << libname_in << " cannot open"; TerminateCode("main"); }
    DeceCheckMT(mtin);
    if(ein > 0.0) DeceDataPoint(&fpin,mfin,mtin,ein);
    else          DeceFileToTable(&fpin,mfin,mtin);
    fpin.close();
  }

  else if( (mfin == 0) && (mtin == 0) ){
    /*** scan the file, and store all the MF,MT sections */
    int ctl = ENDFScanLibrary(libname_in,&dict);
    if     (ctl == -1){ message << "ENDF file " << libname_in << " cannot open"; TerminateCode("main"); }
    else if(ctl == -2){ message << "too many sections"; TerminateCode("main"); }

    /*** reconstruct pointwise cross sections from resonances */
    if(reconr){
      ifstream  fpin;
      fpin.open(&libname_in[0]);
      gfrScanThermal(&fpin,&dict,ein);
      fpin.close();
    }
    /*** enter the main part */
    else{
      DeceMain(libname_in,libname_out,&dict);
    }
  }
  else{ message << "MF or MT number not provided"; TerminateCode("main"); }

  return(0);
}


/**********************************************************/
/*      Main Loop for ENDF Scanner / Reformatter          */
/**********************************************************/
void DeceMain(string libin, string libout, ENDFDict *dict)
{
  ifstream  fpin;
  string    ope,arg;

  fpin.open(&libin[0]);

  /*** read ENDF tape in lib[], determine resonance boundary */
  DeceStoreData(dict,&fpin);


  /*** not in interactive mode */
  if(justquit){
    /*** just print file contents */
    if(filescan){
      DeceShowHeaders(dict);
      DeceScanIndex(dict);
    }
  }
  /*** for all operation commands */
  else{
    do{
      if( CmdFgetOneline() < 0 ) break;
      ope = CmdExtractArgument();

      /*** END: terminate code */
      if( (ope == "end") || (ope == "quit") || (ope == "exit") ) break;

      /*** perform each operation */
      else DeceOperation(dict,lib,&fpin);

    }while(!cin.eof());
  }


  /*** put all sections in a temporal file */
  if(libout != ""){
    /*** redirect STDOUT to a tempfile, if output file is given */
    streambuf *save = cout.rdbuf();
    ofstream  fpout;

    fpout.open(&tempfile[0]);
    cout.rdbuf(fpout.rdbuf());

    DeceOutput(&fpin,dict,lib);

    fpout.close();
    cout.rdbuf(save);

    /*** renumber, fix dictionary, etc */
    DeceRenumber(tempfile,libout,dict);

    /*** remove temoral file */
    remove(&tempfile[0]);
  }

  fpin.close();
  DeceFreeMemory();
}


/**********************************************************/
/*      Read Each MF                                      */
/**********************************************************/
void DeceStoreData(ENDFDict *dict, ifstream *fp)
{
  for(int i=0 ; i<dict->getSEC() ; i++){

    /*** skip text data */
    if( (dict->mf[i] == 1) && (dict->mt[i] == 451) ) continue;

    /*** instance object, and read ENDF data */
    lib[newsec] = new ENDF();
    ENDFRead(fp,lib[newsec],dict->mf[i],dict->mt[i]);

    /*** find resonance boundary energies */
    if(dict->mf[i] == 2) ENDFMF2boundary(dict,lib[newsec]);

    DeceReadMonitor(lib[newsec]->getENDFmat(),dict->mf[i],dict->mt[i],newsec,lib[newsec]->getPOS(),lib[newsec]->getNI(),lib[newsec]->getNX());
#ifdef PeekObject
    ENDFLibPeek(lib[newsec]);
#endif

    dict->setID(i,newsec++);

    if(newsec >= MAX_SECTION){ message << "Too many sections, " << newsec;  TerminateCode("DeceStoreData"); }
  }
}


/**********************************************************/
/*      Data Read Monitor                                 */
/**********************************************************/
void DeceReadMonitor(const int mat, const int mf, const int mt, const int sec, const int n, const int ni, const int nx)
{
  message << "MAT:" << setw(5) << mat;
  message << " MF:" << setw(3) << mf;
  message << " MT:" << setw(4) << mt;
  message << "  assigned for Section " << setw(4) << sec;
  message << " sub-blocks " << setw(4) << n;
  message << " : int data " << setw(8) << ni;
  message << " : dbl data " << setw(8) << nx;
  Notice("DeceReadMonitor");
}


/**********************************************************/
/*      Check If Valid MT                                 */
/**********************************************************/
void DeceCheckMT(int mt)
{
  if( mt <= 0 || 1000 <= mt ){ message << "invalid MT number, " << mt; TerminateCode("DeceCheckMT"); }
}


/**********************************************************/
/*      Allocate New Section                              */
/**********************************************************/
void DeceCreateLib(ENDFDict *dict, int mf, int mt)
{
  DeceCheckMT(mt);

  /*** if already exists */
  int k = dict->getID(mf,mt);
  if( k >= 0 ){
    lib[k]->resetPOS();
    message.str("");
    message << "MF" << mf << ":MT" << mt << " exists, reuse it. assigned section " << dict->getID(mf,mt);
    Notice("DeceCreateLib");
    return;
  }

  try{
    lib[newsec] = new ENDF;
  }
  catch(bad_alloc &e){
    message << "memory allocation error"; TerminateCode("DeceCreateLib");
  }

  lib[newsec]->setENDFmat(dict->getMAT());
  lib[newsec]->setENDFmf(mf);
  lib[newsec]->setENDFmt(mt);

  int i = dict->scanDict(mf,mt);
  if(i < 0){
    if(dict->addDict(mf,mt,0,newsec)){ message << "Too many sections"; TerminateCode("DeceCreateLib"); }
  }
  else{
    dict->setID(i,newsec);
  }

  message << "MF" << mf << ":MT" << mt << " assigned for section " << newsec;
  Notice("DeceCreateLib");

  newsec++;
}


/**********************************************************/
/*     Free Allocated Memory                              */
/**********************************************************/
void  DeceFreeMemory()
{
  for(int i=0 ; i<newsec ; i++) delete lib[i];
}


/**********************************************************/
/*      Help                                              */
/**********************************************************/
void DeceHelp()
{
  DeceBanner();

  cout <<
    "Command Line Mode\n"
    "  % dece : ENDF-6 format read / write program\n"
    "  % dece -o ENDF_out.DAT ENDF_in.DAT\n"
    "      ENDF_in.DAT  : input ENDF-6 formatted file\n"
    "                     output only specified sections\n"
    "      ENDF_out.DAT : output ENDF-6 formatted file\n"
    "                     all sections will be reformatted\n"
    "  % dece -f MF -t MT ENDF_in.DAT\n"
    "             MF MT : MF and MT numbers to be extracted\n"
    "                     convert the section into a tabulated format\n"
    "  % dece -f MF -t MT -e ELAB ENDF_in.DAT\n"
    "              ELAB : laboratory energy\n"
    "                     print cross section at ELAB\n"
    "  % dece -r ENDF_in.DAT\n"
    "                     print thermal cross sections\n"
    "  % dece -s ENDF_in.DAT\n"
    "                     print data contents\n"
    "Interactive Mode\n"
    "  % dece ENDF_in.DAT\n"
    "  command ... (see manual)\n"
    "  quit\n";
  cout << endl;
  exit(0);
}


/**********************************************************/
/*     Warning Message                                    */
/**********************************************************/
void WarningMessage()
{
  cerr << "WARNING   : " << message.str()      << endl;
  message.str("");
}

void Notice(string module){
  if(module == "NOTE"){
    cerr << " (._.) " << message.str() << endl;
  }
  else{
    if(verbflag) cerr << " (@_@) [" << module << "] " << message.str() << endl;
  }
  message.str("");
}


/**********************************************************/
/*     Emergency Stop                                     */
/**********************************************************/
int TerminateCode(string module)
{
  DeceFreeMemory();
  cerr << "ERROR     :[" << module << "] " << message.str() << endl;
  exit(-1);
}


/**********************************************************/
/*     Banner                                             */
/**********************************************************/
void DeceBanner()
{
  for(int i=55 - strlen(&version[0]) ; i>0 ; i--) cout << " ";
  cout << version << "\n";
  cout
    <<"     oooooooooo                   oooooo    ooooooooooo      Am Brunnen vor dem Thore\n"        
    <<"      888     Y8b               d8P    88b   888      8      Da steht ein Lindenbaum\n"         
    <<"      888      888    ooooo    888           888             Ich traeauumt in seinem Schatten\n" 
    <<"      888      888  o88   888  888           888oooo8        So manchen suessen Traum.\n"       
    <<"      888      888  888ooo888  888           888             \n"
    <<"      888     d88   Y88     ,   88b     oo   888      8      Ich schnitt in seine Rinde\n"      
    <<"     oooobood8P      Y8bod8P     Y8bood8P   ooooboooooo      So manches liebe Wort\n"           
    <<"                                                             Es zog in Freud und Leide\n"       
    <<"                                                             Zu ihm mich immer fort.\n"         
    <<"\n"
    <<"                                            Der Lindenbaum, Wilhelm Mueller / Franz Schubert\n";
}


/**********************************************************/
/*      ENDF Help                                         */
/**********************************************************/
void DeceHelpENDF()
{
  cout <<
"MF Numbers\n"
"   1: General Information\n"
"   2:Resonance Parameters\n"
"   3:Reaction Cross Sections\n"
"   4:Angular Distributions of Secondary Particles\n"
"   5:Energy Distributions of Secondary Particles\n"
"   6:Product Energy - Angle Distributions\n"
"   7:Thermal Neutron Scattering Law Data\n"
"   8:Radioactive Decay and Fission Product Data\n"
"   9:Multiplicities for Production of Radioactive Nuclides\n"
"  10:Cross Sections for Production of Radioactive Nuclides \n"
"\n"
"  12:Photon Production Multiplicities and Transition probability arrays\n"
"  13:Photon Production Cross Sections\n"
"  14:Photon Angular Distributions\n"
"  15:Continuous Photon Energy Spectra\n"
"\n"
"  30:Data Covariances Obtained from Parameter Covariances and Sensitivities\n"
"  31:Covariances of the Average Number of Neutrons per Fission\n"
"  32:Covariances of Resonance Parameters\n"
"  33:Covariances of Neutron Cross Sections \n"
"  34:Covariances for Angular Distributions of Secondary Particles\n"
"  35:Covariances for Energy Distributions of Secondary Particles\n"
"  40:Covariances for Production of Radioactive Nuclei\n"
"\n"
"\n"
"MT numbers\n"
"\n"
"  1: Neutron total cross sections\n"
"  2: Elastic scattering cross section for incident particles\n"
"  3: Nonelastic neutron cross section\n"
"  4: Inelastic scattering, Sum of MT = 50 - 91\n"
"  5: Sum of all reactions not given explicitly in another MT number\n"
"\n"
" 11: (z, 2n d)\n"
" 16: (z, 2n)\n"
" 17: (z, 3n)\n"
" 18: (z, fission)\n"
" 19: (n, f) First-chance neutron-induced fission\n"
" 20: (n, n f) Second-chance neutron-induced fission\n"
" 21: (n, 2n f) Third-chance neutron-induced fission\n"
" 22: (z, n alpha)\n"
" 23: (n, n3 alpha) \n"
" 24: (z, 2n alpha)\n"
" 25: (z, 3n alpha)\n"
" 28: (z, n p) \n"
" 29: (z, n 2alpha)\n"
" 30: (z, 2n 2alpha)\n"
" 32: (z, n d)\n"
" 33: (z, n t)\n"
" 34: (z, n 3He)\n"
" 35: (z, n d 2alpha)\n"
" 36: (z, n t 2alpha)\n"
" 37: (z, 4n)\n"
" 38: (n, 3n f) Fourth-chance fission cross section\n"
" 41: (z, 2n p)\n"
" 42: (z, 3n p)\n"
" 44: (z, n 2p) \n"
" 45: (z, n p alpha)\n"
" 50: (y, n0) Production of a neutron, leaving the residual nucleus in the ground state\n"
" 51 - 90: (z, nx) Production of a neutron, with residual in the x-th excited state\n"
" 91: (z, nc) Production of a neutron in the continuum\n"
"\n"
"102: (z, gamma) Radiative capture\n"
"103: (z, p)\n"
"104: (z, d)\n"
"105: (z, t)\n"
"106: (z, 3He)\n"
"107: (z, alpha)\n"
"108: (z, 2alpha)\n"
"109: (z, 3alpha)\n"
"111: (z,2p)\n"
"112: (z, p alpha)\n"
"113: (z, t 2alpha)\n"
"114: (z, d 2alpha)\n"
"115: (z, p d)\n"
"116: (z, p t)\n"
"117: (z, d alpha)\n"
"\n"
"151: Resonance parameters\n"
"\n"
"451: Heading or title information\n"
"452: nu_T , average total (prompt plus delayed) number of neutrons released per fission\n"
"454: Independent fission product yield data\n"
"455: nu_d, average number of delayed neutrons released per fission\n"
"456: nu_p, average number of prompt neutrons released per fission\n"
"457: Radioactive decay data\n"
"458: Energy release in fission for incident neutrons\n"
"459: Cumulative fission product yield data\n"
"460: Delayed fission photons\n"
"\n"
"600 - 648: (z, p_x) Production of a proton, with residual in the x-th excited state\n"
"649: (z, p_c) Production of a proton in the continuum\n"
"650 - 698: (z, d_x) Production of a deuteron, with residual in the x-th excited state\n"
"699: (z, d_c) Production of a deuteron in the continuum\n"
"700 - 748: (z, t_x) Production of a triton, with residual in the x-th excited state\n"
"749: (z, t_c) Production of a triton in the continuum\n"
"750 - 798: (z, 3He_x) Production of a 3He, with residual in the x-th excited state\n"
"799: (z, 3He_c) Production of a 3He in the continuum\n"
"800 - 848: (z, alpha_x) Production of an alpha partcle, with residual in the x-th excited state\n"
"849: (z, alpha_c) Production of an alpha in the continuum\n"
"\n";

  exit(0);
}
