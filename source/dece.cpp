/******************************************************************************/
/**                                                                          **/
/**     DeCE : The Descriptive Correction of ENDF-6 Format Code              **/
/**                                          1.2.5 (Jadeite)  Nov.  2020     **/
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

static string version  = "1.2.5 Jadeite (Nov. 2020)";
static bool   verbflag = false;
static bool   justquit = false;
static bool   filescan = false;
static int    newsec   = 0;
static ENDF   *lib[MAX_SECTION];

static void DeceMain          (string, string, ENDFDict *);
static void DeceStoreData     (ENDFDict *, ifstream *);
static void DeceReadMonitor   (const int, const int, const int, const int, const int, const int, const int);
static void DeceHelp          (void);
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
  while((p=getopt(argc,argv,"o:f:t:e:rqsvh"))!=-1){
    switch(p){
    case 'o':  libname_out = optarg;   break;
    case 'f':  mfin = atoi(optarg);    break;
    case 't':  mtin = atoi(optarg);    break;
    case 'e':  ein  = atof(optarg);    break;
    case 'r':  reconr   = true;        break;
    case 'q':  justquit = true;        break;
    case 's':  filescan = true;
               justquit = true;        break;
    case 'v':  verbflag = true;        break;
    case 'h':  DeceHelp();             break;
    default:                           break;
    }
  }

  /*** check if ENDF file is given */
  if(optind < argc) libname_in = argv[optind];
  if(libname_in == ""){ message << "ENDF-6 formattted file not given"; TerminateCode("main"); }
  if(libname_in == libname_out){ message << "same in/out file names"; TerminateCode("main"); }

  /*** tabular output from one section */
  if((mfin == 2) || (mfin == 32)) mtin = 151;

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
    if(filescan) DeceScanIndex(dict);
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


















