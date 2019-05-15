/******************************************************************************/
/**                                                                          **/
/**     DeCE : The Descriptive Correction of ENDF-6 Format Code              **/
/**                                          1.2.2 (Pyrite)  March  2019     **/
/**                                                            T. Kawano     **/
/**                                       Los Alamos National Laboratory     **/
/******************************************************************************/

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

using namespace std;

#include "dece.h"
#include "command.h"
#include "global.h"
#include "terminate.h"

static string version  = "1.2.2 Pyrite (March 2019)";
static bool   verbflag = false;
static bool   justquit = false;
static bool   renumber = false;
static bool   filescan = false;
static int    newsec   = 0;
static ENDF   *lib[MAX_SECTION];

static void DeceMain          (string, string, ENDFDict *);
static void DeceStoreData     (ENDFDict *, ifstream *);
static void DeceInitOptions   (void);
static void DeceReadMonitor   (int, int, int, int, int);
static void DeceHelp          (void);
static void DeceFreeMemory    (void);
static void DeceBanner        (void);


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
  while((p=getopt(argc,argv,"o:f:t:e:rqnsvh"))!=-1){
    switch(p){
    case 'o':  libname_out = optarg;   break;
    case 'f':  mfin = atoi(optarg);    break;
    case 't':  mtin = atoi(optarg);    break;
    case 'e':  ein  = atof(optarg);    break;
    case 'r':  reconr   = true;        break;
    case 'q':  justquit = true;        break;
    case 'n':  renumber = true;        break;
    case 's':  filescan = true;
               justquit = true;        break;
    case 'v':  verbflag = true;        break;
    case 'h':  DeceHelp();             break;
    default:                           break;
    }
  }

  /*** check if ENDF file is given */
  if(optind < argc) libname_in = argv[optind];
  if(libname_in == "") TerminateCode("ENDF-6 formattted file not given");
  if(libname_in == libname_out) TerminateCode("same in/out file names");

  /*** initialize global options */
  DeceInitOptions();

  /*** tabular output from one section */
  if(mfin == 2) mtin = 151;

  if((mfin > 0) && (mtin > 0)){
    ifstream  fpin;
    fpin.open(&libname_in[0]);  if(!fpin) TerminateCode("ENDF file cannot open",libname_in);
    DeceCheckMT(mtin);
    if(ein > 0.0) DeceDataPoint(&fpin,mfin,mtin,ein);
    else          DeceFileToTable(&fpin,mfin,mtin);
    fpin.close();
  }

  else if( (mfin == 0) && (mtin == 0) ){
    /*** scan the file, and store all the MF,MT sections */
    int ctl = ENDFScanLibrary(libname_in,&dict);
    if     (ctl == -1) TerminateCode("ENDF file cannot open",libname_in);
    else if(ctl == -2) TerminateCode("too many sections");

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
  else TerminateCode("MF or MT number not given");

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
  DataSize size = M;

  for(int i=0 ; i<dict->getSEC() ; i++){

    if(dict->mf[i] ==  1){
      if( (dict->mt[i] == 452) || (dict->mt[i] == 455) || (dict->mt[i] == 456) ) size = M;
      else if(dict->mt[i] == 460) size = L;
      else continue;
    }
    else if(dict->mf[i] ==  2) size = L;
    else if(dict->mf[i] ==  3) size = M;
    else if(dict->mf[i] ==  4) size = L;
    else if(dict->mf[i] ==  5) size = L;
    else if(dict->mf[i] ==  6) size = L;
    else if(dict->mf[i] ==  8) size = M;
    else if(dict->mf[i] ==  9) size = M;
    else if(dict->mf[i] == 10) size = M;
    else if(dict->mf[i] == 12) size = (dict->mt[i] == 460) ? L : M;
    else if(dict->mf[i] == 13) size = M;
    else if(dict->mf[i] == 14) size = M;
    else if(dict->mf[i] == 15) size = M;
    else if(dict->mf[i] == 31) size = M;
    else if(dict->mf[i] == 32) size = L;
    else if(dict->mf[i] == 33) size = L;
    else if(dict->mf[i] == 34) size = M;
    else if(dict->mf[i] == 35) size = L;
    else continue;

    lib[newsec] = new ENDF(size);
    ENDFRead(fp,lib[newsec],dict->mf[i],dict->mt[i]);
    if(dict->mf[i] == 2) ENDFMF2boundary(dict,lib[newsec]);

    DeceReadMonitor(lib[newsec]->getENDFmat(),dict->mf[i],dict->mt[i],newsec,lib[newsec]->getPOS());
    dict->setID(i,newsec++);

    if(newsec >= MAX_SECTION) TerminateCode("Too many sections",newsec);
  }
}


/**********************************************************/
/*      Initialize Global Options                         */
/**********************************************************/
void DeceInitOptions()
{
  ENDFPrintLineNumber(renumber);
}


/**********************************************************/
/*      Data Read Monitor                                 */
/**********************************************************/
void DeceReadMonitor(int mat, int mf, int mt, int sec, int n)
{
  message << "MAT:" << mat <<" ";
  message << "MF:" << mf;
  message << "MT:" << mt;
  message << " assigned for Section " << sec;
  message << " sub-blocks " << n;
  Notice("DeceReadMonitor");
}


/**********************************************************/
/*      Check If Valid MT                                 */
/**********************************************************/
void DeceCheckMT(int mt)
{
  if( mt <= 0 || 1000 <= mt ) TerminateCode("invalid MT number",mt);
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
    if(mf >= 900){
      lib[newsec] = new ENDF(L); // above 900 used for temporal MTs
    }
    else{
      if( (mf == 2) || (mf == 4) || (mf == 6) )
        lib[newsec] = new ENDF(L);
      else if( (mf == 12) || (mf == 14) )
        lib[newsec] = new ENDF(S);
      else
        lib[newsec] = new ENDF(M);
    }
  }
  catch(bad_alloc){
    TerminateCode("memory allocation error");
    return;
  }

  lib[newsec]->setENDFmat(dict->getMAT());
  lib[newsec]->setENDFmf(mf);
  lib[newsec]->setENDFmt(mt);

  int i = dict->scanDict(mf,mt);
  if(i < 0){
    if(dict->addDict(mf,mt,0,newsec)) TerminateCode("Too many sections");
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
int TerminateCode(string msg)
{
  DeceFreeMemory();
  cerr << "ERROR     : " << msg << endl;
  exit(-1);
}

int TerminateCode(string msg, int n)
{
  DeceFreeMemory();
  cerr << "ERROR     : " << msg << " : " << n << endl;
  exit(-1);
}

int TerminateCode(string msg, double x)
{
  DeceFreeMemory();
  cerr << "ERROR     : " << msg << " : " << x << endl;
  exit(-1);
}

int TerminateCode(string msg, string x)
{
  DeceFreeMemory();
  cerr << "ERROR     : " << msg << " : " << x << endl;
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
    <<"     oooooooooo                   oooooo    ooooooooooo       Ich schlief, ich schlief-\n"
    <<"      888     Y8b               d8P    88b   888      8       Aus tiefem Traum bin ich erwacht:-\n"
    <<"      888      888    ooooo    888           888              Die Welt ist tief,\n"
    <<"      888      888  o88   888  888           888oooo8         Und tiefer als der Tag gedacht.\n"
    <<"      888      888  888ooo888  888           888              Tief ist ihr Weh-\n"
    <<"      888     d88   Y88     ,   88b     oo   888      8       Lust-tiefer noch als Herzeleid:\n"
    <<"     oooobood8P      Y8bod8P     Y8bood8P   ooooboooooo       Weh spricht: Vergeh!\n"
    <<"                                                              Doch alle Lust will Ewigkeit-\n"
    <<"                                                              - will tiefe, tiefe Ewigkeit!\n"
    <<"             O Mensch! Gib acht\n"
    <<"             Was spricht die tiefe Mitternacht?\n"
    <<"                                                   Also sprach Zarathustras, Friedrich Nietzsche\n";
}

