/******************************************************************************/
/**                                                                          **/
/**     DeCE : The Descriptive Correction of ENDF-6 Format Code              **/
/**                                          1.2.1 (Adularia)  May  2016     **/
/**                                                            T. Kawano     **/
/**                                       Los Alamos National Laboratory     **/
/******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

using namespace std;

#include "dece.h"
#include "command.h"
#include "terminate.h"

static string version  = "1.2.1 Adularia (May 2015)";
static bool   verbflag = false;
static bool   justquit = false;

static void DeceMain          (string, string, ENDFDict *);
static void DeceStoreData     (ENDFDict *, ifstream *);
static void DeceHelp          (void);
static void DeceCreateLib     (ENDFDict *, int, int);
static void DeceFreeMemory    (void);
static void DeceBanner        (void);

static inline void DeceReadMonitor   (int, int, int, int, int);
static inline void DeceCheckMT       (int);


/**********************************************************/
/*      Global Parameters                                 */
/**********************************************************/

string   tempfile = "DECETempfile.dat";
int      newsec = 0;
ENDF     *lib[MAX_SECTION];

/*** defined in command.cpp */
extern CLine cmd;


/**********************************************************/
/*      DeCE Main                                         */
/**********************************************************/
int main(int argc, char *argv[])
{
  ENDFDict dict;
  int      p,mfin=0,mtin=0;
  double   ein=0.0;
  string   libname_in = "",  libname_out = "";
  bool     reconr = false;

  /*** command line options */
  while((p=getopt(argc,argv,"o:f:t:e:rqvh"))!=-1){
    switch(p){
    case 'o':  libname_out = optarg;   break;
    case 'f':  mfin = atoi(optarg);    break;
    case 't':  mtin = atoi(optarg);    break;
    case 'e':  ein  = atof(optarg);    break;
    case 'r':  reconr   = true;        break;
    case 'q':  justquit = true;        break;
    case 'v':  verbflag = true;        break;
    case 'h':  DeceHelp();             break;
    default:                           break;
    }
  }

  /*** check if ENDF file is given */
  if(optind < argc) libname_in = argv[optind];
  if(libname_in == "") TerminateCode("ENDF-6 formattted file not given");
  if(libname_in == libname_out) TerminateCode("same in/out file names");
  
  /*** tabular output from one section */
  if(mfin == 2) mtin = 151;

  if((mfin > 0) && (mtin > 0)){
    ifstream  fpin;
    fpin.open(&libname_in[0]);  if(!fpin) TerminateCode("ENDF file cannot open",libname_in);
    DeceCheckMT(mtin);
    if(ein > 0.0) DeceDataPoint(&fpin,mfin,mtin,ein);
    else          DeceFileToTable(&fpin,mfin,mtin,0);
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
      gfrScanThermal(&fpin,&dict);
      fpin.close();
    }

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

  /*** for all operation commands */
  do{
    if( justquit ) break;
    if( CmdFgetOneline()<0 ) break;
    ope = CmdExtractArgument();

    /*** END: terminate code */
    if(ope == "end" || ope == "quit" || ope == "exit") break;

    /*** CALC: manipulate TAB1 record */
    else if(ope == "calc"){
      DeceCreateLib(dict,3,cmd.mt);
      DeceCalc(dict,lib,cmd.mt,cmd.opt1,cmd.opt2,cmd.text[0]);
    }

    /*** MAKE4: generate / reconstruct total inelastic (abbrev. of CALC) */
    else if(ope == "make4"){
      DeceCreateLib(dict,3,4);
      DeceDuplicate(dict,lib,3,51,4);
      DeceCalc(dict,lib,4,51,91,':');
    }

    /*** DUPLICATE: copy data */
    else if(ope == "duplicate"){
      DeceCreateLib(dict,cmd.mf,cmd.opt1);
      DeceDuplicate(dict,lib,cmd.mf,cmd.mt,cmd.opt1);
    }

    /*** READ: read tabulated cross section data file, and replace section */
    else if( (ope == "read") || (ope == "multiread") || (ope == "mergeread") ){
      for(int mt=cmd.mt ; mt <= cmd.mtend ; mt++){
        DeceCreateLib(dict,cmd.mf,mt);
        if(ope == "mergeread")
          DeceRead(dict,lib[dict->getID(cmd.mf,mt)],cmd.mf,mt,cmd.text,cmd.opt1,true);
        else
          DeceRead(dict,lib[dict->getID(cmd.mf,mt)],cmd.mf,mt,cmd.text,cmd.opt1,false);
      }
    }

    /*** ANGDIST: read tabulated angular distribution data file */
    else if( (ope == "angdist") || (ope == "multiangdist") ){
      for(int mt=cmd.mt ; mt <= cmd.mtend ; mt++){
        DeceCreateLib(dict,cmd.mf,mt);
        DeceAngdist(dict,lib,cmd.mf,mt,cmd.text,cmd.opt1);
      }
    }

    /*** LIBREAD: read a section from another ENDF file, and replace */
    else if( (ope == "libread") || (ope == "multilibread") ){
      for(int mt=cmd.mt ; mt <= cmd.mtend ; mt++){
        DeceCreateLib(dict,cmd.mf,mt);
        DeceLibRead(dict,lib[dict->getID(cmd.mf,mt)],cmd.text);
      }
    }

    /*** TABLE: tabulate MF3, MF4, MF6 data */
    else if(ope == "table"){
      DeceTable(dict,lib,&fpin,cmd.mf,cmd.mt,cmd.opt1);
    }

    /*** EXTRACT: dead copy section */
    else if(ope == "extract" ){
      DeceExtract(dict,lib,&fpin,cmd.mf,cmd.mt);
    }

    /*** ADDPOINT, DELPOINT: add / remove a point in TAB1 record */
    else if( (ope == "addpoint") || (ope == "delpoint") ){
      DeceCheckMT(cmd.mt);
      DecePoint(dict,lib,cmd.mf,cmd.mt,cmd.x,cmd.y,ope);
    }

    /*** FACTOR, NORMALIZE: multiply by a factor */
    else if(ope == "factor" || ope == "normalize"){
      DeceCheckMT(cmd.mt);
      DeceFactor(dict,lib,cmd.mf,cmd.mt,cmd.x,cmd.y,cmd.xmin,cmd.xmax);
    }

    /*** FUNC1, FUNC2: multiply by a factor that is calculated with a function */
    else if(ope == "applyfunc1" || ope == "applyfunc2"){
      DeceCheckMT(cmd.mt);
      DeceApplyFunc(dict,lib,cmd.mf,cmd.mt,cmd.opt1,cmd.x,cmd.y,cmd.xmin);
    }

    /*** DELETE: delete section */
    else if( (ope == "delete") || (ope == "multidelete") ){
      for(int mt=cmd.mt ; mt <= cmd.mtend ; mt++){
        DeceDelete(dict,cmd.mf,mt);
      }
    }

    /*** CHANGEINT: change interpolation scheme in MF3 */
    else if(ope == "changeint"){
      DeceCheckMT(cmd.mt);
      DeceChangeInt(dict,lib,cmd.mt,cmd.opt1,cmd.opt2,cmd.opt3);
    }

    /*** CHANGEQVAL: change Q-vales in MF3 */
    else if(ope == "changeqval"){
      DeceCheckMT(cmd.mt);
      DeceChangeQvalue(dict,lib,cmd.mt,cmd.x,cmd.y);
    }

    /*** FIXAWR: fix AWR in at the top of the file */
    else if(ope == "fixawr"){
      DeceFixAWR(dict);
    }

    /*** miscellaneous MF1 manipulations */
    /*** NUTOTAL: generate / reconstruct total nu-bar as the sum of 455 and 456 */
    else if(ope == "nutotal"){
      DeceCreateLib(dict,1,452);
      DeceCalc452(dict,lib);
    }

    /*** miscellaneous MF6 manipulations */
    /*** BOUNDCORRECT: correct energy boundary in continuum in MF6 */
    else if(ope == "boundcorrect"){
      DeceCheckMT(cmd.mt);
      DeceBoundCorrect(dict,lib,cmd.mt);
    }

    /*** DUPLICATEPOINT: duplicate the last point in MF6 */
    else if(ope == "duplicatepoint"){
      DeceCheckMT(cmd.mt);
      DeceDuplicatePoint(dict,lib,cmd.mt,cmd.x);
    }

    /*** GENPROD: generate production cross section from MF6 MT5 */
    else if(ope == "genprod"){
      DeceCreateLib(dict,3,cmd.mt);
      DeceGenProdCS(dict,lib,cmd.mt,cmd.opt1);
    }

    /*** resonance manipulation */
    /*** RECONSTRUCT: reconstruct cross sections from resonances */
    else if(ope == "reconstruct"){
      if(dict->getID(2,151) >= 0){
        gfrPtCross(dict,lib,cmd.xmin,cmd.xmax,cmd.x);
      }
    }

    /*** RECONANGDIST: calculate Legendre coefficients from resonance parameters */
    else if(ope == "reconangdist"){
      if(dict->getID(2,151) >= 0){
        gfrAngDist(dict,lib,cmd.xmin,cmd.xmax,cmd.x);
      }
    }

    /*** SMOOTHANGDIST: calculate energy averaged Legendre coefficients in RRR */
    else if(ope == "smoothangdist"){
      if(dict->getID(2,151) >= 0){
        gfrAngDistSmooth(dict,lib,cmd.xmin);
      }
    }

    /*** TPID: replace TPID */
    else if(ope == "tpid")     CmdExtractString(dict->tpid);

    /*** INDEX: print section index */
    else if(ope == "index")    DeceScan(dict);

    /*** data processing */
    /*** POINTWISE: create pointwise cross section */
    else if(ope == "pointwise"){
      if(dict->getID(2,151) >= 0){
        DeceCreateLib(dict,3,901);
        DeceCreateLib(dict,3,902);
        DeceCreateLib(dict,3,903);
        DeceCreateLib(dict,3,904);
        DeceGeneratePointwise(dict,lib);
      }
    }
    
    /*** Unknown command */
    else TerminateCode("command not found",ope);

  }while(!cin.eof());


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

  for(int i=0 ; i<dict->sec ; i++){

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
//  else if(dict->mf[i] == 32) size = L; // MF32 temporarily disabled
    else if(dict->mf[i] == 33) size = L;
    else continue;

    lib[newsec] = new ENDF(size);
    ENDFRead(fp,lib[newsec],dict->mf[i],dict->mt[i]);
    if(dict->mf[i] == 2) ENDFMF2boundary(dict,lib[newsec]);

    if(verbflag) DeceReadMonitor(lib[newsec]->getENDFmat(),dict->mf[i],dict->mt[i],newsec,lib[newsec]->getPOS());
    dict->setID(i,newsec++);

    if(newsec >= MAX_SECTION) TerminateCode("Too many sections",newsec);
  }
}


/**********************************************************/
/*      Data Read Monitor                                 */
/**********************************************************/
void DeceReadMonitor(int mat, int mf, int mt, int sec, int n)
{
  cerr << " (@_@) <";
  cerr << " MAT:" << mat;
  cerr << " MF:" << mf;
  cerr << " MT:" << mt;
  cerr << " assigned for Section " << sec;
  cerr << " sub-blocks " << n << endl;
}


/**********************************************************/
/*      Check If Valid MT                                 */
/**********************************************************/
inline void DeceCheckMT(int mt)
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
  if( dict->getID(mf,mt)>=0 ) return;

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

  lib[newsec]->setENDFmat(dict->mat);
  lib[newsec]->setENDFmf(mf);
  lib[newsec]->setENDFmt(mt);

  int i = dict->scanDict(mf,mt);
  if(i<0){
    dict->addDict(mf,mt,0,newsec);
  }
  else{
    dict->setID(i,newsec);
  }

  newsec++;
  if(newsec >= MAX_SECTION) TerminateCode("Too many sections");
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
void WarningMessage(string msg)
{
  if(verbflag) cerr << "WARNING   :" << msg << endl;
}

void WarningMessage(string msg, int n)
{
  if(verbflag) cerr << "WARNING   :" << msg << n << endl;
}

void WarningMessage(string msg, double x)
{
  if(verbflag) cerr << "WARNING   :" << msg << x << endl;
}

void WarningMessage(string msg, string x)
{
  if(verbflag) cerr << "WARNING   :" << msg << x << endl;
}


/**********************************************************/
/*     Emergency Stop                                     */
/**********************************************************/
int TerminateCode(string msg)
{
  DeceFreeMemory();
  cerr << "ERROR     :" << msg << endl;
  exit(-1);
}

int TerminateCode(string msg, int n)
{
  DeceFreeMemory();
  cerr << "ERROR     :" << msg << " : " << n << endl;
  exit(-1);
}

int TerminateCode(string msg, double x)
{
  DeceFreeMemory();
  cerr << "ERROR     :" << msg << " : " << x << endl;
  exit(-1);
}

int TerminateCode(string msg, string x)
{
  DeceFreeMemory();
  cerr << "ERROR     :" << msg << " : " << x << endl;
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

