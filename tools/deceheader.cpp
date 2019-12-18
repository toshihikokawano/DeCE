/******************************************************************************/
/**                                                                          **/
/**     DeCE Tools : Print / Replace ENDF Header Data                        **/
/**                                                                          **/
/******************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <iomanip>
#include <unistd.h>

using namespace std;

#include "../source/endflib.h"
#include "deceheader.h"


int main(int, char *[]);
int main(int argc, char *argv[])
{
  ENDFDict dict;
  ifstream libin, datin;
  string   libname = "", datname = "";
  char     *cbuf = NULL, **line = NULL;

  int p = 0;
  while((p = getopt(argc,argv,"d:")) != -1){
    switch(p){
    case 'd': // data file when patched
      datname = optarg;
      break;
    default:
      break;
    }
  }
  if(optind < argc) libname = argv[optind];

  if(libname == ""){
    cerr << "ENDF file not given" << endl;
    exit(-1);
  }

  /*** scan the ENDF file and initiate DICT object */
  if(ENDFScanLibrary(libname,&dict) < 0){
    cerr << "ENDF file cannot open " << libname << endl;
    exit(-1);
  }


  /*** check MF1 to see how many text lines in the comment section */
  libin.open(libname.c_str());
  int nline = DeceHeaderScanMF1(&libin,line,false);

  /*** allocate 1-dim buffer to store all the text data */
  cbuf = new char [(MAX_FIELD + 1) * nline];
  line = new char * [nline];
  for(int i=0 ; i<nline ; i++) line[i] = &cbuf[(MAX_FIELD +1) * i];

  /*** this time, copy the comment section into the line array */
  DeceHeaderScanMF1(&libin,line,true);

  /*** check if the file has the standard text header */
  if( ! DeceHeaderCheckStandard(nline,line) ){
    cerr << "Provided ENDF file is not in a standard header format" << endl;
    exit(-1);
  }


  /*** copy all the text data into the TextField objects */
  DeceHeaderCopyData(nline,line);

  /*** main process*/ 
  /*** when data file is given, replace text by those in the data file */
  if(datname != ""){
    datin.open(datname.c_str());
    if(!datin){
      cerr << "Data file cannot open " << datname << endl;
      exit(-1);
    }
    DeceHeaderReadData(&datin);
    datin.close();

    /*** replace text data in Dictionary by inputs */
    DeceHeaderReplaceData(&dict);

    /*** print entire data file */
    DeceHeaderCreateLib(&libin,&dict,nline,line);
  }
  /*** if not given, print the data */
  else DeceHeaderPrint();


  libin.close();

  delete [] line;
  delete [] cbuf;

  return(0);
}


/*******************************************************************************
 Scan MF1/MT451 and return the total number of lines.
 When boolean copy is true, copy all the text fields into line array. */
int DeceHeaderScanMF1(ifstream *fpin, char *line[], bool copy)
{
  ENDF lib(S);
  string data, s;

  /*** skip first 4 lines */
  fpin->seekg(0,ios_base::beg);
  ENDFSeekHead(fpin,&lib,1,451);
  for(int i=0 ; i < 3 ; i++) getline(*fpin,data);


  int m = 0; // number of lines
  while( getline(*fpin,data) ){

    int n = data.length();
    if(n >= 75){
      s = data.substr(72,3);
      int mt = atoi(s.c_str());
      if(mt == 0) break;
    }

    /*** when copy flag is set, copy the text contents into the line array */
    if(copy){
      for(int i=0 ; i<66 ; i++) line[m][i] = (i < n) ? data[i] : ' ';
      line[m][66] = '\0';
    }

    m ++;
  }

  return m;
}


/*******************************************************************************
  Test the file provided to see if this is in the standard format. */
bool DeceHeaderCheckStandard(const int nline, char *line[])
{
  if(nline < 5) return false;

  char dash1[5], dash2[6], dash3[7];

  /*** first 3 lines should be there */
  strncpy(dash1, &line[2][0],4); dash1[4] = '\0';
  strncpy(dash2, &line[3][0],5); dash2[5] = '\0';
  strncpy(dash3, &line[4][0],6); dash3[6] = '\0';

  /*** when these lines are given, maybe this is in the standard format */
  bool test1 = strncmp(dash1,"----",4);
  bool test2 = strncmp(dash2,"-----",5);
  bool test3 = strncmp(dash3,"------",6);
  bool test4 = (line[2][4] == '-');
  bool test5 = (line[3][5] == '-');
  bool test6 = (line[4][6] == '-');

  if( !test1 && !test2 && !test3 && !test4 && !test5 && !test6) return true;
  else return false;
}


/*******************************************************************************
  copy text filed data into objects */
void DeceHeaderCopyData(const int nline, char *line[])
{
  if(nline < 2) return;

  zsymam.copy( MAX_FIELD,line[ zsymam.ypos]);
  alab.copy(   MAX_FIELD,line[   alab.ypos]);
  edate.copy(  MAX_FIELD,line[  edate.ypos]);
  auth.copy(   MAX_FIELD,line[   auth.ypos]);
  refer.copy(  MAX_FIELD,line[  refer.ypos]);
  ddate.copy(  MAX_FIELD,line[  ddate.ypos]);
  rdate.copy(  MAX_FIELD,line[  rdate.ypos]);
  endate.copy( MAX_FIELD,line[ endate.ypos]);

  if(nline < 5) return;

  libname.copy(MAX_FIELD,line[libname.ypos]);
  sublib.copy( MAX_FIELD,line[ sublib.ypos]);
  format.copy( MAX_FIELD,line[ format.ypos]);
}


/*******************************************************************************
 Read text data to be used for substitution */
void DeceHeaderReadData(ifstream *fpin)
{
  string data;
  char s1[WORD_LENGTH+1], s2[MAX_FIELD+1];

  while( getline(*fpin,data) ){
    int n = data.length();

    if(data[0] == '#') continue;

    if(n < WORD_LENGTH) continue;
    else if(n > WORD_LENGTH + MAX_FIELD) n = WORD_LENGTH + MAX_FIELD;

    /*** copy key word to s1 */
    for(int i=0 ; i<WORD_LENGTH ; i++) s1[i] = data[i];
    s1[WORD_LENGTH] = '\0';

    int i0 = WORD_LENGTH, i1;
    /*** when data are not given */
    if(i0 == n) i1 = n;
    /*** when text is given, remove all extra spaces after the text data, last point is i1 */
    else{
      for(i1=n-1 ; i1>=-i0 ; i1--) if(data[i1] != ' ') break;
    }

    /*** copy the actual data into s2 */
    for(int i=i0 ; i<=i1 ; i++) s2[i-i0] = data[i];

    /*** fill right side of s2 by spaces */
    for(int i=i1-i0+1 ; i<MAX_FIELD ; i++) s2[i] = ' ';
    s2[MAX_FIELD] = '\0';

    if(     !strncmp(s1,"  ZSYMAM:",WORD_LENGTH)){
      zsymam.read(MAX_FIELD,s2);
    }
    else if(!strncmp(s1,"    ALAB:",WORD_LENGTH)){ alab.read(MAX_FIELD,s2); }
    else if(!strncmp(s1,"    AUTH:",WORD_LENGTH)){ auth.read(MAX_FIELD,s2); }
    else if(!strncmp(s1,"   REFER:",WORD_LENGTH)){ refer.read(MAX_FIELD,s2); }
    else if(!strncmp(s1,"   EDATE:",WORD_LENGTH)){ edate.read(MAX_FIELD,s2); }
    else if(!strncmp(s1,"   DDATE:",WORD_LENGTH)){ ddate.read(MAX_FIELD,s2); }
    else if(!strncmp(s1,"   RDATE:",WORD_LENGTH)){ rdate.read(MAX_FIELD,s2); }
    else if(!strncmp(s1,"  ENDATE:",WORD_LENGTH)){ endate.read(MAX_FIELD,s2); }
    else if(!strncmp(s1," LIBNAME:",WORD_LENGTH)){ libname.read(MAX_FIELD,s2); }
    else if(!strncmp(s1,"  SUBLIB:",WORD_LENGTH)){ sublib.read(MAX_FIELD,s2); }
    else if(!strncmp(s1,"  FORMAT:",WORD_LENGTH)){ format.read(MAX_FIELD,s2); }
    else{
      cerr << " Keyword " << s1 << "not defined" << endl;
    }
  }
}


/*******************************************************************************
 replace text data in DICT object by given data */
void DeceHeaderReplaceData(ENDFDict *dict)
{
  zsymam.paste( MAX_FIELD,dict->text[zsymam.ypos]);
  alab.paste(   MAX_FIELD,dict->text[  alab.ypos]);
  edate.paste(  MAX_FIELD,dict->text[ edate.ypos]);
  auth.paste(   MAX_FIELD,dict->text[  auth.ypos]);

  refer.paste(  MAX_FIELD,dict->text[ refer.ypos]);
  ddate.paste(  MAX_FIELD,dict->text[ ddate.ypos]);
  rdate.paste(  MAX_FIELD,dict->text[ rdate.ypos]);
  endate.paste( MAX_FIELD,dict->text[endate.ypos]);

  for(int i=44 ; i<55 ; i++) dict->text[1][i] = ' ';

  libname.paste(MAX_FIELD,dict->text[libname.ypos]);
  sublib.paste( MAX_FIELD,dict->text[ sublib.ypos]);
  format.paste( MAX_FIELD,dict->text[ format.ypos]);

  for(int i=0 ; i<4 ; i++) dict->text[2][i] = '-';
  for(int i=0 ; i<5 ; i++) dict->text[3][i] = '-';
  for(int i=0 ; i<6 ; i++) dict->text[4][i] = '-';

  sprintf(&dict->text[2][22],"MATERIAL "); dict->text[2][31] = ' ';
  sprintf(&dict->text[2][31],"%4d",dict->getMAT());
  for(int i=35 ; i<66 ; i++) dict->text[2][i] = ' ';

  for(int i=0 ; i<5 ; i++){
    dict->text[i][66] = '\0';
//  cout << dict->text[i] << "<" << endl;
  }
}


/*******************************************************************************
 Print entire file */
void DeceHeaderCreateLib(ifstream *fpin, ENDFDict *dict, const int nline, char *line[])
{
  string data;
  ENDF   lib(S);
  int    mf = 1, mt = 451;

  ENDFWriteTPID(dict);

  /*** header part */
  lib.setENDFmat(dict->getMAT());
  lib.setENDFmf(mf);
  lib.setENDFmt(mt);
  lib.setENDFhead(dict->getDICThead());
  ENDFWriteHEAD(&lib);

  for(int i=0 ; i<3 ; i++){
    Record r = dict->getDICTcont(i);
    ENDFWriteRecord(r);
    ENDFPrintRight(dict->getMAT(),mf,mt);
  }

  for(int i=0 ; i<nline ; i++){
    if(i < 5) cout << dict->text[i];
    else      cout << line[i];
    ENDFPrintRight(dict->getMAT(),mf,mt);
  }

  ENDFWriteSEND(&lib);
  ENDFWriteFEND(dict->getMAT());

  for(int mf=2 ; mf <= 40 ; mf++){
    bool printed = false;
    for(int i=0 ; i<dict->getSEC() ; i++){
      if(dict->mf[i] == mf){
        ENDFExtract(fpin,mf,dict->mt[i]);
        printed = true;
      }
    }
    if(printed) ENDFWriteFEND(dict->getMAT());
  }
  ENDFWriteFEND(0);
  ENDFWriteFEND(-1);
}


/*******************************************************************************
  show text data */
void DeceHeaderPrint()
{
  cout << "  ZSYMAM:"; zsymam.print(); cout << endl;
  cout << "    ALAB:"; alab.print();   cout << endl;
  cout << "    AUTH:"; auth.print();   cout << endl;
  cout << "   REFER:"; refer.print();  cout << endl;
  cout << "   EDATE:"; edate.print();  cout << endl;
  cout << "   DDATE:"; ddate.print();  cout << endl;
  cout << "   RDATE:"; rdate.print();  cout << endl;
  cout << "  ENDATE:"; endate.print(); cout << endl;

  cout << " LIBNAME:"; libname.print(); cout << endl;
  cout << "  SUBLIB:"; sublib.print();  cout << endl;
  cout << "  FORMAT:"; format.print();  cout << endl;
}
