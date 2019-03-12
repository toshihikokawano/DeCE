/******************************************************************************/
/**     Command Analyzer                                                     **/
/******************************************************************************/
#include <cstring>
#include <iostream>
#include <cstdlib>

using namespace std;

#include "command.h"

static double getval  (string);
static void   gettext (string, char *);
static int argc = 0;

CLine cmd;
string operation = "";

/**********************************************************/
/*      Skip comment-line and read 1-line                 */
/**********************************************************/
int CmdFgetOneline(void)
{
  /*** read one line, skip comment or blank line */
  while(1){
    cin.getline(cmd.line,MAX_TEXTLENGTH-1);
    if(cin.eof() !=0 ) return(-1);
    if(cmd.line[0] == '#') continue;
    if(strlen(cmd.line) <= 1) continue;
    break;
  }
  return(0);
}


/**********************************************************/
/*      Extract Parameters from Line                      */
/**********************************************************/
string CmdExtractArgument(void)
{
  char   work[MAX_TEXTLENGTH];
  string arg;
  string ope;
  char   *tok = NULL;
  string d1 = " ";
  string d2 = " +-*/=:";

  strcpy(work,cmd.line);
  tok = strtok(work,d2.c_str());

  if(tok!=NULL) ope = (string)tok;
  else          ope = "";

  for(char *p = &ope[0] ; *p ; p++) *p = tolower(*p);

  argc = 1;
  if(ope == "calc"){
    cmd.mt      = (int)getval(d2);
    cmd.opt1    = (int)getval(d2);
    cmd.opt2    = (int)getval(d2);
    cmd.text[0] =  CmdCheckChar('+') * '+'
                 + CmdCheckChar('-') * '-'
                 + CmdCheckChar('*') * '*'
                 + CmdCheckChar('/') * '/'
                 + CmdCheckChar(':') * ':';
  }
  else if(ope == "tpid"){
    CmdExtractString(cmd.text);
  }
  else if(ope == "read" || ope == "angdist" || ope == "libread" || ope == "mergeread"){
    cmd.mf    = (int)getval(d1);
    cmd.mt    = (int)getval(d1);
    CmdExtractString(cmd.text);
    cmd.opt1  = (int)getval(d1);
    cmd.mtend = cmd.mt;
  }
  else if(ope == "multiread" || ope == "multiangdist" || ope == "multilibread"){
    cmd.mf    = (int)getval(d1);
    cmd.mt    = (int)getval(d1);
    cmd.mtend = (int)getval(d1);
    CmdExtractString(cmd.text);
    cmd.opt1  = (int)getval(d1);
  }
  else if(ope == "addpoint" || ope == "delpoint"){
    cmd.mf    = (int)getval(d1);
    cmd.mt    = (int)getval(d1);
    cmd.x     =      getval(d1);
    cmd.y     =      getval(d1);
  }
  else if(ope == "factor"){
    cmd.mf    = (int)getval(d1);
    cmd.mt    = (int)getval(d1);
    cmd.x     =      0.0;
    cmd.y     =      getval(d1);
    cmd.xmin  =      getval(d1);
    cmd.xmax  =      getval(d1);
  }
  else if(ope == "normalize"){
    cmd.mf    = (int)getval(d1);
    cmd.mt    = (int)getval(d1);
    cmd.x     =      getval(d1);
    cmd.y     =      getval(d1);
    cmd.xmin  =      getval(d1);
    cmd.xmax  =      getval(d1);
  }
  else if(ope == "applyfunc1" || ope == "applyfunc2" || ope == "applyfunc3"){
    cmd.mf    = (int)getval(d1);
    cmd.mt    = (int)getval(d1);
    cmd.x     =      getval(d1);
    cmd.y     =      getval(d1);
    cmd.xmin  =      getval(d1);
    if(     ope == "applyfunc1") cmd.opt1 = 1;
    else if(ope == "applyfunc2") cmd.opt1 = 2;
    else if(ope == "applyfunc3") cmd.opt1 = 3;
  }
  else if(ope == "readjust"){
    cmd.mt    = (int)getval(d1);
  }
  else if(ope == "multidelete"){
    cmd.mf    = (int)getval(d1);
    cmd.mt    = (int)getval(d1);
    cmd.mtend = (int)getval(d1);
    cmd.opt1  = (int)getval(d1);
    cmd.opt2  = (int)getval(d1);
  }
  else if(ope == "duplicate" || ope == "copy"){
    cmd.mf    = (int)getval(d1);
    cmd.mt    = (int)getval(d1);
    cmd.opt1  = (int)getval(d1);
  }
  else if(ope == "changeint"){
    cmd.mt    = (int)getval(d1);
    cmd.opt1  = (int)getval(d1);
    cmd.opt2  = (int)getval(d1);
    cmd.opt3  = (int)getval(d1);
  }
  else if(ope == "changeqval"){
    cmd.mt    = (int)getval(d1);
    cmd.x     =      getval(d1);
    cmd.y     =      getval(d1);
  }
  else if(ope == "boundcorrect"){
    cmd.mt    = (int)getval(d1);
  }
  else if(ope == "duplicatepoint"){
    cmd.mt    = (int)getval(d1);
    cmd.x     =      getval(d1);
  }
  else if(ope == "genprod"){
    cmd.mt    = (int)getval(d1);
    cmd.opt1  = (int)getval(d1);
  }
  else if(ope == "reconstruct" || ope == "reconangdist" || ope == "smoothangdist"){
    cmd.xmin  =      getval(d1);
    cmd.xmax  =      getval(d1);
    cmd.x     =      getval(d1);
  }
  else if(ope == "set" || ope == "unset"){
    gettext(d1,cmd.text);
  }
  else{
    cmd.mf    = (int)getval(d1);
    cmd.mt    = (int)getval(d1);
    cmd.mtend = cmd.mt;
    cmd.opt1  = (int)getval(d1);
    cmd.opt2  = (int)getval(d1);
  }

  operation = ope;
  return(ope);
}


double getval(string delim)
{
  char   work[MAX_TEXTLENGTH];
  double x = 0.0;
  char   *tok = NULL;

  strcpy(work,cmd.line);
  tok = strtok(work,delim.c_str());

  for(int i=0 ; i<argc ; i++) tok = strtok(NULL,delim.c_str());

  if(tok != NULL) x = atof(tok);

  argc ++;

  return(x);
}


void gettext(string delim, char *str)
{
  char   work[MAX_TEXTLENGTH];
  char   *tok = NULL;

  strcpy(work,cmd.line);
  tok = strtok(work,delim.c_str());

  for(int i=0 ; i<argc ; i++) tok = strtok(NULL,delim.c_str());

  if(tok != NULL) strcpy(str,tok);

  argc ++;
}


/**********************************************************/
/*      Return Command                                    */
/**********************************************************/
string CmdGetOperation()
{
  return operation;
}


/**********************************************************/
/*      Extract Comma Quoted Text from Line               */
/**********************************************************/
void CmdExtractString(char *d)
{
  char *s = cmd.line;
  int lens = strlen(s);

  argc++;

  int i0=0, i1=0;
  for(int i=0 ; i<lens ; i++){
    if(s[i] == '"'){
      i0=i+1;
      break;
    }
  }
  if(i0 == 0) return;

  for(int i=i0 ; i<lens ; i++){
    if(s[i] == '"'){
      i1=i-1;
      break;
    }
  }

  if(i1-i0+1 < MAX_TEXTLENGTH){
    for(int i=i0 ; i<=i1 ; i++)  d[i-i0] = s[i];
    d[i1-i0+1] = '\0';
  }
}


/**********************************************************/
/*      Scan Characte in Line                             */
/**********************************************************/
bool CmdCheckChar(char c)
{
  if(strchr(cmd.line,c)) return true;
  else                   return false;
}
