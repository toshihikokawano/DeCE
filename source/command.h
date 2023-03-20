/******************************************************************************/
/**                                                                          **/
/**     Commands and Options                                                 **/
/**             extract given comand and its options                         **/
/**                                                                          **/
/******************************************************************************/

#define MAX_TEXTLENGTH 256

/**************************************/
/*      Class CLine : Command Line    */
/**************************************/
class CLine{
 public:
  char   line[MAX_TEXTLENGTH];
  char   parm[MAX_TEXTLENGTH];
  char   text[MAX_TEXTLENGTH];
  char   topt[MAX_TEXTLENGTH];
  int    mf;
  int    mt;
  int    mtend;
  int    opt1;
  int    opt2;
  int    opt3;
  double x;
  double y;
  double xmin;
  double xmax;

  CLine(){
    line[0] = '\0';
    parm[0] = '\0';
    text[0] = '\0';
    topt[0] = '\0';
    mf      = 0;
    mt      = 0;
    mtend   = 0;
    opt1    = 0;
    opt2    = 0;
    opt3    = 0;
    x       = 0.0;
    y       = 0.0;
    xmin    = 0.0;
    xmax    = 0.0;
  }
};


string   CmdGetOperation (void);
int      CmdFgetOneline (void);
string   CmdExtractArgument (void);
void     CmdExtractString (char *);
void     CmdExtractSecondString (char *);
bool     CmdCheckChar (char);

