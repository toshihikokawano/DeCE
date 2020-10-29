/******************************************************************************/
/**                                                                          **/
/**     ENDF TEXT : TEXT Data Definition                                     **/
/**                                                                          **/
/******************************************************************************/

#ifndef TEXT_WIDTH
#define TEXT_WIDTH 66
#endif

#ifndef __CSTRING_H__
#include <cstring>
#define __CSTRING_H__
#endif

/**************************************/
/*      Class ENDF : Text Field       */
/**************************************/
class ENDFText{
  private:
    bool allocated;     // flag for memory allocation
    int  length;        // text data length
    char *text;         // data content
  public:
    int xpos, ypos;     // location (X-line,Y-column) of the text data

  ENDFText(const int m, const int x, const int y){
    length = m;
    text = new char [m+1];
    for(int i=0 ; i<length ; i++) text[i] = ' ';
    text[length] = '\0';
    xpos = x;
    ypos = y;
    allocated = true;
  }

  ~ENDFText(){
    if(allocated){
      delete [] text;
    }
    allocated = false;
  }

  int getlen(){ return length; }

  void read(char *src){
    int c = strlen(src);
    if(c >= length){
      for(int i=0 ; i<length ; i++) text[i] = src[i];
    }
    else{
      for(int i=0 ; i<c ; i++) text[i] = src[i];
      for(int i=c ; i<length ; i++) text[i] = ' ';
    }
  }

  void copy(char *src){
    for(int i=0 ; i<length ; i++) text[i] = src[i+xpos];
  }

  void paste(int c, char *dst){
    int ix = length;
    if(c < length + xpos) ix = c - xpos;
    for(int i=0 ; i<ix ; i++){
      char c = text[i];
      if(isprint(c) == 0) c = ' ';
      dst[i+xpos] = c;
    }
  }

  void print(){
    int ix = length-1;
    for(int i=ix ; i>=0 ; i--){
      if(text[i] != ' '){ ix = i; break; }
    }
    for(int i=0 ; i<=ix ; i++) cout << text[i];
  }

  void getText(){
    for(int i=0 ; i<length ; i++) cout << text[i];
  }
};



