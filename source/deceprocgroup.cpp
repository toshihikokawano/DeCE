/******************************************************************************/
/**     DeCE Proccessing: Group Cross Section                                **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <complex>

using namespace std;

#include "dece.h"
#include "gfr.h"
#include "terminate.h"
#include "groupstructure_SANDIIa.h"

static const int ncx = 14;
static int mtr[ncx] = {1, 2, 4, 16, 17, 18, 22, 28, 102, 103, 104, 105, 106, 107};

static void DeceGroupAverage (ENDF *lib, const int weight);

/**********************************************************/
/*      Generate Group Cross Section                      */
/*      --------                                          */
/*      */
/**********************************************************/
void DeceGenerateGroup(ENDFDict *dict, ENDF *lib[], const int weight)
{
  double **xdat;
  xdat = new double * [ncx];
  for(int j=0 ; j<ncx ; j++){
    xdat[j] = new double [grpEnergyPoint];
  }

  int id = dict->getID(3,1);
  if(id < 0){
    message << "MF3 MT1 should exit for processing";
    WarningMessage();
    return;
  }

  /*** excluded channels' MT number negative */
  for(int j=1 ; j<ncx ; j++){
    if(dict->getID(3,mtr[j]) < 0) mtr[j] *= -1;
  }

  /*** group average */
  for(int j=0 ; j<ncx ; j++){
    id = dict->getID(3,mtr[j]);
    if(id > 0) DeceGroupAverage(lib[id],weight);
  }


  for(int j=0 ; j<ncx ; j++){
    delete [] xdat[j];
  }
  delete [] xdat;
}


void DeceGroupAverage(ENDF *lib, const int w)
{
  for(int i=0 ; i<grpEnergyPoint-1 ; i++){
    double e0 = grpEnergyGrid[i];
    double e1 = grpEnergyGrid[i+1];
  }
}
