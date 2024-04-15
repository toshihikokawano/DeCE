/*
   masstable.h : 
        class and prototypes of mass excess data
 */

class MassExcess{
 public:
  unsigned int za;    // Z*1000 + A
  float        mass;  // mass excess
};

/**************************************/
/*      masstable.cpp                 */
/**************************************/
double  mass_excess           (int, int);
double  mass_qvalue           (int, int, int);
double  mass_threshold        (int, double);
