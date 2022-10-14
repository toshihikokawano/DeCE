/**********************************************************/
/*   Z and A numbers of Nucleus                           */
/**********************************************************/
class ZAnumber{ 
 private:
    unsigned int Z;
    unsigned int A;
 public:
    ZAnumber(){
      Z = 0;
      A = 0;
    }
    ZAnumber(int z, int a){
      Z = z;
      A = a;
    }
    void setZA(int z, int a){
      Z = z;
      A = a;
    }
    unsigned int getZ(){ return (Z); }
    unsigned int getA(){ return (A); }
    unsigned int getN(){ return (A-Z); }
    ZAnumber operator+(ZAnumber x){
      ZAnumber y;
      y.Z = Z + x.Z;
      y.A = A + x.A;
      return y;
    }
    ZAnumber operator-(ZAnumber x){
      ZAnumber y;
      y.Z = Z - x.Z;
      y.A = A - x.A;
      return y;
    }
    bool operator==(ZAnumber x){
      bool z = false;
      if( (Z == x.Z) && (A == x.A) ) z = true;
      return z;
    }
    bool operator!=(ZAnumber x){
      bool z = false;
      if( (Z != x.Z) || (A != x.A) ) z = true;
      return z;
    }
};



/*** kalbach.cpp */
void   ddxKalbach (const int, const double, const double, const double, const double, double *, double *);
void   ddxKalbachSetParm (const double, const double, const double);

