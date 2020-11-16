/******************************************************************************/
/**                                                                          **/
/**     GFR :   Resonance Parameter Utilities                                **/
/**                                                                          **/
/******************************************************************************/

//const double gcKfactor = 0.21908;
const double gcKfactor = 0.2196771;

static inline double gfrENDFChannelRadius(int a)
{
  return(1.23*pow((double)a/1.00865,1.0/3.0) + 0.8);
}



const int MAX_RESONANCE = 5000;
const int NRANGE        =   10;
const int LMAX          =    5;
const int MAX_PAIRS     =   10;
const int MAX_XSEC      =  100;


/**************************************/
/*      Class Resonance               */
/**************************************/
class Resonance{
 public:
  int      l ;      /* resonance angular momentum */
  int     j2 ;      /* resonance spin x2 */
  double  er ;      /* resonance energy */
  double  gn ;      /* neutron channel width */
  double  gg ;      /* capture channel width */
};

/*** Breit-Wigner */
class BWResonance : public Resonance{
 public:
  double  gt ;      /* neutron channel width */
  double  gf ;      /* fission channel width */
  double  gx ;      /* other channel width */
  double   s ;      /* SL(e0) */
  double   p ;      /* PL(e0) */
};

/*** Reich-Moore */
class RMResonance : public Resonance{
 public:
  double gf1 ;      /* first fission channel width */
  double gf2 ;      /* seconf fission channel width */
  double   s ;      /* SL(e0) */
  double   p ;      /* PL(e0) */
};

/*** Unresolved Resonances */
class UBWResonance : public Resonance{
 public:
  double  gf ;      /* fission channel width */
  double  gx ;      /* other channel width */
  double   d ;      /* average energy spacing */
};

class URResonance {
 private:
  bool   allocated;
 public:
  int    l   ;
  int    j2  ;
  int    ne  ;
  int    itp ;
  double dfn ;      /* degree-of-freedom for neutron width */
  double dfg ;      /* degree-of-freedom for neutron width */
  double dff ;      /* degree-of-freedom for neutron width */
  double dfx ;      /* degree-of-freedom for neutron width */
  UBWResonance *bw;

  URResonance(){
    l      = 0;
    j2     = 0;
    ne     = 0;
    itp    = 0;
    allocated = false;
  }
  ~URResonance(){
    if(allocated){
      delete [] bw;
      allocated = false;
    }
  }

  void memalloc(int m){
    if(!allocated){
      bw = new UBWResonance [m];
      allocated = true;
    }
  }
};


/**************************************/
/*      Two Particle Pair             */
/**************************************/
class ParPair{
 public:
  double mass[2]  ; /* pair mass numbers */ 
  int    znum[2]  ; /* pair atomic numbers */
  int    spin2[2] ; /* spin x2 */
  int    parity[2]; /* parity */
  double qvalue   ; /* Q-value */
  int    fpen     ; /* flag for penetrability */
  int    fsft     ; /* flag for shift factor */
  int    mt       ; /* MT number */

  ParPair(){
    for(int i=0 ; i<2 ; i++){
      znum[i] = spin2[i] = parity[i] = 0;
      mass[i] = 0.0;
    }
    qvalue = 0.0;
    fpen   = 0;
    fsft   = 0;
    mt     = 0;
  }
};


/**************************************/
/*      R-Matrix Limited              */
/**************************************/
class RMLParameter{
 private:
  bool   allocated;
 public:
  int     j2;                 // resonance spin x2
  int     parity;             // parity
  int     nchannel;           // number of channels
  int     nresonance;         // number of resonances
  int     *l;                 // channel angular momentum
  int     *s2;                // channel spin x2
  int     *pidx;              // particle pair index
  double  *radius_effective;  // effective channel radius
  double  *radius_true;       // true channel radius
  double  *energy;            // resonance energy
  double **gamma;             // partial width

  RMLParameter(){
    j2         = 0;
    parity     = 0;
    nchannel   = 0;
    nresonance = 0;
    allocated  = false;
  }

  ~RMLParameter(){
    if(allocated){
      delete [] l;
      delete [] s2;
      delete [] pidx;
      delete [] radius_effective;
      delete [] radius_true;
      delete [] energy;
      for(int i=0 ; i<nchannel ; i++) delete [] gamma[i];
      delete [] gamma;

      allocated = false;
    }
  }

  void memalloc(int n, int m){
    if(!allocated){
      nchannel   = n;
      nresonance = m;
      l    = new int [n];
      s2   = new int [n];
      pidx = new int [n];
      radius_effective = new double [n];
      radius_true      = new double [n];
      energy = new double [m];
      gamma  = new double * [n];
      for(int i=0 ; i<nchannel ; i++) gamma[i] = new double [m];

      allocated = true;
    }
  }
};


/**************************************/
/*      RML Energy-Dependent Data     */
/**************************************/
class RMLChannel{
 public:
  double reduced_mass;        // reduced mass
  double ecms;                // CMS energy
  double wave_number;         // wave number
  double coulomb;             // Coulomb parmaeter
  double alpha_effective;     // k a(effective) 
  double alpha_true;          // k a(true)
  double mratio;              // mass ratio, (M+m)/M to convert into LAB
  bool   charge;              // flag for charged particle channel
  bool   open;                // flag for open channel

  RMLChannel(){
    reduced_mass    = 0.0;
    ecms            = 0.0;
    wave_number     = 0.0;
    coulomb         = 0.0;
    alpha_effective = 0.0;
    alpha_true      = 0.0;
    mratio          = 1.0;
    open            = true;
    charge          = false;
  }
};


/**************************************/
/*      Class System                  */
/**************************************/
class System{
 private:
  bool   first_call;          // flag for multiple energy, start new calculation
  bool   last_call;           // reset loaded data
 public:
  unsigned int target_A;
  unsigned int target_Z;
  int    target_spin2;        // target spin x 2
  int    target_parity;       // target parity
  int    incident_spin2;      // incident particle spin x 2
  double reduced_mass;        // reduced mass
  double radius;              // channel_radius [fm]
  double ecms;                // CMS energy
  double wave_number;         // wave number
  double alpha;               // alpha = kR
  int    nrange;              // number of energy range
  int    nl;                  // number of orbital angular momentum
  int    nj;                  // number of J states
  int    npair;               // number of two-particle pairs
  int    format;              // format, should be 3
  int    avefission_flag;     // averaged fission width given
  int    selfshield_flag;     // self-shielding flag
  int    gammaunit_flag;      // IFG= 0: in eV, 1: in sqrt(eV)
  int    relativ_flag;        // relativistic flag

  int    *idx;                // pointer to the data block
  int    *lru;                // resolved or unresolved
  int    *lrf;                // resonance formula
  int    *naps;               // NAPS channel radius control
  int    *nro ;               // NRO energy dependent radii flag
  double *emin;               // Emin for the range
  double *emax;               // Emax for the range

  int     nfw;                // energy-dependent fission width, energy points
  double *fwx;                // fission widths

  System(){
    first_call   = true;
    last_call    = false;
    target_A     = 0;
    target_Z     = 0;
    target_spin2 = 0;
    target_parity= 0;
    incident_spin2 = 0;
    reduced_mass = 0.0;
    radius       = 0.0;
    ecms         = 0.0;
    wave_number  = 0.0;
    alpha        = 0.0;
    nrange       = 0;
    nl           = 0;
    nj           = 0;
    npair        = 0;
    format       = 0;
    avefission_flag = 0;
    selfshield_flag = 0;
    gammaunit_flag  = 0;
    relativ_flag    = 0;

    idx  = new int [NRANGE];
    lru  = new int [NRANGE];
    lrf  = new int [NRANGE];
    naps = new int [NRANGE];
    nro  = new int [NRANGE];
    emin = new double [NRANGE];
    emax = new double [NRANGE];

    for(int i=0 ; i<NRANGE ; i++){
      idx[i] = lru[i] = lrf[i] = naps[i] = nro[i] = 0;
      emin[i] = emax[i] = 0.0;
    }
  }

  ~System(){
    delete [] idx;
    delete [] lru;
    delete [] lrf;
    delete [] naps;
    delete [] nro;
    delete [] emin;
    delete [] emax;
  }

  bool isFirstCall(){ return first_call; }
  bool isLastCall(){ return last_call; }
  void FirstCall(){ first_call = true; }
  void OnceCalled(){ first_call = false; }
  void LastCall(){ last_call = true; }
};


/**************************************/
/*      Class Channel Wave Function   */
/**************************************/
class ChannelWaveFunc{
 public:
  complex<double> H;       // Hanlel function, G+iF
  complex<double> D;       // derivative Hanlel function, G'+iF'
  complex<double> L;       // L-function, (G'+iF')/(G+iF) * rho
  complex<double> phase;   // hard-sphare phase, exp(-phi)
  complex<double> phase2;  // 2 x phase, exp(-2 phi)
  complex<double> phaseC;  // Coulomb phase, exp(-phi_c)
  double          a;       // channel radius
  double          p;       // phase

  ChannelWaveFunc(){
    H      = complex<double>(0.0,0.0);
    D      = complex<double>(0.0,0.0);
    L      = complex<double>(1.0,0.0);
    phase  = complex<double>(0.0,0.0);
    phase2 = complex<double>(0.0,0.0);
    phaseC = complex<double>(0.0,0.0);
    a      = 0.0;
    p      = 0.0;
  }

  /*** set all functions by a, H = G+iF, and H' = G'+iF' */
  void setData(const double x, complex<double> y1, complex<double> y2){
    a = x;
    H = y1;
    D = y2;
    L = a * D / H;
  }

  /*** set phase factors by H */
  void setPhase(complex<double> y){
    double x = abs(y);
    if(y != 0.0) p = acos(y.real() / x); // atan(F / G) = acos(G / sqrt(G^2+F^2))
    else         p = 0.0;
    phase  = complex<double>(cos( -p), -sin(  p)); // exp^{-ip}
    phase2 = complex<double>(cos(2*p), -sin(2*p)); // exp^{-2ip}
  }

  void setCoulombPhase(double z){
    phaseC = complex<double>(cos(  z),  sin(  z)); // exp^{iz}
  }

  /*** get penetrabilitiy and shift factor */
  double P(){ return L.imag(); }
  double S(){ return L.real(); }
};


/**************************************/
/*      Resonance Cross Sections      */
/**************************************/
class GFRcross{
 private:
  int    nch;       // number of reaction channels
  bool   allocated; // memory allocation flag
 public:
  double energy;    // LAB energy
  int    *type;     // reaction type ( = MT number )
  double *xsec;     // cross section

  GFRcross(){
    nch = 0;
    energy = 0.0;
    allocated = false;
  }

  GFRcross(const int n){
    nch = n;
    if(nch > MAX_XSEC) nch = MAX_XSEC;
    energy = 0.0;
    allocated = false;
    memalloc(nch);
    clear();
  }

  ~GFRcross(){
    memfree();
  }

  void memalloc(const int n){
    if(!allocated){
      nch = n;
      type = new int [n];
      xsec = new double [n];
      allocated = true;
      clear();
    }
  }

  void memfree(){
    if(allocated){
      delete [] type;
      delete [] xsec;
      allocated = false;
    }
  }

  void clear(){
    if(allocated){
      for(int i=0 ; i<nch ; i++){
        type[i] = 0;
        xsec[i] = 0.0;
      }
    }
  }

  void zero(){
    if(allocated){
      for(int i=0 ; i<nch ; i++) xsec[i] = 0.0;
    }
  }

  void set(const int t, const double d){
    for(int i=0 ; i<nch ; i++){
      if(type[i] == t) xsec[i] = d;
    }
  }

  void add(const int t, const double d){
    for(int i=0 ; i<nch ; i++){
      if(type[i] == t) xsec[i] += d;
    }
  }

  double get(const int t){
    double d = 0.0;
    for(int i=0 ; i<nch ; i++){
      if(type[i] == t){ d = xsec[i]; break; }
    }
    return d;
  }

  int getNch(){ return nch; }

  double sum(){
    double sum = 0.0;
    for(int i=0 ; i<nch ; i++) sum += xsec[i];
    return sum;
  }
};


/**************************************/
/*      One Point Cross Sections      */
/**************************************/
class Pcross{
 public:
  double energy;
  double total;
  double elastic;
  double capture;
  double fission;
  double proton;
  double alpha;
  double other;

  Pcross(){
    clear();
  }

  void clear(){
    energy   = 0.0;
    total    = 0.0;
    elastic  = 0.0;
    capture  = 0.0;
    fission  = 0.0;
    proton   = 0.0;
    alpha    = 0.0;
    other    = 0.0;
  }

  void setTotal(){
    total = elastic + capture + fission + proton + alpha + other;
  }

  Pcross operator+(Pcross x){
    Pcross y;
    y.total   = total   + x.total;
    y.elastic = elastic + x.elastic;
    y.capture = capture + x.capture;
    y.fission = fission + x.fission;
    y.proton  = proton  + x.proton;
    y.alpha   = alpha   + x.alpha;
    y.other   = other   + y.other;
    return y;
  }

  Pcross operator-(Pcross x){
    Pcross y;
    y.total   = total   - x.total;
    y.elastic = elastic - x.elastic;
    y.capture = capture - x.capture;
    y.fission = fission - x.fission;
    y.proton  = proton  - x.proton;
    y.alpha   = alpha   - x.alpha;
    y.other   = other   - y.other;
    return y;
  }
};


/**************************************/
/*      Scattering Channel S-matrix   */
/**************************************/
class SElement{
 public:
  int l  ;
  int s2 ;
  int j2 ;
  complex<double> S;

  SElement(){
    l  = 0;
    s2 = 0;
    j2 = 0;
    S  = complex<double>(0.0,0.0);
  }
};

class Smatrix{
 private:
  int  size;
  int  index;
  bool allocated;
 public:
  SElement *element;

  Smatrix(){
    allocated = false;
    size  = 0;
    index = 0;
  }
  ~Smatrix(){
    if(allocated){ memfree(); }
  }

  void memalloc(int n){
    if(!allocated){
      element = new SElement [n];
      allocated = true;
      for(int i=0 ; i<n ; i++){
        element[i].l = element[i].s2 = element[i].j2 = 0;
        element[i].S = complex<double>(0.0,0.0);
      }
      size = n;
    }
  }

  void memfree(){
    if(allocated){
      delete [] element;
      allocated = false;
    }
  }

  void setIndex   (int k){ index = k; }
  int  getIndex   ( )    { return(index); }
  void resetIndex ( )    { index = 0; }
  void inclIndex  ( )    { index ++;  }

  int findIndex(int pl, int pj2, int ps2){
    int id = -1;
    if(allocated){
      for(int i=0 ; i<size ; i++){
        if(element[i].s2 == 0) break;
        if( (element[i].l  == pl ) &&
            (element[i].s2 == ps2) &&
            (element[i].j2 == pj2) ){
          id = i;
          break;
        }
      }
    }
    return(id);
  }

  void setElement(int pl, int pj2, int ps2, complex<double> z){
    if(allocated && (index < size)){
      element[index].l  = pl;
      element[index].j2 = pj2;
      element[index].s2 = ps2;
      element[index].S  = z;
    }
  }

  complex<double> getElement(int k){
    complex<double> s(0.0,0.0);
    if(allocated && (k < size)){
      s = element[k].S;
    }
    return(s);
  }

  complex<double> getElement(int k, int *pl, int *pj2, int *ps2){
    complex<double> s(0.0,0.0);
    if(allocated && (k < size)){
      s = element[k].S;
      *pl  = element[k].l;
      *pj2 = element[k].j2;
      *ps2 = element[k].s2;
    }
    return(s);
  }
};


/**************************************/
/*      gfr.cpp                       */
/**************************************/
void    gfrReadHEADData           (System *, ENDF *);


/**************************************/
/*      gfrcross.cpp                  */
/**************************************/
Pcross  gfrCrossSection           (const int, const double, System *, ENDF *);
void    gfrSetEnergy              (const double, System *);
void    gfrPenetrability          (const int, const double, ChannelWaveFunc *);
complex<double> gfrLfunction      (const int, const double, const double);


/**************************************/
/*      gfrcs1.cpp                    */
/**************************************/
Pcross  gfrCrossSection1          (const int, const int, const double, System *, ENDF *);


/**************************************/
/*      gfrcs3.cpp                    */
/**************************************/
Pcross  gfrCrossSection3          (const int, const double, System *, ENDF *);


/**************************************/
/*      gfrcs7.cpp                    */
/**************************************/
Pcross  gfrCrossSection7         (const int, const double, System *, ENDF *);


/**************************************/
/*      gfrcsurr.cpp                  */
/**************************************/
Pcross  gfrCrossSectionURR        (const int, const double, System *, ENDF *);


/**************************************/
/*      gfrformula.cpp                */
/**************************************/
Pcross  gfrSLBreitWigner          (const int, const int, const int, const double, ChannelWaveFunc *, BWResonance *);
Pcross  gfrMLBreitWignerENDF      (const int, const int, const int, const int, const double, ChannelWaveFunc *, BWResonance *);
Pcross  gfrMLBreitWigner          (const int, const int, const int, const double, ChannelWaveFunc *, BWResonance *);

Pcross  gfrBreitWignerUmatrix     (const int, const int, const int, const int, const double, ChannelWaveFunc *, BWResonance *, const int);
Pcross  gfrReichMoore             (const int, const int, const int, const int, const int, const double, ChannelWaveFunc *, RMResonance *);


/**************************************/
/*      gfrlefcoef.cpp                */
/**************************************/
void    gfrLegendreCoefficient    (System *, double *);
double  gfrCompoundReaction       (System *);


/**************************************/
/*      gfrenergy.cpp                 */
/**************************************/
int     gfrAutoEnergyRRR          (System *, ENDF *, double *, const double);
int     gfrAutoEnergyURR          (double *, const double, const double);
int     gfrFixedEnergyRRR         (double, double, double, double *, const double, const double);
