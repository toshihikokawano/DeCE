/******************************************************************************/
/**                                                                          **/
/**     Common Group Structures                                              **/
/**             0: SAND-IIa dosimetry 640 group                              **/
/**             1: LANL 70 group                                             **/
/**             2: VITAMIN-J 175 group                                       **/
/**             3: SAND-IIa dosimetry 725 group                              **/
/**             4: LANL 618 group                                            **/
/**                                                                          **/
/******************************************************************************/

/*** SAND-IIA */
static const string grpStructureName0 = "SAND-IIa 640";
static const int grpEnergyPoint0 = 641;
static double grpEnergyGrid0[grpEnergyPoint0] = {
 1.0000e-04, 1.0500e-04, 1.1000e-04, 1.1500e-04, 1.2000e-04, 1.2750e-04, 1.3500e-04, 1.4250e-04, 1.5000e-04, 1.6000e-04,
 1.7000e-04, 1.8000e-04, 1.9000e-04, 2.0000e-04, 2.1000e-04, 2.2000e-04, 2.3000e-04, 2.4000e-04, 2.5500e-04, 2.7000e-04,
 2.8000e-04, 3.0000e-04, 3.2000e-04, 3.4000e-04, 3.6000e-04, 3.8000e-04, 4.0000e-04, 4.2500e-04, 4.5000e-04, 4.7500e-04,
 5.0000e-04, 5.2500e-04, 5.5000e-04, 5.7500e-04, 6.0000e-04, 6.3000e-04, 6.6000e-04, 6.9000e-04, 7.2000e-04, 7.6000e-04,
 8.0000e-04, 8.4000e-04, 8.8000e-04, 9.2000e-04, 9.6000e-04, 1.0000e-03, 1.0500e-03, 1.1000e-03, 1.1500e-03, 1.2000e-03,
 1.2750e-03, 1.3500e-03, 1.4250e-03, 1.5000e-03, 1.6000e-03, 1.7000e-03, 1.8000e-03, 1.9000e-03, 2.0000e-03, 2.1000e-03,
 2.2000e-03, 2.3000e-03, 2.4000e-03, 2.5500e-03, 2.7000e-03, 2.8000e-03, 3.0000e-03, 3.2000e-03, 3.4000e-03, 3.6000e-03,
 3.8000e-03, 4.0000e-03, 4.2500e-03, 4.5000e-03, 4.7500e-03, 5.0000e-03, 5.2500e-03, 5.5000e-03, 5.7500e-03, 6.0000e-03,
 6.3000e-03, 6.6000e-03, 6.9000e-03, 7.2000e-03, 7.6000e-03, 8.0000e-03, 8.4000e-03, 8.8000e-03, 9.2000e-03, 9.6000e-03,
 1.0000e-02, 1.0500e-02, 1.1000e-02, 1.1500e-02, 1.2000e-02, 1.2750e-02, 1.3500e-02, 1.4250e-02, 1.5000e-02, 1.6000e-02,
 1.7000e-02, 1.8000e-02, 1.9000e-02, 2.0000e-02, 2.1000e-02, 2.2000e-02, 2.3000e-02, 2.4000e-02, 2.5500e-02, 2.7000e-02,
 2.8000e-02, 3.0000e-02, 3.2000e-02, 3.4000e-02, 3.6000e-02, 3.8000e-02, 4.0000e-02, 4.2500e-02, 4.5000e-02, 4.7500e-02,
 5.0000e-02, 5.2500e-02, 5.5000e-02, 5.7500e-02, 6.0000e-02, 6.3000e-02, 6.6000e-02, 6.9000e-02, 7.2000e-02, 7.6000e-02,
 8.0000e-02, 8.4000e-02, 8.8000e-02, 9.2000e-02, 9.6000e-02, 1.0000e-01, 1.0500e-01, 1.1000e-01, 1.1500e-01, 1.2000e-01,
 1.2750e-01, 1.3500e-01, 1.4250e-01, 1.5000e-01, 1.6000e-01, 1.7000e-01, 1.8000e-01, 1.9000e-01, 2.0000e-01, 2.1000e-01,
 2.2000e-01, 2.3000e-01, 2.4000e-01, 2.5500e-01, 2.7000e-01, 2.8000e-01, 3.0000e-01, 3.2000e-01, 3.4000e-01, 3.6000e-01,
 3.8000e-01, 4.0000e-01, 4.2500e-01, 4.5000e-01, 4.7500e-01, 5.0000e-01, 5.2500e-01, 5.5000e-01, 5.7500e-01, 6.0000e-01,
 6.3000e-01, 6.6000e-01, 6.9000e-01, 7.2000e-01, 7.6000e-01, 8.0000e-01, 8.4000e-01, 8.8000e-01, 9.2000e-01, 9.6000e-01,
 1.0000e+00, 1.0500e+00, 1.1000e+00, 1.1500e+00, 1.2000e+00, 1.2750e+00, 1.3500e+00, 1.4250e+00, 1.5000e+00, 1.6000e+00,
 1.7000e+00, 1.8000e+00, 1.9000e+00, 2.0000e+00, 2.1000e+00, 2.2000e+00, 2.3000e+00, 2.4000e+00, 2.5500e+00, 2.7000e+00,
 2.8000e+00, 3.0000e+00, 3.2000e+00, 3.4000e+00, 3.6000e+00, 3.8000e+00, 4.0000e+00, 4.2500e+00, 4.5000e+00, 4.7500e+00,
 5.0000e+00, 5.2500e+00, 5.5000e+00, 5.7500e+00, 6.0000e+00, 6.3000e+00, 6.6000e+00, 6.9000e+00, 7.2000e+00, 7.6000e+00,
 8.0000e+00, 8.4000e+00, 8.8000e+00, 9.2000e+00, 9.6000e+00, 1.0000e+01, 1.0500e+01, 1.1000e+01, 1.1500e+01, 1.2000e+01,
 1.2750e+01, 1.3500e+01, 1.4250e+01, 1.5000e+01, 1.6000e+01, 1.7000e+01, 1.8000e+01, 1.9000e+01, 2.0000e+01, 2.1000e+01,
 2.2000e+01, 2.3000e+01, 2.4000e+01, 2.5500e+01, 2.7000e+01, 2.8000e+01, 3.0000e+01, 3.2000e+01, 3.4000e+01, 3.6000e+01,
 3.8000e+01, 4.0000e+01, 4.2500e+01, 4.5000e+01, 4.7500e+01, 5.0000e+01, 5.2500e+01, 5.5000e+01, 5.7500e+01, 6.0000e+01,
 6.3000e+01, 6.6000e+01, 6.9000e+01, 7.2000e+01, 7.6000e+01, 8.0000e+01, 8.4000e+01, 8.8000e+01, 9.2000e+01, 9.6000e+01,
 1.0000e+02, 1.0500e+02, 1.1000e+02, 1.1500e+02, 1.2000e+02, 1.2750e+02, 1.3500e+02, 1.4250e+02, 1.5000e+02, 1.6000e+02,
 1.7000e+02, 1.8000e+02, 1.9000e+02, 2.0000e+02, 2.1000e+02, 2.2000e+02, 2.3000e+02, 2.4000e+02, 2.5500e+02, 2.7000e+02,
 2.8000e+02, 3.0000e+02, 3.2000e+02, 3.4000e+02, 3.6000e+02, 3.8000e+02, 4.0000e+02, 4.2500e+02, 4.5000e+02, 4.7500e+02,
 5.0000e+02, 5.2500e+02, 5.5000e+02, 5.7500e+02, 6.0000e+02, 6.3000e+02, 6.6000e+02, 6.9000e+02, 7.2000e+02, 7.6000e+02,
 8.0000e+02, 8.4000e+02, 8.8000e+02, 9.2000e+02, 9.6000e+02, 1.0000e+03, 1.0500e+03, 1.1000e+03, 1.1500e+03, 1.2000e+03,
 1.2750e+03, 1.3500e+03, 1.4250e+03, 1.5000e+03, 1.6000e+03, 1.7000e+03, 1.8000e+03, 1.9000e+03, 2.0000e+03, 2.1000e+03,
 2.2000e+03, 2.3000e+03, 2.4000e+03, 2.5500e+03, 2.7000e+03, 2.8000e+03, 3.0000e+03, 3.2000e+03, 3.4000e+03, 3.6000e+03,
 3.8000e+03, 4.0000e+03, 4.2500e+03, 4.5000e+03, 4.7500e+03, 5.0000e+03, 5.2500e+03, 5.5000e+03, 5.7500e+03, 6.0000e+03,
 6.3000e+03, 6.6000e+03, 6.9000e+03, 7.2000e+03, 7.6000e+03, 8.0000e+03, 8.4000e+03, 8.8000e+03, 9.2000e+03, 9.6000e+03,
 1.0000e+04, 1.0500e+04, 1.1000e+04, 1.1500e+04, 1.2000e+04, 1.2750e+04, 1.3500e+04, 1.4250e+04, 1.5000e+04, 1.6000e+04,
 1.7000e+04, 1.8000e+04, 1.9000e+04, 2.0000e+04, 2.1000e+04, 2.2000e+04, 2.3000e+04, 2.4000e+04, 2.5500e+04, 2.7000e+04,
 2.8000e+04, 3.0000e+04, 3.2000e+04, 3.4000e+04, 3.6000e+04, 3.8000e+04, 4.0000e+04, 4.2500e+04, 4.5000e+04, 4.7500e+04,
 5.0000e+04, 5.2500e+04, 5.5000e+04, 5.7500e+04, 6.0000e+04, 6.3000e+04, 6.6000e+04, 6.9000e+04, 7.2000e+04, 7.6000e+04,
 8.0000e+04, 8.4000e+04, 8.8000e+04, 9.2000e+04, 9.6000e+04, 1.0000e+05, 1.0500e+05, 1.1000e+05, 1.1500e+05, 1.2000e+05,
 1.2750e+05, 1.3500e+05, 1.4250e+05, 1.5000e+05, 1.6000e+05, 1.7000e+05, 1.8000e+05, 1.9000e+05, 2.0000e+05, 2.1000e+05,
 2.2000e+05, 2.3000e+05, 2.4000e+05, 2.5500e+05, 2.7000e+05, 2.8000e+05, 3.0000e+05, 3.2000e+05, 3.4000e+05, 3.6000e+05,
 3.8000e+05, 4.0000e+05, 4.2500e+05, 4.5000e+05, 4.7500e+05, 5.0000e+05, 5.2500e+05, 5.5000e+05, 5.7500e+05, 6.0000e+05,
 6.3000e+05, 6.6000e+05, 6.9000e+05, 7.2000e+05, 7.6000e+05, 8.0000e+05, 8.4000e+05, 8.8000e+05, 9.2000e+05, 9.6000e+05,
 1.0000e+06, 1.1000e+06, 1.2000e+06, 1.3000e+06, 1.4000e+06, 1.5000e+06, 1.6000e+06, 1.7000e+06, 1.8000e+06, 1.9000e+06,
 2.0000e+06, 2.1000e+06, 2.2000e+06, 2.3000e+06, 2.4000e+06, 2.5000e+06, 2.6000e+06, 2.7000e+06, 2.8000e+06, 2.9000e+06,
 3.0000e+06, 3.1000e+06, 3.2000e+06, 3.3000e+06, 3.4000e+06, 3.5000e+06, 3.6000e+06, 3.7000e+06, 3.8000e+06, 3.9000e+06,
 4.0000e+06, 4.1000e+06, 4.2000e+06, 4.3000e+06, 4.4000e+06, 4.5000e+06, 4.6000e+06, 4.7000e+06, 4.8000e+06, 4.9000e+06,
 5.0000e+06, 5.1000e+06, 5.2000e+06, 5.3000e+06, 5.4000e+06, 5.5000e+06, 5.6000e+06, 5.7000e+06, 5.8000e+06, 5.9000e+06,
 6.0000e+06, 6.1000e+06, 6.2000e+06, 6.3000e+06, 6.4000e+06, 6.5000e+06, 6.6000e+06, 6.7000e+06, 6.8000e+06, 6.9000e+06,
 7.0000e+06, 7.1000e+06, 7.2000e+06, 7.3000e+06, 7.4000e+06, 7.5000e+06, 7.6000e+06, 7.7000e+06, 7.8000e+06, 7.9000e+06,
 8.0000e+06, 8.1000e+06, 8.2000e+06, 8.3000e+06, 8.4000e+06, 8.5000e+06, 8.6000e+06, 8.7000e+06, 8.8000e+06, 8.9000e+06,
 9.0000e+06, 9.1000e+06, 9.2000e+06, 9.3000e+06, 9.4000e+06, 9.5000e+06, 9.6000e+06, 9.7000e+06, 9.8000e+06, 9.9000e+06,
 1.0000e+07, 1.0100e+07, 1.0200e+07, 1.0300e+07, 1.0400e+07, 1.0500e+07, 1.0600e+07, 1.0700e+07, 1.0800e+07, 1.0900e+07,
 1.1000e+07, 1.1100e+07, 1.1200e+07, 1.1300e+07, 1.1400e+07, 1.1500e+07, 1.1600e+07, 1.1700e+07, 1.1800e+07, 1.1900e+07,
 1.2000e+07, 1.2100e+07, 1.2200e+07, 1.2300e+07, 1.2400e+07, 1.2500e+07, 1.2600e+07, 1.2700e+07, 1.2800e+07, 1.2900e+07,
 1.3000e+07, 1.3100e+07, 1.3200e+07, 1.3300e+07, 1.3400e+07, 1.3500e+07, 1.3600e+07, 1.3700e+07, 1.3800e+07, 1.3900e+07,
 1.4000e+07, 1.4100e+07, 1.4200e+07, 1.4300e+07, 1.4400e+07, 1.4500e+07, 1.4600e+07, 1.4700e+07, 1.4800e+07, 1.4900e+07,
 1.5000e+07, 1.5100e+07, 1.5200e+07, 1.5300e+07, 1.5400e+07, 1.5500e+07, 1.5600e+07, 1.5700e+07, 1.5800e+07, 1.5900e+07,
 1.6000e+07, 1.6100e+07, 1.6200e+07, 1.6300e+07, 1.6400e+07, 1.6500e+07, 1.6600e+07, 1.6700e+07, 1.6800e+07, 1.6900e+07,
 1.7000e+07, 1.7100e+07, 1.7200e+07, 1.7300e+07, 1.7400e+07, 1.7500e+07, 1.7600e+07, 1.7700e+07, 1.7800e+07, 1.7900e+07,
 1.8000e+07, 1.8100e+07, 1.8200e+07, 1.8300e+07, 1.8400e+07, 1.8500e+07, 1.8600e+07, 1.8700e+07, 1.8800e+07, 1.8900e+07,
 1.9000e+07, 1.9100e+07, 1.9200e+07, 1.9300e+07, 1.9400e+07, 1.9500e+07, 1.9600e+07, 1.9700e+07, 1.9800e+07, 1.9900e+07,
 2.0000e+07};

/*** LANL 70 */
static const string grpStructureName1 = "LANL70";
static const int grpEnergyPoint1 =  71;
static double grpEnergyGrid1[grpEnergyPoint1] = {
      10.677,      61.4421,    101.301,     130.073,     167.017,     214.454,     275.365,     353.575,     453.999,     582.947,
     748.518,     961.117,    1089.090,    1234.100,    1398.420,    1584.610,    1795.600,    2034.680,    2305.600,    2612.590,
    2960.450,    3354.630,    3801.290,    4307.430,    4880.950,    5530.840,    6267.270,    7101.740,    8047.330,    9118.820,
   10333.300,   11708.800,   13267.800,   15034.400,   17036.200,   19304.500,   21874.900,   24787.500,   28087.900,   31827.800,
   40867.700,   52475.200,   67379.500,   86517.000,  111090.000,  142642.000,  183156.000,  235178.000,  301974.000,  387742.000,
  439369.000,  497871.000,  564161.000,  639279.000,  724398.000,  820850.000,  930145.000, 1053990.000, 1194330.000, 1353350.000,
 1737740.000, 2231300.000, 2865050.000, 3678790.000, 4723670.000, 6065310.000, 7788010.000, 1.00000e+07, 1.28403e+07, 1.64872e+07,
 2.00000e+07};

/*** VITAMINE-J */
static const string grpStructureName2 = "VITAMINE-J 175";
static const int grpEnergyPoint2 =  176;
static double grpEnergyGrid2[grpEnergyPoint2] = {
 1.0000E-05, 1.0000E-01, 4.1399E-01, 5.3158E-01, 6.8256E-01, 8.7642E-01, 1.1254E+00, 1.4450E+00, 1.8554E+00, 2.3824E+00,
 3.0590E+00, 3.9279E+00, 5.0435E+00, 6.4760E+00, 8.3153E+00, 1.0677E+01, 1.3710E+01, 1.7603E+01, 2.2603E+01, 2.9023E+01,
 3.7267E+01, 4.7851E+01, 6.1442E+01, 7.8893E+01, 1.0130E+02, 1.3007E+02, 1.6702E+02, 2.1445E+02, 2.7536E+02, 3.5358E+02,
 4.5400E+02, 5.8295E+02, 7.4852E+02, 9.6112E+02, 1.2341E+03, 1.5846E+03, 2.0347E+03, 2.2487E+03, 2.4852E+03, 2.6126E+03,
 2.7465E+03, 3.0354E+03, 3.3546E+03, 3.7074E+03, 4.3074E+03, 5.5308E+03, 7.1017E+03, 9.1188E+03, 1.0595E+04, 1.1709E+04,
 1.5034E+04, 1.9305E+04, 2.1875E+04, 2.3579E+04, 2.4176E+04, 2.4788E+04, 2.6058E+04, 2.7000E+04, 2.8500E+04, 3.1828E+04,
 3.4307E+04, 4.0868E+04, 4.6309E+04, 5.2475E+04, 5.6562E+04, 6.7379E+04, 7.2000E+04, 7.9500E+04, 8.2500E+04, 8.6517E+04,
 9.8037E+04, 1.1109E+05, 1.1679E+05, 1.2277E+05, 1.2907E+05, 1.3569E+05, 1.4264E+05, 1.4996E+05, 1.5764E+05, 1.6573E+05,
 1.7422E+05, 1.8316E+05, 1.9255E+05, 2.0242E+05, 2.1280E+05, 2.2371E+05, 2.3518E+05, 2.4724E+05, 2.7324E+05, 2.8725E+05,
 2.9452E+05, 2.9720E+05, 2.9850E+05, 3.0197E+05, 3.3373E+05, 3.6883E+05, 3.8774E+05, 4.0762E+05, 4.5049E+05, 4.9787E+05,
 5.2340E+05, 5.5023E+05, 5.7844E+05, 6.0810E+05, 6.3928E+05, 6.7206E+05, 7.0651E+05, 7.4274E+05, 7.8082E+05, 8.2085E+05,
 8.6294E+05, 9.0718E+05, 9.6164E+05, 1.0026E+06, 1.1080E+06, 1.1648E+06, 1.2246E+06, 1.2873E+06, 1.3534E+06, 1.4227E+06,
 1.4957E+06, 1.5724E+06, 1.6530E+06, 1.7377E+06, 1.8268E+06, 1.9205E+06, 2.0190E+06, 2.1225E+06, 2.2313E+06, 2.3069E+06,
 2.3457E+06, 2.3653E+06, 2.3852E+06, 2.4660E+06, 2.5924E+06, 2.7253E+06, 2.8650E+06, 3.0119E+06, 3.1664E+06, 3.3287E+06,
 3.6788E+06, 4.0657E+06, 4.4933E+06, 4.7237E+06, 4.9659E+06, 5.2205E+06, 5.4881E+06, 5.7695E+06, 6.0653E+06, 6.3763E+06,
 6.5924E+06, 6.7032E+06, 7.0469E+06, 7.4082E+06, 7.7880E+06, 8.1873E+06, 8.6071E+06, 9.0484E+06, 9.5123E+06, 1.0000E+07,
 1.0513E+07, 1.1052E+07, 1.1618E+07, 1.2214E+07, 1.2523E+07, 1.2840E+07, 1.3499E+07, 1.3840E+07, 1.4191E+07, 1.4550E+07,
 1.4918E+07, 1.5683E+07, 1.6487E+07, 1.6905E+07, 1.7333E+07, 1.9640E+07};

/*** SAND-IIA 725 */
static const string grpStructureName3 = "SAND-IIa 725";
static const int grpEnergyPoint3 = 726;
static double grpEnergyGrid3[grpEnergyPoint3] = {
 1.0000e-05, 1.0500e-05, 1.1000e-05, 1.1500e-05, 1.2000e-05, 1.2750e-05, 1.3500e-05, 1.4250e-05, 1.5000e-05, 1.6000e-05, 
 1.7000e-05, 1.8000e-05, 1.9000e-05, 2.0000e-05, 2.1000e-05, 2.2000e-05, 2.3000e-05, 2.4000e-05, 2.5500e-05, 2.7000e-05, 
 2.8000e-05, 3.0000e-05, 3.2000e-05, 3.4000e-05, 3.6000e-05, 3.8000e-05, 4.0000e-05, 4.2500e-05, 4.5000e-05, 4.7500e-05, 
 5.0000e-05, 5.2500e-05, 5.5000e-05, 5.7500e-05, 6.0000e-05, 6.3000e-05, 6.6000e-05, 6.9000e-05, 7.2000e-05, 7.6000e-05, 
 8.0000e-05, 8.4000e-05, 8.8000e-05, 9.2000e-05, 9.6000e-05, 1.0000e-04, 1.0500e-04, 1.1000e-04, 1.1500e-04, 1.2000e-04, 
 1.2750e-04, 1.3500e-04, 1.4250e-04, 1.5000e-04, 1.6000e-04, 1.7000e-04, 1.8000e-04, 1.9000e-04, 2.0000e-04, 2.1000e-04, 
 2.2000e-04, 2.3000e-04, 2.4000e-04, 2.5500e-04, 2.7000e-04, 2.8000e-04, 3.0000e-04, 3.2000e-04, 3.4000e-04, 3.6000e-04, 
 3.8000e-04, 4.0000e-04, 4.2500e-04, 4.5000e-04, 4.7500e-04, 5.0000e-04, 5.2500e-04, 5.5000e-04, 5.7500e-04, 6.0000e-04, 
 6.3000e-04, 6.6000e-04, 6.9000e-04, 7.2000e-04, 7.6000e-04, 8.0000e-04, 8.4000e-04, 8.8000e-04, 9.2000e-04, 9.6000e-04, 
 1.0000e-03, 1.0500e-03, 1.1000e-03, 1.1500e-03, 1.2000e-03, 1.2750e-03, 1.3500e-03, 1.4250e-03, 1.5000e-03, 1.6000e-03, 
 1.7000e-03, 1.8000e-03, 1.9000e-03, 2.0000e-03, 2.1000e-03, 2.2000e-03, 2.3000e-03, 2.4000e-03, 2.5500e-03, 2.7000e-03, 
 2.8000e-03, 3.0000e-03, 3.2000e-03, 3.4000e-03, 3.6000e-03, 3.8000e-03, 4.0000e-03, 4.2500e-03, 4.5000e-03, 4.7500e-03, 
 5.0000e-03, 5.2500e-03, 5.5000e-03, 5.7500e-03, 6.0000e-03, 6.3000e-03, 6.6000e-03, 6.9000e-03, 7.2000e-03, 7.6000e-03, 
 8.0000e-03, 8.4000e-03, 8.8000e-03, 9.2000e-03, 9.6000e-03, 1.0000e-02, 1.0500e-02, 1.1000e-02, 1.1500e-02, 1.2000e-02, 
 1.2750e-02, 1.3500e-02, 1.4250e-02, 1.5000e-02, 1.6000e-02, 1.7000e-02, 1.8000e-02, 1.9000e-02, 2.0000e-02, 2.1000e-02, 
 2.2000e-02, 2.3000e-02, 2.4000e-02, 2.5500e-02, 2.7000e-02, 2.8000e-02, 3.0000e-02, 3.2000e-02, 3.4000e-02, 3.6000e-02, 
 3.8000e-02, 4.0000e-02, 4.2500e-02, 4.5000e-02, 4.7500e-02, 5.0000e-02, 5.2500e-02, 5.5000e-02, 5.7500e-02, 6.0000e-02, 
 6.3000e-02, 6.6000e-02, 6.9000e-02, 7.2000e-02, 7.6000e-02, 8.0000e-02, 8.4000e-02, 8.8000e-02, 9.2000e-02, 9.6000e-02, 
 1.0000e-01, 1.0500e-01, 1.1000e-01, 1.1500e-01, 1.2000e-01, 1.2750e-01, 1.3500e-01, 1.4250e-01, 1.5000e-01, 1.6000e-01, 
 1.7000e-01, 1.8000e-01, 1.9000e-01, 2.0000e-01, 2.1000e-01, 2.2000e-01, 2.3000e-01, 2.4000e-01, 2.5500e-01, 2.7000e-01, 
 2.8000e-01, 3.0000e-01, 3.2000e-01, 3.4000e-01, 3.6000e-01, 3.8000e-01, 4.0000e-01, 4.2500e-01, 4.5000e-01, 4.7500e-01, 
 5.0000e-01, 5.2500e-01, 5.5000e-01, 5.7500e-01, 6.0000e-01, 6.3000e-01, 6.6000e-01, 6.9000e-01, 7.2000e-01, 7.6000e-01, 
 8.0000e-01, 8.4000e-01, 8.8000e-01, 9.2000e-01, 9.6000e-01, 1.0000e+00, 1.0500e+00, 1.1000e+00, 1.1500e+00, 1.2000e+00, 
 1.2750e+00, 1.3500e+00, 1.4250e+00, 1.5000e+00, 1.6000e+00, 1.7000e+00, 1.8000e+00, 1.9000e+00, 2.0000e+00, 2.1000e+00, 
 2.2000e+00, 2.3000e+00, 2.4000e+00, 2.5500e+00, 2.7000e+00, 2.8000e+00, 3.0000e+00, 3.2000e+00, 3.4000e+00, 3.6000e+00, 
 3.8000e+00, 4.0000e+00, 4.2500e+00, 4.5000e+00, 4.7500e+00, 5.0000e+00, 5.2500e+00, 5.5000e+00, 5.7500e+00, 6.0000e+00, 
 6.3000e+00, 6.6000e+00, 6.9000e+00, 7.2000e+00, 7.6000e+00, 8.0000e+00, 8.4000e+00, 8.8000e+00, 9.2000e+00, 9.6000e+00, 
 1.0000e+01, 1.0500e+01, 1.1000e+01, 1.1500e+01, 1.2000e+01, 1.2750e+01, 1.3500e+01, 1.4250e+01, 1.5000e+01, 1.6000e+01, 
 1.7000e+01, 1.8000e+01, 1.9000e+01, 2.0000e+01, 2.1000e+01, 2.2000e+01, 2.3000e+01, 2.4000e+01, 2.5500e+01, 2.7000e+01, 
 2.8000e+01, 3.0000e+01, 3.2000e+01, 3.4000e+01, 3.6000e+01, 3.8000e+01, 4.0000e+01, 4.2500e+01, 4.5000e+01, 4.7500e+01, 
 5.0000e+01, 5.2500e+01, 5.5000e+01, 5.7500e+01, 6.0000e+01, 6.3000e+01, 6.6000e+01, 6.9000e+01, 7.2000e+01, 7.6000e+01, 
 8.0000e+01, 8.4000e+01, 8.8000e+01, 9.2000e+01, 9.6000e+01, 1.0000e+02, 1.0500e+02, 1.1000e+02, 1.1500e+02, 1.2000e+02, 
 1.2750e+02, 1.3500e+02, 1.4250e+02, 1.5000e+02, 1.6000e+02, 1.7000e+02, 1.8000e+02, 1.9000e+02, 2.0000e+02, 2.1000e+02, 
 2.2000e+02, 2.3000e+02, 2.4000e+02, 2.5500e+02, 2.7000e+02, 2.8000e+02, 3.0000e+02, 3.2000e+02, 3.4000e+02, 3.6000e+02, 
 3.8000e+02, 4.0000e+02, 4.2500e+02, 4.5000e+02, 4.7500e+02, 5.0000e+02, 5.2500e+02, 5.5000e+02, 5.7500e+02, 6.0000e+02, 
 6.3000e+02, 6.6000e+02, 6.9000e+02, 7.2000e+02, 7.6000e+02, 8.0000e+02, 8.4000e+02, 8.8000e+02, 9.2000e+02, 9.6000e+02, 
 1.0000e+03, 1.0500e+03, 1.1000e+03, 1.1500e+03, 1.2000e+03, 1.2750e+03, 1.3500e+03, 1.4250e+03, 1.5000e+03, 1.6000e+03, 
 1.7000e+03, 1.8000e+03, 1.9000e+03, 2.0000e+03, 2.1000e+03, 2.2000e+03, 2.3000e+03, 2.4000e+03, 2.5500e+03, 2.7000e+03, 
 2.8000e+03, 3.0000e+03, 3.2000e+03, 3.4000e+03, 3.6000e+03, 3.8000e+03, 4.0000e+03, 4.2500e+03, 4.5000e+03, 4.7500e+03, 
 5.0000e+03, 5.2500e+03, 5.5000e+03, 5.7500e+03, 6.0000e+03, 6.3000e+03, 6.6000e+03, 6.9000e+03, 7.2000e+03, 7.6000e+03, 
 8.0000e+03, 8.4000e+03, 8.8000e+03, 9.2000e+03, 9.6000e+03, 1.0000e+04, 1.0500e+04, 1.1000e+04, 1.1500e+04, 1.2000e+04, 
 1.2750e+04, 1.3500e+04, 1.4250e+04, 1.5000e+04, 1.6000e+04, 1.7000e+04, 1.8000e+04, 1.9000e+04, 2.0000e+04, 2.1000e+04, 
 2.2000e+04, 2.3000e+04, 2.4000e+04, 2.5500e+04, 2.7000e+04, 2.8000e+04, 3.0000e+04, 3.2000e+04, 3.4000e+04, 3.6000e+04, 
 3.8000e+04, 4.0000e+04, 4.2500e+04, 4.5000e+04, 4.7500e+04, 5.0000e+04, 5.2500e+04, 5.5000e+04, 5.7500e+04, 6.0000e+04, 
 6.3000e+04, 6.6000e+04, 6.9000e+04, 7.2000e+04, 7.6000e+04, 8.0000e+04, 8.4000e+04, 8.8000e+04, 9.2000e+04, 9.6000e+04, 
 1.0000e+05, 1.0500e+05, 1.1000e+05, 1.1500e+05, 1.2000e+05, 1.2750e+05, 1.3500e+05, 1.4250e+05, 1.5000e+05, 1.6000e+05, 
 1.7000e+05, 1.8000e+05, 1.9000e+05, 2.0000e+05, 2.1000e+05, 2.2000e+05, 2.3000e+05, 2.4000e+05, 2.5500e+05, 2.7000e+05, 
 2.8000e+05, 3.0000e+05, 3.2000e+05, 3.4000e+05, 3.6000e+05, 3.8000e+05, 4.0000e+05, 4.2500e+05, 4.5000e+05, 4.7500e+05, 
 5.0000e+05, 5.2500e+05, 5.5000e+05, 5.7500e+05, 6.0000e+05, 6.3000e+05, 6.6000e+05, 6.9000e+05, 7.2000e+05, 7.6000e+05, 
 8.0000e+05, 8.4000e+05, 8.8000e+05, 9.2000e+05, 9.6000e+05, 1.0000e+06, 1.1000e+06, 1.2000e+06, 1.3000e+06, 1.4000e+06, 
 1.5000e+06, 1.6000e+06, 1.7000e+06, 1.8000e+06, 1.9000e+06, 2.0000e+06, 2.1000e+06, 2.2000e+06, 2.3000e+06, 2.4000e+06, 
 2.5000e+06, 2.6000e+06, 2.7000e+06, 2.8000e+06, 2.9000e+06, 3.0000e+06, 3.1000e+06, 3.2000e+06, 3.3000e+06, 3.4000e+06, 
 3.5000e+06, 3.6000e+06, 3.7000e+06, 3.8000e+06, 3.9000e+06, 4.0000e+06, 4.1000e+06, 4.2000e+06, 4.3000e+06, 4.4000e+06, 
 4.5000e+06, 4.6000e+06, 4.7000e+06, 4.8000e+06, 4.9000e+06, 5.0000e+06, 5.1000e+06, 5.2000e+06, 5.3000e+06, 5.4000e+06, 
 5.5000e+06, 5.6000e+06, 5.7000e+06, 5.8000e+06, 5.9000e+06, 6.0000e+06, 6.1000e+06, 6.2000e+06, 6.3000e+06, 6.4000e+06, 
 6.5000e+06, 6.6000e+06, 6.7000e+06, 6.8000e+06, 6.9000e+06, 7.0000e+06, 7.1000e+06, 7.2000e+06, 7.3000e+06, 7.4000e+06, 
 7.5000e+06, 7.6000e+06, 7.7000e+06, 7.8000e+06, 7.9000e+06, 8.0000e+06, 8.1000e+06, 8.2000e+06, 8.3000e+06, 8.4000e+06, 
 8.5000e+06, 8.6000e+06, 8.7000e+06, 8.8000e+06, 8.9000e+06, 9.0000e+06, 9.1000e+06, 9.2000e+06, 9.3000e+06, 9.4000e+06, 
 9.5000e+06, 9.6000e+06, 9.7000e+06, 9.8000e+06, 9.9000e+06, 1.0000e+07, 1.0100e+07, 1.0200e+07, 1.0300e+07, 1.0400e+07, 
 1.0500e+07, 1.0600e+07, 1.0700e+07, 1.0800e+07, 1.0900e+07, 1.1000e+07, 1.1100e+07, 1.1200e+07, 1.1300e+07, 1.1400e+07, 
 1.1500e+07, 1.1600e+07, 1.1700e+07, 1.1800e+07, 1.1900e+07, 1.2000e+07, 1.2100e+07, 1.2200e+07, 1.2300e+07, 1.2400e+07, 
 1.2500e+07, 1.2600e+07, 1.2700e+07, 1.2800e+07, 1.2900e+07, 1.3000e+07, 1.3100e+07, 1.3200e+07, 1.3300e+07, 1.3400e+07, 
 1.3500e+07, 1.3600e+07, 1.3700e+07, 1.3800e+07, 1.3900e+07, 1.4000e+07, 1.4100e+07, 1.4200e+07, 1.4300e+07, 1.4400e+07, 
 1.4500e+07, 1.4600e+07, 1.4700e+07, 1.4800e+07, 1.4900e+07, 1.5000e+07, 1.5100e+07, 1.5200e+07, 1.5300e+07, 1.5400e+07, 
 1.5500e+07, 1.5600e+07, 1.5700e+07, 1.5800e+07, 1.5900e+07, 1.6000e+07, 1.6100e+07, 1.6200e+07, 1.6300e+07, 1.6400e+07, 
 1.6500e+07, 1.6600e+07, 1.6700e+07, 1.6800e+07, 1.6900e+07, 1.7000e+07, 1.7100e+07, 1.7200e+07, 1.7300e+07, 1.7400e+07, 
 1.7500e+07, 1.7600e+07, 1.7700e+07, 1.7800e+07, 1.7900e+07, 1.8000e+07, 1.8100e+07, 1.8200e+07, 1.8300e+07, 1.8400e+07, 
 1.8500e+07, 1.8600e+07, 1.8700e+07, 1.8800e+07, 1.8900e+07, 1.9000e+07, 1.9100e+07, 1.9200e+07, 1.9300e+07, 1.9400e+07, 
 1.9500e+07, 1.9600e+07, 1.9700e+07, 1.9800e+07, 1.9900e+07, 2.0000e+07, 2.0500e+07, 2.1000e+07, 2.1500e+07, 2.2000e+07, 
 2.2500e+07, 2.3000e+07, 2.3500e+07, 2.4000e+07, 2.4500e+07, 2.5000e+07, 2.5500e+07, 2.6000e+07, 2.6500e+07, 2.7000e+07, 
 2.7500e+07, 2.8000e+07, 2.8500e+07, 2.9000e+07, 2.9500e+07, 3.0000e+07, 3.1000e+07, 3.2000e+07, 3.3000e+07, 3.4000e+07, 
 3.5000e+07, 3.6000e+07, 3.7000e+07, 3.8000e+07, 3.9000e+07, 4.0000e+07, 4.2000e+07, 4.4000e+07, 4.6000e+07, 4.8000e+07, 
 5.0000e+07, 5.2000e+07, 5.4000e+07, 5.6000e+07, 5.8000e+07, 6.0000e+07};

/*** LANL 618 */
static const string grpStructureName4 = "LANL 618";
static const int grpEnergyPoint4 = 619;
static double grpEnergyGrid4[grpEnergyPoint4] = {
 1.00000000000000e-05, 2.56901129797510e-05, 4.23558357164050e-05, 6.98329672839171e-05, 1.15135098557100e-04, 1.39000000000000e-04, 1.89825685995247e-04, 2.43741005558083e-04, 3.12969646225607e-04, 4.01860980405450e-04,
 5.15999712815652e-04, 6.62556746258873e-04, 8.50739702194323e-04, 1.09237140060287e-03, 1.23781896276759e-03, 1.40263264283687e-03, 1.58939100945164e-03, 1.80101596367844e-03, 2.04081845319089e-03, 2.31255027322349e-03,
 2.62046276474246e-03, 2.96937332818714e-03, 3.36474079341315e-03, 3.81275082502696e-03, 4.32041269930856e-03, 4.89566896683177e-03, 5.54751971649269e-03, 6.28616338510141e-03, 7.12315631555298e-03, 8.07159355992206e-03,
 9.14631375620984e-03, 1.03641312841130e-02, 1.10325603234355e-02, 1.17440993319742e-02, 1.25015286638674e-02, 1.33078079906897e-02, 1.41660878664320e-02, 1.50797220383603e-02, 1.60522805518561e-02, 1.70875637004458e-02,
 1.81896168755305e-02, 1.93627463738410e-02, 2.06115362243856e-02, 2.19408661006432e-02, 2.33559303879934e-02, 2.48622584808902e-02, 2.64657363890912e-02, 2.81726297373683e-02, 2.99896082485731e-02, 3.19237718057234e-02,
 3.39826781949507e-02, 3.61743726377138e-02, 3.85074192276762e-02, 4.09909343950883e-02, 4.36346225294370e-02, 4.64488138995581e-02, 4.94445050193864e-02, 5.26334016170732e-02, 5.60279643753727e-02, 5.96414576220314e-02,
 6.34880011604368e-02, 6.75826254430556e-02, 7.19413303032538e-02, 7.65811474749932e-02, 8.15202071447017e-02, 8.67778087953710e-02, 9.23744966197059e-02, 9.83321397970035e-02, 1.04674017947447e-01, 1.11424912097725e-01,
 1.18611201513438e-01, 1.26260966776645e-01, 1.34404099511350e-01, 1.43072419185677e-01, 1.52000000000000e-01, 1.62122290476778e-01, 1.72578279879602e-01, 1.83708622661412e-01, 1.95556810878505e-01, 2.08169141583816e-01,
 2.21594897733660e-01, 2.35886540761951e-01, 2.51099915574398e-01, 2.67294468763689e-01, 2.84533480898340e-01, 3.02884313792894e-01, 3.22418673725673e-01, 3.43212891632624e-01, 3.65348221372105e-01, 3.88911157226101e-01,
 4.14000000000000e-01, 4.40694076191185e-01, 4.69116402183442e-01, 4.99371810711756e-01, 5.31578525442442e-01, 5.65862394813206e-01, 6.02357383788648e-01, 6.41206097331274e-01, 6.82560337633487e-01, 7.26581697287950e-01,
 7.73442190714156e-01, 8.23324926308510e-01, 8.76424821944364e-01, 9.32949366617847e-01, 9.93119431215624e-01, 1.05717013157269e+00, 1.13000000000000e+00, 1.19793069922005e+00, 1.27519059148733e+00, 1.35743331870246e+00,
 1.44498024610924e+00, 1.53817346522906e+00, 1.63737713059081e+00, 1.74297888267274e+00, 1.85539136261598e+00, 1.97505382462877e+00, 2.10243385238185e+00, 2.23802918610180e+00, 2.38236966750182e+00, 2.53601931014967e+00,
 2.69957850336301e+00, 2.87368635824370e+00, 3.06000000000000e+00, 3.25631325144309e+00, 3.46632741266196e+00, 3.68988632357374e+00, 3.92786354548104e+00, 4.18118897955002e+00, 4.45085250041942e+00, 4.73790782415717e+00,
 5.04347662567888e+00, 5.36875292171691e+00, 5.71500773646672e+00, 6.08359406814152e+00, 6.47595217584221e+00, 6.89361520740109e+00, 7.33821519019035e+00, 7.81148940830449e+00, 8.32000000000000e+00, 8.85157713916813e+00,
 9.42245481732848e+00, 1.00301509424501e+01, 1.06770401003478e+01, 1.13656500244640e+01, 1.20986714730416e+01, 1.28789687433204e+01, 1.37095908638408e+01, 1.45937835085895e+01, 1.55350016795403e+01, 1.65369232071503e+01,
 1.76034631215617e+01, 1.87387889506673e+01, 1.99473370048166e+01, 2.12338297117944e+01, 2.26000000000000e+01, 2.40610812906042e+01, 2.56128877094204e+01, 2.72647770435633e+01, 2.90232040865040e+01, 3.08950399301257e+01,
 3.28875988136648e+01, 3.50086667042598e+01, 3.72665317207867e+01, 3.96700165198641e+01, 4.22285127705753e+01, 4.49520178526194e+01, 4.78511739212901e+01, 5.09373094919281e+01, 5.42224837063415e+01, 5.77195334541645e+01,
 6.14000000000000e+01, 6.54048000453254e+01, 6.96230472348795e+01, 7.41133479945056e+01, 7.88932482720022e+01, 8.39814256315774e+01, 8.93977622368364e+01, 9.51634225407686e+01, 1.01300935986307e+02, 1.07834285040617e+02,
 1.14788998907105e+02, 1.22192253281342e+02, 1.30072976540676e+02, 1.38461962782503e+02, 1.47391992152865e+02, 1.56897958935589e+02, 1.67000000000000e+02, 1.77788679457205e+02, 1.89255064140519e+02, 2.01460967099726e+02,
 2.14454083165892e+02, 2.28285183222401e+02, 2.43008312593295e+02, 2.58681002226541e+02, 2.75364493497472e+02, 2.93123977510781e+02, 3.12028849836190e+02, 3.32152981673137e+02, 3.53575008504100e+02, 3.76378637364449e+02,
 4.00652973929511e+02, 4.26492870696926e+02, 4.54000000000000e+02, 4.83279736674251e+02, 5.14448601797023e+02, 5.47627686010971e+02, 5.82946637308688e+02, 6.20543465259898e+02, 6.60565080286848e+02, 7.03167867719981e+02,
 7.48518298877006e+02, 7.96793581553195e+02, 8.48182352464692e+02, 9.02885414350579e+02, 9.61116520613947e+02, 1.02310321056796e+03, 1.05557999269466e+03, 1.08908769855066e+03, 1.12365905316802e+03, 1.15932782038279e+03,
 1.19612883581024e+03, 1.23500000000000e+03, 1.27327251787187e+03, 1.31369052626409e+03, 1.35539153996700e+03, 1.39841628594101e+03, 1.44280678395902e+03, 1.48860638764470e+03, 1.53585982681347e+03, 1.58461325115751e+03,
 1.63491427531748e+03, 1.68681202538499e+03, 1.74035718688117e+03, 1.79560205425833e+03, 1.85260058197288e+03, 1.91140843717952e+03, 1.97208305409813e+03, 2.03468369010644e+03, 2.09927148361327e+03, 2.16590951376885e+03,
 2.23466286207060e+03, 2.30559867592442e+03, 2.37878623422368e+03, 2.45429701500989e+03, 2.53220476528119e+03, 2.61258557301668e+03, 2.69551794148722e+03, 2.78108286592499e+03, 2.86936391262682e+03, 2.96044730056855e+03,
 3.05442198561012e+03, 3.15137974737356e+03, 3.25141527887886e+03, 3.35000000000000e+03, 3.46111354800741e+03, 3.57098108576248e+03, 3.68433619353942e+03, 3.80128957869464e+03, 3.92195546281323e+03, 4.04645169326264e+03,
 4.17489985828732e+03, 4.30742540575688e+03, 4.44415776568380e+03, 4.58523047663021e+03, 4.73078131612718e+03, 4.88095243523415e+03, 5.03589049736953e+03, 5.19574682154838e+03, 5.36067753016696e+03, 5.53084370147834e+03,
 5.70641152690821e+03, 5.88755247336443e+03, 6.07444345069879e+03, 6.26726698448458e+03, 6.46621139427874e+03, 6.67147097754267e+03, 6.88324619940125e+03, 7.10174388842549e+03, 7.32717743863004e+03, 7.55976701788271e+03,
 7.79973978292963e+03, 8.04733010124613e+03, 8.30277977992978e+03, 8.56633830185941e+03, 8.83826306935050e+03, 9.12000000000000e+03, 9.40828206378196e+03, 9.70693299519909e+03, 1.00150641248322e+04, 1.03329763864764e+04,
 1.06609802665909e+04, 1.09993961075332e+04, 1.13485544204187e+04, 1.17087962079117e+04, 1.20804732972634e+04, 1.24639486839205e+04, 1.28595968860421e+04, 1.32678043102699e+04, 1.36889696291092e+04, 1.41235041702888e+04,
 1.45718323184816e+04, 1.50343919297757e+04, 1.55116347593038e+04, 1.60040269024456e+04, 1.65120492500366e+04, 1.70361979580257e+04, 1.75769849320427e+04, 1.81349383273462e+04, 1.87106030646422e+04, 1.93045413622771e+04,
 1.99173332853231e+04, 2.05495773120946e+04, 2.12018909186467e+04, 2.18749111818289e+04, 2.25692954014803e+04, 2.32857217423771e+04, 2.40248898965561e+04, 2.48000000000000e+04, 2.55743621709957e+04, 2.63861795709192e+04,
 2.72237668213834e+04, 2.80879419452551e+04, 2.89795489322345e+04, 2.98994585631306e+04, 3.08485692603026e+04, 3.18278079650967e+04, 3.28381310431359e+04, 3.38805252183471e+04, 3.49560085366367e+04, 3.60656313601573e+04,
 3.72104773931352e+04, 3.83916647402616e+04, 3.96103469986807e+04, 4.08677143846407e+04, 4.21649948959093e+04, 4.35034555110877e+04, 4.48844034269952e+04, 4.63091873353325e+04, 4.77791987398702e+04, 4.92958733154505e+04,
 5.08606923101270e+04, 5.24751839918138e+04, 5.41409251408564e+04, 5.58595425899810e+04, 5.76327148131282e+04, 5.94621735647209e+04, 6.13497055709682e+04, 6.32971542748575e+04, 6.53064216365378e+04, 6.76000000000000e+04,
 6.95183239638479e+04, 7.17250724500870e+04, 7.40018706527728e+04, 7.63509421885996e+04, 7.87745812594328e+04, 8.12751548929221e+04, 8.38551052542408e+04, 8.65169520312063e+04, 8.92632948951132e+04, 9.20968160396814e+04,
 9.50202828005989e+04, 9.80365503582183e+04, 1.01148564526046e+05, 1.04359364627745e+05, 1.07672086465471e+05, 1.11089965382423e+05, 1.14616339422619e+05, 1.18254652590966e+05, 1.22008458216826e+05, 1.25881422424340e+05,
 1.29877327712922e+05, 1.34000076651408e+05, 1.38253695689465e+05, 1.42642339089993e+05, 1.47170292986351e+05, 1.51841979568379e+05, 1.56661961401289e+05, 1.61634945881659e+05, 1.66765789834876e+05, 1.72059504258514e+05,
 1.77521259216285e+05, 1.84000000000000e+05, 1.88970396775857e+05, 1.94968961085980e+05, 2.01157940267409e+05, 2.07543378736997e+05, 2.14131512781987e+05, 2.20928776650624e+05, 2.27941808836123e+05, 2.35177458560091e+05,
 2.42642792461767e+05, 2.50345101499601e+05, 2.58291908071908e+05, 2.66490973363555e+05, 2.74950304925867e+05, 2.83678164497131e+05, 2.92683076071361e+05, 3.03000000000000e+05, 3.11559512696998e+05, 3.21449473268761e+05,
 3.31653374889103e+05, 3.42181183116660e+05, 3.53043179850858e+05, 3.64249973373642e+05, 3.75812508709979e+05, 3.87742078317220e+05, 4.00050333113793e+05, 4.12749293857975e+05, 4.25851362887876e+05, 4.39369336234074e+05,
 4.53316416116763e+05, 4.67706223839590e+05, 4.82552813092796e+05, 5.00000000000000e+05, 5.13674795672507e+05, 5.29980584033558e+05, 5.46803973679148e+05, 5.64161395037774e+05, 5.82069800095720e+05, 6.00546678953079e+05,
 6.19610076905320e+05, 6.39278612067076e+05, 6.59571493555382e+05, 6.80508540250102e+05, 7.02110200149880e+05, 7.24397570342515e+05, 7.47392417609257e+05, 7.71117199683167e+05, 7.95595087182277e+05, 8.23000000000000e+05,
 8.46906561847805e+05, 8.73790261954204e+05, 9.01527342308164e+05, 9.30144892106635e+05, 9.59670860449985e+05, 9.90134083638263e+05, 1.02156431333394e+06, 1.05399224561864e+06, 1.08744955097221e+06, 1.12196890520344e+06,
 1.15758402136263e+06, 1.19432968266720e+06, 1.23224177647237e+06, 1.27135732932036e+06, 1.31171454310194e+06, 1.35300000000000e+06, 1.39631286281399e+06, 1.44063659101453e+06, 1.48636730538123e+06, 1.53354966844928e+06,
 1.58222976049498e+06, 1.63245512453958e+06, 1.68427481278184e+06, 1.73800000000000e+06, 1.79290120550119e+06, 1.84981399907304e+06, 1.90853339864316e+06, 1.96911675204194e+06, 2.03162322751532e+06, 2.09611387151098e+06,
 2.16265166829887e+06, 2.23200000000000e+06, 2.30213071747361e+06, 2.37520819095458e+06, 2.45060539245526e+06, 2.52839595804746e+06, 2.60865586126285e+06, 2.69146348729184e+06, 2.77689970953790e+06, 2.86500000000000e+06,
 2.95599435377371e+06, 3.04982768711059e+06, 3.14663961018459e+06, 3.24652467358350e+06, 3.34958042925295e+06, 3.45590752576975e+06, 3.56560980663947e+06, 3.68000000000000e+06, 3.79557188183090e+06, 3.91605626676799e+06,
 4.04036523663342e+06, 4.10399173096370e+06, 4.16862019678508e+06, 4.23426641285263e+06, 4.30094640640062e+06, 4.36867645705557e+06, 4.43747310081080e+06, 4.50735313406362e+06, 4.57833361771614e+06, 4.65043188134056e+06,
 4.72366552741015e+06, 4.79805243559677e+06, 4.87361076713619e+06, 4.95035896926199e+06, 5.02831577970941e+06, 5.10750023129011e+06, 5.18793165653889e+06, 5.26962969243371e+06, 5.35261428518990e+06, 5.43690569513000e+06,
 5.52252450163020e+06, 5.60949160814471e+06, 5.69782824730923e+06, 5.78755598612484e+06, 5.87869673122347e+06, 5.97127273421627e+06, 6.07000000000000e+06, 6.16082127790678e+06, 6.25784009604591e+06, 6.35638673826052e+06,
 6.45648526427892e+06, 6.55816011271502e+06, 6.66143610703488e+06, 6.76633846161729e+06, 6.87289278790972e+06, 6.98112510068126e+06, 7.09106182437398e+06, 7.20272979955440e+06, 7.31615628946642e+06, 7.43136898668758e+06,
 7.54839601989007e+06, 7.66726596070820e+06, 7.79000000000000e+06, 7.91065110850296e+06, 8.03522573689061e+06, 8.16176213022340e+06, 8.29029118180400e+06, 8.42084427143382e+06, 8.55345327307423e+06, 8.68815056262843e+06,
 8.82496902584595e+06, 8.96394206635151e+06, 9.10510361380034e+06, 9.24848813216205e+06, 9.39413062813476e+06, 9.54206665969188e+06, 9.69233234476344e+06, 9.84496437005408e+06, 1.00000000000000e+07, 1.01250000000000e+07,
 1.02500000000000e+07, 1.03750000000000e+07, 1.05000000000000e+07, 1.06250000000000e+07, 1.07500000000000e+07, 1.08750000000000e+07, 1.10000000000000e+07, 1.11250000000000e+07, 1.12500000000000e+07, 1.13750000000000e+07,
 1.15000000000000e+07, 1.16250000000000e+07, 1.17500000000000e+07, 1.18750000000000e+07, 1.20000000000000e+07, 1.21250000000000e+07, 1.22500000000000e+07, 1.23750000000000e+07, 1.25000000000000e+07, 1.26250000000000e+07,
 1.27500000000000e+07, 1.28750000000000e+07, 1.30000000000000e+07, 1.31250000000000e+07, 1.32500000000000e+07, 1.33750000000000e+07, 1.35000000000000e+07, 1.36250000000000e+07, 1.37500000000000e+07, 1.38750000000000e+07,
 1.40000000000000e+07, 1.41250000000000e+07, 1.42500000000000e+07, 1.43750000000000e+07, 1.45000000000000e+07, 1.46250000000000e+07, 1.47500000000000e+07, 1.48750000000000e+07, 1.50000000000000e+07, 1.51250000000000e+07,
 1.52500000000000e+07, 1.53750000000000e+07, 1.55000000000000e+07, 1.56250000000000e+07, 1.57500000000000e+07, 1.58750000000000e+07, 1.60000000000000e+07, 1.61250000000000e+07, 1.62500000000000e+07, 1.63750000000000e+07,
 1.65000000000000e+07, 1.66250000000000e+07, 1.67500000000000e+07, 1.68750000000000e+07, 1.70000000000000e+07, 1.71250000000000e+07, 1.72500000000000e+07, 1.73750000000000e+07, 1.75000000000000e+07, 1.76250000000000e+07,
 1.77500000000000e+07, 1.78750000000000e+07, 1.80000000000000e+07, 1.81250000000000e+07, 1.82500000000000e+07, 1.83750000000000e+07, 1.85000000000000e+07, 1.86250000000000e+07, 1.87500000000000e+07, 1.88750000000000e+07,
 1.90000000000000e+07, 1.91250000000000e+07, 1.92500000000000e+07, 1.93750000000000e+07, 1.95000000000000e+07, 1.96250000000000e+07, 1.97500000000000e+07, 1.98750000000000e+07, 2.00000000000000e+07};
