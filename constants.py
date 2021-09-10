### define constants
#nGores = 280;
#nTendons = 280;
#
#gravity            =9.8066 ; # m/sec^2
#Newtons2PoundForce =0.224809;
#
#LinDenT =0.015205; # kg/m	# linear density of the tendon
#E1=212.5e+06;
#E2=246.6e+06;
#SuspendedPayload=2273 ; # kg
#NadirWeight = SuspendedPayload*gravity ;
#tauVol = 0.0125  ;
#
#BalloonMass     =2331 ; # kg
#GrossMass  = 4606;
##TopMass = 91.4*0  ;	# top plate mass in kg
#TopMass = 91.4  ;	# top plate mass in kg
#BotMass = 64.11 ;	# bottom plate mass in kg
#ApexWeight = gravity * TopMass ;
#NadirWeight = SuspendedPayload*gravity + BotMass*gravity;
#buoyancySEA_LEVEL = 1.0555; # kg/m^3 ; 
#
### this is b
##buoyancy_Float = 0.084887 *1.125;  # N/m^3
#
## debug
#buoyancy_Float = -0.0849;
#
#b_d = buoyancy_Float;
#buoyancy_spool=  0.084887/tauVol  ; # P0=-170
## buoyancy_spool=  buoyancy_Float * 10 # for 1# vol
#buoyancy = buoyancy_spool;
#ethickness = 38.1e-06;
#BalloonFilmWeightArealDensity = ethickness * 920 * gravity;	# weight of the film per unit area
#BalloonFilmMassArealDensity = ethickness * 920 ;
#CapLength = 52.9 ;
#
#
#ETape = 761560;
## spring constant for the film-net
##ETape_net = ETape/400;
#ETape_net = ETape/4;
#
#RBottomPlate=0.10;
#RTopPlate   =0.35;
#ny=87;
#
#Spool_Pntr=62; 
#Collar_Pntr=58;
#Top_Pntr = ny;
#nStations = ny;
#
### P_0
##P_0 = 100; # Pascal
#P_0 = -100; # Pascal

#######################################################################

gravity            =9.8066 ; # m/sec^2
Newtons2PoundForce =0.224809;


buoyancySEA_LEVEL = 1.0555; # kg/m^3 ; 

PressureAtNadir = -100; # negative means inside gas pushes outward
PressureAtNadir = -2.8208   # KSTEP = 5
V_d =  523152  ;         # design volume


TargetVolume = V_d ;
TargetVolume =  4.2572e+05  # for KSTEP = 5
tauVol = V_d/TargetVolume ;      #  0.0125  ;
b_d = 0.084887   # N/m^3
buoyancy = b_d * tauVol  ;
buoyancy_spool=  b_d/ (0.0125)   ; # P0=170



LinDenT =0.015205;
TendonMassDensity = 0.015205 ; # kg/m
    wTape         =  0.015205  *gravity
E1=212.5e+06;
E2=246.6e+06;
YoungsModulus = mean([E1, E2]) ;
PoissonsRatio = 0.5 * (0.58 + 0.673);
SuspendedPayload=2273 ; # kg




BalloonMass     =2331 ; # kg
GrossMass  = 4606;
TopMass = 91.4 ;
BotMass = 64.1

ApexWeight = gravity * TopMass ;
NadirWeight = SuspendedPayload*gravity + BotMass*gravity;


ethickness = 38.1e-06;
BalloonFilmWeightArealDensity = ethickness * 920 * gravity;
BalloonFilmMassArealDensity = ethickness * 920 ;
CapLength = 52.9 *0  ;


ETape = 761560 ;
ETape_net = ETape/4;	#filmnet spring constant
#ETape_net = ETape*100;	#filmnet spring constant
RBottomPlate=0.10;
RTopPlate   =0.35;
ny=87

Spool_Pntr=62; 
Collar_Pntr=58;
Top_Pntr = ny
nStations = ny;

P_0 = PressureAtNadir;
nGores = 280;
buoyancy_Float = -buoyancy;
