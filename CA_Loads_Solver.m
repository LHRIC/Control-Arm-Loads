clc, clearvars
%Assumptions: 
%Driveshaft inboard braking reduction of moment in the rear
%Isn’t accounted for, therefore the actual factor of safety on the rear
%under braking should be higher than calculated
%Force is applied to the inboard side of the hardpoint

%% Inputs
side = ("front"); %"front" or "rear"
WB = 1.535; %wheelbase(m)
Tf = 1.241; %front track(m)
Tr = 1.241; %rear track(m)
lateralG = 1.7; %max lateral G's
maxAccel = 1.5; %max forward acceleration G's
maxBrake = 1.5; %max braking G's
maxBump = 2.5; %Max G's in bump, static = 0
Fy = 0; %Max lateral grip (N) % 2471 front
Fx = 0; %Max longitudinal grip(N) %2330 front
LLTD = .5; %Lateral load transfer distribution towards rear
maxVel = 27.5; %max velocity (m/s)
mass = 295; %car + driver mass(kg)
WD = .55; %Static weight distribution towards rear
cgHeight = .32; %CG height(m)
CoP = .6; %center of pressure distribution to rear
CLA = 3.7; %cl*A
rho = 1.2754; %air density
mechanicalTrail = .0127; %mechanical trail
scrubRadius = .00635; %scrub radius(m) assuming negative scrub
maxTireMz = 60; %Tire self aligning moment (Nm)
Mx = 50; %overturning moment??
My = 80; %rolling resistance


%% Inputting Hardpoint Coordinate
if side == "rear"
    xlRange = 'B2:D12'; %input front or rear spreadsheet cells
elseif side == "Front"
     xlRange = 'B2:D12';
else
    xlRange = 'B19:D29';
end

hp_geom = readmatrix('hardpoint_geometry_test.xls', 'range', xlRange);
%Assume you're looking at right corner from front
%Positive x-axis is forwards
%Positive y-axis is rightwards(outboard)
%Positive z-axis is upwards
UAA_f = hp_geom(4,1:3)/1000; %inboard upper fore x,y,z
UAA_a = hp_geom(5,1:3)/1000; %inboard upper aft x,y,z
UAA_bj = hp_geom(6,1:3)/1000; %outboard balljoint x,y,z
LAA_f = hp_geom(1,1:3)/1000;
LAA_a = hp_geom(2,1:3)/1000;
LAA_bj = hp_geom(3,1:3)/1000;
TIE_in = hp_geom(10,1:3)/1000;
TIE_out = hp_geom(9,1:3)/1000;
PR_in = hp_geom(8,1:3)/1000; %push/pullrod
PR_out = hp_geom(7,1:3)/1000;
wheelCenter = hp_geom(11, 1:3)/1000;


%% Generating Geometries
%Assume you're looking at right corner from front
%Start by assuming all members are in tension
%Assumes z vectors are point down, change if not
%Positive x-axis is backwards
%Positive y-axis is rightwards(outboard)
%Positive z-axis is upwards
%aft upper a arm geometry
UCAA_vector = (UAA_a - UAA_bj); %calculating x, y, and z vector lengths;
UCAA_Lt = (UCAA_vector(1)^2 + UCAA_vector(2)^2 + UCAA_vector(3)^2)^0.5; %total length
UCAA_unitX = (UCAA_vector(1)/UCAA_Lt);
UCAA_unitY = (UCAA_vector(2)/UCAA_Lt);
UCAA_unitZ = (UCAA_vector(3)/UCAA_Lt);
UCAA_Mx = -(UCAA_unitY)*(UAA_bj(3)-wheelCenter(3)) + (UCAA_unitZ)*(UAA_bj(2)-wheelCenter(2));
%Inwards vector is negative
UCAA_My = (UCAA_unitX)*(UAA_bj(3)-wheelCenter(3)) - (UCAA_unitZ)*(UAA_bj(1)-wheelCenter(1));
%Backwards vector is negative
UCAA_Mz = -(UCAA_unitX)*(UAA_bj(2)-wheelCenter(2)) + (UCAA_unitY)*(UAA_bj(1)-wheelCenter(1));
%Clockwise is negative
%fore upper a arm geometry
UCAF_vector = (UAA_f - UAA_bj);
UCAF_Lt = (UCAF_vector(1)^2 + UCAF_vector(2)^2 + UCAF_vector(3)^2)^0.5;
UCAF_unitX = (UCAF_vector(1)/UCAF_Lt);
UCAF_unitY = (UCAF_vector(2)/UCAF_Lt);
UCAF_unitZ = (UCAF_vector(3)/UCAF_Lt);
UCAF_Mx = -(UCAF_unitY)*(UAA_bj(3)-wheelCenter(3)) + (UCAF_unitZ)*(UAA_bj(2)-wheelCenter(2));
%Inwards vector is negative
UCAF_My = (UCAF_unitX)*(UAA_bj(3)-wheelCenter(3)) - (UCAF_unitZ)*(UAA_bj(1)-wheelCenter(1));
%Backwards vector is negative
UCAF_Mz = -(UCAF_unitX)*(UAA_bj(2)-wheelCenter(2)) + (UCAF_unitY)*(UAA_bj(1)-wheelCenter(1));
%Clockwise is negative
%aft lower a arm geometry
LCAA_vector = (LAA_a - LAA_bj); %calculating x, y, and z vector lengths;
LCAA_Lt = (LCAA_vector(1)^2 + LCAA_vector(2)^2 + LCAA_vector(3)^2)^0.5; %total length
LCAA_unitX = (LCAA_vector(1)/LCAA_Lt);
LCAA_unitY = (LCAA_vector(2)/LCAA_Lt);
LCAA_unitZ = (LCAA_vector(3)/LCAA_Lt);
LCAA_Mx = -(LCAA_unitY)*(LAA_bj(3)-wheelCenter(3)) + (LCAA_unitZ)*(LAA_bj(2)-wheelCenter(2));
LCAA_My = (LCAA_unitX)*(LAA_bj(3)-wheelCenter(3)) - (LCAA_unitZ)*(LAA_bj(1)-wheelCenter(1));
LCAA_Mz = -(LCAA_unitX)*(LAA_bj(2)-wheelCenter(2)) + (LCAA_unitY)*(LAA_bj(1)-wheelCenter(1));
%fore lower a arm geometry
LCAF_vector = (LAA_f - LAA_bj);
LCAF_Lt = (LCAF_vector(1)^2 + LCAF_vector(2)^2 + LCAF_vector(3)^2)^0.5;
LCAF_unitX = (LCAF_vector(1)/LCAF_Lt);
LCAF_unitY = (LCAF_vector(2)/LCAF_Lt);
LCAF_unitZ = (LCAF_vector(3)/LCAF_Lt);
LCAF_Mx = -(LCAF_unitY)*(LAA_bj(3)-wheelCenter(3)) + (LCAF_unitZ)*(LAA_bj(2)-wheelCenter(2));
LCAF_My = (LCAF_unitX)*(LAA_bj(3)-wheelCenter(3)) - (LCAF_unitZ)*(LAA_bj(1)-wheelCenter(1));
LCAF_Mz = -(LCAF_unitX)*(LAA_bj(2)-wheelCenter(2)) + (LCAF_unitY)*(LAA_bj(1)-wheelCenter(1));
%Tie rod geometry geometry
TIE_vector = (TIE_in - TIE_out);
TIE_Lt = (TIE_vector(1)^2 + TIE_vector(2)^2 + TIE_vector(3)^2)^0.5;
TIE_unitX = (TIE_vector(1)/TIE_Lt);
TIE_unitY = (TIE_vector(2)/TIE_Lt);
TIE_unitZ = (TIE_vector(3)/TIE_Lt);
TIE_Mx = -(TIE_unitY)*(TIE_out(3)-wheelCenter(3)) + (TIE_unitZ)*(TIE_out(2) - wheelCenter(2));
TIE_My = (TIE_unitX)*(TIE_out(3)-wheelCenter(3)) - (TIE_unitZ)*(TIE_out(1)-wheelCenter(1));
TIE_Mz = -(TIE_unitX)*(TIE_out(2)-wheelCenter(2)) + (TIE_unitY)*(TIE_out(1)-wheelCenter(1));
%Push/pullrod geometry
PR_vector = (PR_in - PR_out);
PR_Lt = (PR_vector(1)^2 + PR_vector(2)^2 + PR_vector(3)^2)^0.5;
PR_unitX = (PR_vector(1)/PR_Lt);
PR_unitY = (PR_vector(2)/PR_Lt);
PR_unitZ = (PR_vector(3)/PR_Lt);
PR_Mx = -(PR_unitY)*(PR_out(3)-wheelCenter(3)) + (PR_unitZ)*(PR_out(2)-wheelCenter(2));
PR_My = (PR_unitX)*(PR_out(3)-wheelCenter(3)) -(PR_unitZ)*(PR_out(1)-wheelCenter(1));
PR_Mz = -(PR_unitX)*(PR_out(2)-wheelCenter(2)) +(PR_unitY)*(PR_out(1)-wheelCenter(1));
vectorTable = [UCAF_vector; UCAA_vector; LCAF_vector;...
LCAA_vector; TIE_vector; PR_vector];


%% Vehicle State Calculations
%Change weight distribution based on front or rear
if side == "Front"
T = Tf;
WD = 1-WD;
LLTD = 1-LLTD;
CoP = 1-CoP;
Ax = maxBrake;
elseif side == "front"
T = Tf;
WD = 1-WD;
LLTD = 1-LLTD;
CoP = 1-CoP;
Ax = maxBrake;
else
T = Tr;
WD = WD;
LLTD = LLTD;
CoP = CoP;
Ax = maxAccel;
end
weight = mass*9.81*WD/2; %Car+driver weight(N) per tire
downforce = (.5*rho*CLA*CoP*maxVel^2)/2; %Aero downforce per tire(N)
latLoadTransfer = ((mass*9.81*cgHeight/T)*lateralG*LLTD);%(N)
bumpLoad = mass*WD*9.81*maxBump/2; %I'm not sure if this is right
accelWeightTransfer = abs((cgHeight/WB)*(mass*9.81*Ax)/2);
Fz = weight + latLoadTransfer + downforce + bumpLoad + accelWeightTransfer;
steeringMz = Fy * mechanicalTrail; %steering contribution to Mz
scrubMz = -scrubRadius * Fx;
Mz = steeringMz + scrubMz + maxTireMz;
%Assume tie rod takes entire Mz
MzTIE = Mz/ (TIE_unitY * (TIE_out(1) - wheelCenter(1)));


%% Huge Matrix Time!!!
%Ax = B
%x = inv(A) * B
%input all unit vector and moments for each control arm
A = [UCAF_unitX, UCAA_unitX, LCAF_unitX, LCAA_unitX, TIE_unitX, PR_unitX;...
UCAF_unitY, UCAA_unitY, LCAF_unitY, LCAA_unitY, TIE_unitY, PR_unitY;...
UCAF_unitZ, UCAA_unitZ, LCAF_unitZ, LCAA_unitZ, TIE_unitZ, PR_unitZ;...
UCAF_Mx, UCAA_Mx, LCAF_Mx, LCAA_Mx, TIE_Mx, PR_Mx;...
UCAF_My, UCAA_My, LCAF_My, LCAA_My, TIE_My, PR_My;...
UCAF_Mz, UCAA_Mz, LCAF_Mz, LCAA_Mz, TIE_Mz, PR_Mz];
%X = [FUCAF; FUCAA; F3LCAF; F4LCAA; F5TIE; F6PR]
%tire vectors and moments, subtract Tie rod forces
B = [Fx; Fy ; Fz ; 0; 0; 0]
force = A\B;

%Add Mz contribution to force to overall tie rod force
forceTIE =force(5) - MzTIE;
disp( "(+) = compression (-) = tension (N)")


if side == "Front"
    disp('Front Control Arms')
elseif side == "front"
    disp('Front Control Arms')
else
    disp('Rear Control Arms')
end

disp(" Upper fore CA: " + force(1))
disp(" Upper aft CA: " + force(2))
disp(" Lower fore CA: " + force(3))
disp(" Lower aft CA: " + force(4))
disp(" Tie rod: " + forceTIE)
disp(" Pull/pushrod: " + force(6))
