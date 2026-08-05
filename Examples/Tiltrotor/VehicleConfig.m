clear
clc
close all

%% things to add

% throttle forcing
% control surface saturation
% fuselage inertia contribution

% revamped aero models!
% verify correct tangential flow

% fix limited control inputs to aero model

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..'));  % reach the root
setupPaths();

%% Surface Definitions

% wings ===========================================================

wing = StripModel("wing_R");
wing.setPlanform('Span',0.6,'RootChord',0.15,'TaperRatio',0.7, ...
                 'TwistTip',0,'Incidence',0,'SweepQC',-2);
wing.setLocation('RootQC',[-0.05 ; 0.05; 0],'Dihedral',0,'SpanSide',+1);
wing.addControl('Name',"elevon",'EtaInner',0.2,'EtaOuter',0.9,'ChordFraction',0.3);
wing.discretise('nStrips',30);

wing.addMass('AreaDensity',1);
wing.addBoundProps('Radius',4.5*0.0254,'Eta',0.5,'Spin',1);

wingL = wing.mirror();

% vertical tail ==================================================

Vstab = StripModel("fin");
Vstab.setPlanform('Span',0.25,'RootChord',0.15,'TaperRatio',0.6);
Vstab.setLocation('RootQC',[-0.5; 0; 0],'Dihedral',90);   % vertical
Vstab.addControl('Name',"rudder",'EtaInner',0.1,'EtaOuter',0.9,'ChordFraction',0.35);
Vstab.discretise('nStrips',12);

Vstab.addMass('AreaDensity',0.5);


% horizontal tail ================================================

Hstab = StripModel("Hstab");
Hstab.setPlanform('Span',0.25,'RootChord',0.15,'TaperRatio',0.6,'Incidence',-1);
Hstab.setLocation('RootQC',[-0.5; 0; 0],'Dihedral',0);   % horz
Hstab.addControl('Name',"elv",'EtaInner',0.1,'EtaOuter',0.9,'ChordFraction',0.35);
Hstab.discretise('nStrips',12);

Hstab.addMass('AreaDensity',0.5);

HstabL = Hstab.mirror();


% fuselage ============================================

Fusl = StripModel("Hstab");
Fusl.setPlanform('Span',0.10,'RootChord',.6,'TaperRatio',1.2);
Fusl.setLocation('RootQC',[0.1; 0; 0.1],'Dihedral',90);   % horz
Fusl.discretise('nStrips',12);

Fusl.addMass('AreaDensity',0.5);

%% Define the (intial) aerodynamic model 

ModelName = "ViternaCorrigan";  % Model name for the aerodynamic model

%% Propeller:


% ---- 1. set up the prop ----
propparams.D   = 12*0.0254;   % diameter [m]  (12 in)
propparams.P   = 6*0.0254;    % pitch [m]     (6 in)
propparams.B   = 2;           % number of blades
propparams.rho = 1.225;       % air density [kg/m^3]
propparams.V_inf = 0;
propparams.axis = wing.BoundProps.Axis;
propparams.hubloc = wing.BoundProps.Hub;
propparams.spin = 1;


% Remove
% % ---- 2. run BEMT at an operating point (rpm, V_inf) ----
% prop = propellerBEMT(4000, propparams);   % 6000 rpm, hover (V_inf = 0)
% 
% % ---- 3. sample the wake where the wing is ----
% s_wing = 0.25 * propparams.D;      % wing 0.25 diameters behind the disk [m]
% S = prop.develop(s_wing);
% FF  = prop.develop(1e3);
% 
% % radial induced-flow profile for the wing-interaction calc:
% figure
% hold on
% plot(prop.rR, prop.va, prop.rR, prop.vt);
% plot(S.r ./ S.R , S.va, S.r ./ S.R , S.vt);
% plot(prop.rR, FF.va, prop.rR, FF.vt); 
% 
% legend('axial - immediate','swirl','axial - at wing','swirl');
% xlabel('r/R'); ylabel('induced velocity [m/s]');

%disp(prop.T)
%disp(prop.T_axis)
%disp(prop.M_axis)

%% Implementation of Aero



%wing.tiltSurf(70)
%wingL.tiltSurf(70)

%StripPlotter(wing, wingL, Vstab, Hstab, HstabL,Fusl,'Title',"tilt-rotor");




%=====================================================================

% what I want the 'run' script to do:

    % takes in control inuts and u,v,w and p,q,r
    % 1.
    % take in the defined vehicle, which is avaiable after vehicle config
    % runs: includes surfaces, props, extra inertias etc.
    % 2. 
    % tilt the surfaces (if applicable)
    % 3. 
    % evaluate vehicle inertia
    % 4.
    % apply prop thrust
    % 5. 
    % run aerodynamics model, compute and return forcing
    
    

%% Package the vehicle once, after configuration

    plane = struct();
    plane.wing   = wing;
    plane.wingL  = wingL;
    plane.Vstab  = Vstab;
    plane.Hstab  = Hstab;
    plane.HstabL = HstabL;
    plane.Fusl   = Fusl;

    plane.propParams     = propparams;      % D, P, B, rho, spin
    plane.wingWakeDist   = 0.05;          % 0.25*propparams.D, from your config
    plane.rpmFromThrottle = @(t) 50 + 12.6*980*t;   % placeholder - tune this
    plane.addedInertia   = [];                 % or a struct(Mass, MassLoc, Iv)
    plane.ModelName = ModelName;


    % needed for analysis
    plane.refArea = 2*sum(wing.Strips.area);
    plane.refChord = wing.Strips.MAC;



    
    % % test
    % 
    % ctrl      = [0; 0; 1; 0; 0; 0; 0];   % elevonR, elevonL, elv, thr , thrl, rudder, tiltAngle
    % bodyRates = [0; 0; 0];
    % bodyVel   = [15; 0; 0];
    % 
    % 
    % yaw = 0;
    % pitch = 0;
    % roll = 0;
    % orQuat = compact(quaternion([yaw pitch roll],"eulerd","zyx","frame"));
    % 
    % height = -0.1;
    % 
    % [TotalF, TotalM, Mass, CG, Iv] = runVehicle(plane, ctrl, bodyRates, bodyVel, orQuat, height, 1.22, ModelName);
    % 
    % TotalM(1)
    % TotalM(2)
    % TotalF(3)
    