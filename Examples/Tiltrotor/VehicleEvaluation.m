clear 
clc
close all

%% Constructor


addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..'));  % reach the root
setupPaths();

analysis = RunAnalysis("VehicleConfig", "runVehicle"); 

%% Coefficient Sweep


U_inf = 20;
Angle = 'alpha';
Range = -5:25;

opts.Plot = true;
opts.ctrl = zeros(7,1);

optsCell = namedargs2cell(opts);

[F, M, Cf, Cm] = analysis.CoefficientSweep(U_inf, Angle, Range, optsCell{:});
%   or equivalently: CoefficientSweep(analysis, U_inf, Angle, Range, opts)

%% Linearization

U_e = 0;
alpha_e  = deg2rad(2);
theta_e = 0;
ctrl_e = [0;0;0;0.5;0.5;0;90];
opts2 = struct;


% 
[A,x_e,f_e,As,Aa] = analysis.SystemJacobian(U_e,alpha_e, theta_e, ctrl_e);
[B,bx_e,bf_e,Bs,Ba] = analysis.CtrlJacobian(U_e,alpha_e, theta_e, ctrl_e);


AllocSurf = [1:2,4:6];

% trim: 

%[ctrlTrim, A, B, resNorm, f_e, As, Bs, Aa, Ba] = analysis.TrimSolve(AllocSurf,U_e,alpha_e,theta_e,ctrl_e);
%(resNorm)

disp(f_e)

eig(As)
eig(Aa)

%% Stab Check

ActiveSurf = [1:7];
% Decoupled subsystems, only through the surfaces you actually allocate
% (map full-B allocSurf columns onto the partitioned Bs/Ba as needed)
analysis.Controllability(As, Bs(:,ActiveSurf), Label="Longitudinal");
analysis.Controllability(Aa, Ba(:,ActiveSurf), Label="Lateral");

% Full 9-state, all controls
analysis.Controllability(A, B(:,AllocSurf), Label="Full 9-DOF");

%% Eigensolve


%% Eigensweep

% specify variable
