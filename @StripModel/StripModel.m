classdef StripModel < handle
% STRIPMODEL  Geometric definition of a single lifting surface (panel) for
% strip-theory aerodynamic modelling of an aircraft/tailsitter.
%
% ---- Step 1 scope (this file): GEOMETRY ONLY ----
%   - define a lifting surface (one root->tip panel)
%   - planform spec : chord, taper, twist, wing setting angle (+ sweep)
%   - control surfaces : occupied span (eta) + chord fraction
%   - location/orientation : root quarter-chord position rel. vehicle datum
%   - discretise() into spanwise strips (pure geometry)
%
% ---- Conventions ----
%   Body axes : X forward, Y starboard (right), Z down   (standard aerospace)
%   Lengths   : metres
%   Angles    : DEGREES in all properties/inputs (converted internally)
%   Anchor    : root quarter-chord (RootQC), body axes, relative to
%   aircraft datum
%   A panel runs root (eta=0) -> tip (eta=1) along the span axis, which is
%   body +Y rotated by dihedral about body X (tip-up corresponds to -Z).
%   Dihedral = 90 deg  ->  vertical fin.   Use mirror() for symmetric halves.
%
% NOT in this step: aero/force models, control deflection effects, plotting.
%
% ---- File layout (class folder @StripModel) ----
%   This class is split across a class folder. This file holds the
%   properties and constructor; each method body lives in its own file:
%     setPlanform.m  setLocation.m  addControl.m  discretise.m
%     mirror.m  tiltSurf.m  aeromodel.m  summary.m  simplePlate.m (private)
%   MATLAB treats the whole @StripModel folder as one class, so usage is
%   unchanged (e.g. wing = StripModel("main"); wing.discretise()).

    properties
        Name (1,1) string = "surface"

        % ---- Planform (in the surface reference plane) ----
        Span        (1,1) double {mustBePositive} = 1.0   % root->tip length [m]
        RootChord   (1,1) double {mustBePositive} = 0.2   % chord at root [m]
        TaperRatio  (1,1) double {mustBePositive} = 1.0   % c_tip / c_root [-]
        TwistTip    (1,1) double = 0.0                    % geometric twist at tip, washout<0 [deg]
        Incidence   (1,1) double = 0.0                    % wing setting angle at root [deg]
        SweepQC     (1,1) double = 0.0                    % quarter-chord sweep, aft +ve [deg]

        % ---- Location / orientation (relative to body CG) ----
        RootQC      (3,1) double = [0;0;0]                % root quarter-chord pos, body axes [m]
        Dihedral    (1,1) double = 0.0                    % dihedral [deg]; 90 => vertical fin
        SpanSide    (1,1) double {mustBeMember(SpanSide,[-1 1])} = 1  % +1 right, -1 left/mirror

        % ---- Control surfaces (struct array) ----
        Controls (1,:) struct = struct('Name',{},'EtaInner',{},'EtaOuter',{}, ...
                                        'ChordFraction',{},'Deflection',{},'Limits',{},'Direction',{})

        % ---- Discretisation result (filled by discretise) ----
        Strips (1,1) struct = struct()

        % ---- Inertial Data (post discretization) ----
        Inertias (3,3) double = zeros(3)
        PartCG (3,1) double = 0
        PartMass (1,1) double = 0 
        MassSpec  cell   = {}          % replay arguments

        % ---- Propulsive Data (post discretization)  ----
        Prop (1,1) struct = struct()
        BoundProps struct = struct([])   % 1 x nProps prop definitions
        Propwash   double = []           % 1 x n, 1 = strip in propwash
        PropCover  double = []           % 1 x n, fractional coverage


        % ---- Tilt (tiltrotor / all-moving surface) ----
        % Rotation about the body Y-axis by Tilt [deg], applied at discretise().
        % The axis is parallel to body Y and passes through RootQC + TiltOffset,
        % so the panel swings about a pivot offset from its reference point.
        % (Exact when the span axis is parallel to Y, i.e. dihedral ~ 0.)
        Tilt       (1,1) double = 0            % wing setting-angle offset [deg]
        TiltOffset (3,1) double = [0;0;0]      % tilt-axis offset from RootQC [m] (X,Z used)
    end

    % ================= constructor =================
    methods (Access = public)
        function obj = StripModel(name)
            arguments
                name (1,1) string = "surface"
            end
            obj.Name = name;
        end
    end

    % ========= public methods defined in separate files =========
    methods (Access = public)
        
        % initialization
        obj = setPlanform(obj, opts)
        obj = setLocation(obj, opts)
        obj = addControl(obj, opts)

        % discretization yields strip objects
        S   = discretise(obj, opts)
        % specified masses
        addMass(obj,opts)
        refreshMass(obj)
        % add bound propellers
        addBoundProps(obj,opts)
        refreshBoundProps(obj,opts)

        % mirror yields y-mirrored lifting surfaces
        m   = mirror(obj,opts)
        % tilt-rotor command
        obj = tiltSurf(obj, tilt)
       

        % aerodynamics model (here if needed)
        [F,M,maxAlpha,minAlpha] = aeromodel(obj,LocCG,Ctrl,BodyVel,BodyRates,rho,ModelName)
              summary(obj)
    end

    % ========= private methods defined in separate files =========
    methods (Access = private)
        [L,D,M,Cl,Cd,Cm] = simplePlate(obj, alpha, delta, chord, surfChord, qStrip)
        [L,D,M,Cl,Cd,Cm] = ViternaCorrigan(obj, alpha, delta, chord, surfChord, qStrip)
    end
end
