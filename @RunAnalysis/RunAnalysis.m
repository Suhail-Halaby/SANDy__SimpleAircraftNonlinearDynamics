classdef RunAnalysis < handle
    properties
        ConfigScriptName (1,1) string
        RunScriptName    (1,1) string
        RunFcn (1,1) function_handle = @() []
    end

    %% ----------------------------- Class Constructor --------------------
    methods (Access = public)

        function obj = RunAnalysis(ConfigScript, RunScript)
            arguments
                ConfigScript (1,1) string
                RunScript    (1,1) string
            end

            % --- Both scripts must exist on the MATLAB path ---
            if exist(ConfigScript, 'file') ~= 2
                error("RunAnalysis:ConfigNotFound", ...
                      "Config script '%s' not found on the path.", ConfigScript);
            end
            if exist(RunScript, 'file') ~= 2
                error("RunAnalysis:RunNotFound", ...
                      "Run script '%s' not found on the path.", RunScript);
            end

            % --- Config script must actually execute without error ---
            try
                RunAnalysis.runSandboxed(ConfigScript);
            catch ME
                error("RunAnalysis:ConfigInvalid", ...
                      "Config script '%s' failed to run: %s", ...
                      ConfigScript, ME.message);
            end

            % --- Store validated names ---
            obj.ConfigScriptName = ConfigScript;
            obj.RunScriptName    = RunScript;
            obj.RunFcn = str2func(RunScript);   % add a RunFcn property
        end

    end

    %% Analysis functionality (defined in separate files) ----------------
    % the running script takes the form of:
    % runVehicle(plane, ctrl, bodyRates, bodyVel, orQuat, height, rho, ModelName)
    methods (Access = public)
        [F, M, Cf, Cm, Cw]              = CoefficientSweep(obj, U_inf, Angle, Range, opts)
        [A, x_e, f_e, As, Aa]           = SystemJacobian(obj, u_e, alpha_e, theta_e, ctrl_e, opts)
        [B, x_e, f_e, Bs, Ba]           = CtrlJacobian(obj, u_e, alpha_e, theta_e, ctrl_e, opts)
        [ctrlTrim, A, B, resNorm, f_e, As, Bs, Aa, Ba] = ...
            TrimSolve(obj, allocSurf, u_e, alpha_e, theta_e, ctrl_e, solveOpts, linOpts)
        C                               = Controllability(obj, A, B, opts)
    end

    %% Private helpers (defined in separate files) -----------------------
    methods (Access = private)
        plane = loadVehicle(obj)
        plotSweep(obj, Range, Angle, Cf, Cm)
        Fw    = bodyToWind(obj, Fb, alpha, beta)
        plotWind(obj, Range, Angle, Cw, Cm)
        xdot  = stateDeriv(obj, plane, x, ctrl, cfg)
        q     = euler2quat(obj, phi, theta, psi)
    end

    %% Private Static Methods (defined in separate files) ----------------
    methods (Access = private, Static)
        runSandboxed(scriptName)
    end

    % functionality
    % 1. create running file
    %    takes in control inputs + tilts (vehicle must be defined with vehicle config)
    %    spits out forces and moments
    % 2. coefficients (Cl,Cd,Cm) vs. (alpha,beta sweep)
    % 3. linearization (populate stability derivatives, 9-dof linear matrices)
end