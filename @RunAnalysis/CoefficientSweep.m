function [F, M, Cf, Cm, Cw] = CoefficientSweep(obj, U_inf, Angle, Range, opts)
    arguments
        obj
        U_inf (1,1) double {mustBePositive}
        Angle (1,1) string {mustBeMember(Angle, ["alpha","beta"])}
        Range (1,:) double
        opts.ctrl      (:,1) double  = zeros(7,1)   % [elevonR elevonL elv thr thrl rudder tilt]
        opts.bodyRates (3,1) double  = zeros(3,1)   % [p q r]
        opts.height    (1,1) double  = 100
        opts.rho       (1,1) double  = 1.225
        opts.subtractGravity (1,1) logical = true
        opts.Plot      (1,1) logical = false
        opts.WindPlot  (1,1) logical = true
        opts.DegreeInput (1,1) logical = true
    end

    % --- 1. Build the vehicle once (isolated workspace) ---
    plane = obj.loadVehicle();
    assert(isfield(plane,'ModelName'), "plane.ModelName missing after buildVehicle");

    % --- 2. Reference quantities ---
    S = plane.refArea;
    c = plane.refChord;
    if isfield(plane, 'refSpan')
        b = plane.refSpan;
    else
        b = c;   % fallback
        warning("RunAnalysis:NoRefSpan", ...
            "plane.refSpan not defined; using refChord for roll/yaw coefficients.");
    end

    % --- Normalise the sweep range ONCE ---
    % Range_rad drives the trig; Range_deg is what the plots display.
    if opts.DegreeInput
        Range_deg = Range;
        Range_rad = deg2rad(Range);
    else
        Range_rad = Range;
        Range_deg = rad2deg(Range);
    end

    g    = 9.81;
    qInf = 0.5 * opts.rho * U_inf^2;

    % Sweep is imposed through the velocity vector; vehicle stays level.
    orQuat = [1 0 0 0];   % identity: level, wings-level  [w x y z]

    % --- 3. Sweep ---
    N = numel(Range_rad);
    F = zeros(3, N);
    M = zeros(3, N);
    alphas = zeros(1, N);
    betas  = zeros(1, N);

    for k = 1:N
        ang = Range_rad(k);
        if Angle == "alpha"
            alpha = ang;   beta = 0;
        else
            alpha = 0;     beta = ang;
        end

        % Freestream -> body velocity [u v w] at fixed speed U_inf
        u = U_inf * cos(alpha) * cos(beta);
        v = U_inf * sin(beta);
        w = U_inf * sin(alpha) * cos(beta);
        bodyVel = [u; v; w];

        [Ftot, Mtot, Mass, ~, ~] = obj.RunFcn(plane, opts.ctrl, opts.bodyRates, bodyVel, ...
             orQuat, opts.height, opts.rho, plane.ModelName);


        % ---------- TEMP DIAGNOSTIC ----------
        % fprintf('RunScriptName = "%s"  (class: %s)\n', obj.RunScriptName, class(obj.RunScriptName));
        % fprintf('ctrl      %s %s\n', class(opts.ctrl),       mat2str(size(opts.ctrl)));
        % fprintf('bodyRates %s %s\n', class(opts.bodyRates),  mat2str(size(opts.bodyRates)));
        % fprintf('bodyVel   %s %s\n', class(bodyVel),         mat2str(size(bodyVel)));
        % fprintf('orQuat    %s %s\n', class(orQuat),          mat2str(size(orQuat)));
        % fprintf('height    %s %s\n', class(opts.height),     mat2str(size(opts.height)));
        % fprintf('rho       %s %s\n', class(opts.rho),        mat2str(size(opts.rho)));
        % fprintf('ModelName %s %s\n', class(plane.ModelName), mat2str(size(plane.ModelName)));
        % which('runVehicle', '-all')
        % 
        % % Bypass feval — call runVehicle DIRECTLY with the same variables:
        % [Ftot, Mtot, Mass, ~, ~] = runVehicle(plane, opts.ctrl, opts.bodyRates, bodyVel, ...
        %     orQuat, opts.height, opts.rho, plane.ModelName);
        % ---------- END DIAGNOSTIC ----------

        Ftot = Ftot(:);
        Mtot = Mtot(:);

        % Remove weight from the forces (level -> gravity is +z body)
        if opts.subtractGravity
            Ftot = Ftot - [0; 0; Mass * g];
        end

        F(:,k)    = Ftot;
        M(:,k)    = Mtot;
        alphas(k) = alpha;
        betas(k)  = beta;
    end

    % --- 4. Non-dimensionalise (body axes) ---
    Cf = F ./ (qInf * S);
    Cm = [ M(1,:) ./ (qInf * S * b);    % C_l  (roll,  span)
           M(2,:) ./ (qInf * S * c);    % C_m  (pitch, chord)
           M(3,:) ./ (qInf * S * b) ];  % C_n  (yaw,   span)

    % --- 5. Body-axis plot ---
    if opts.Plot
        obj.plotSweep(Range_deg, Angle, Cf, Cm);
    end

    % --- Wind-axis coefficients: [CD; CY; CL] ---
    Cw = zeros(3, N);
    for k = 1:N
        Fw      = obj.bodyToWind(F(:,k), alphas(k), betas(k));  % [D; Yw; L]
        Cw(:,k) = Fw ./ (qInf * S);
    end

    if opts.WindPlot
        obj.plotWind(Range_deg, Angle, Cw, Cm);
    end
end