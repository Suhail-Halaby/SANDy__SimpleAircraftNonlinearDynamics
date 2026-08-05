function [F,M,maxAlpha,minAlpha] = aeromodel(obj,LocCG,Ctrl,BodyVel,BodyRates,rho,ModelName,Va,Vt,rR)
% aeromodel  Outputs forces and moments in the body frame.
    arguments
        obj
        LocCG (3,1) double
        Ctrl (1,1) double
        BodyVel (3,1) double
        BodyRates (3,1) double
        rho (1,1) double
        ModelName string {mustBeMember(ModelName,["simplePlate","ViternaCorrigan"])} = "simplePlate"
        Va (1,:) double = []    
        Vt (1,:) double = []
        rR (1,:) double = [] 
    end


    % 0. Check actuator saturation
    for i = length(obj.Controls)
        if Ctrl < obj.Controls(i).Limits(1) || Ctrl > obj.Controls(i).Limits(2)
            warning(['Control Surface Input Saturated. Affected surface: ', char(obj.Controls.Name)]);
            if  Ctrl < obj.Controls(i).Limits(1)
                Ctrl =  obj.Controls(i).Limits(1);
                warning([char(obj.Controls.Name),'  saturated at Min'])
            else 
                Ctrl =  obj.Controls(i).Limits(2);
                warning([char(obj.Controls.Name),'  saturated at Max'])
            end
        end
    end
    % 1. Velocities from rigid body motion

    % distance from cg:
    r = obj.Strips.Colloc - LocCG;
    % induced velocity (add trhustwash later)
    adjBodyRates = BodyRates.*ones(size(r));
    vLocal = BodyVel + cross(adjBodyRates,r,1);

    % 2. Velocities from propwash
% 2. Velocities from propwash
    if isprop(obj,"BoundProps") && ~isempty(obj.BoundProps) && ~isempty(Va)

        Pc = obj.Strips.Colloc;                  % 3 x N strip collocation pts (body frame)
        N  = size(Pc,2);
        vWash = zeros(3,N);
    
        for k = 1:numel(obj.BoundProps)  % useful if multiple props per elem introduced
            bp   = obj.BoundProps(k);
            hub  = bp.Hub(:);
            axs  = bp.Axis(:)/norm(bp.Axis);     % unit thrust/forward axis
            Rk   = bp.Radius;
            spin = bp.Spin;                      % +1 / -1
            cov  = bp.Cover;
            rRg  = rR(:).';                    % 1 x M  r/R stations for Va,Vt
    
            d    = Pc - hub;                     % 3 x N  hub -> strip
            sAx  = axs.' * d;                    % 1 x N  axial coord along +axis
            dRad = d - axs*sAx;                  % 3 x N  radial component
            rMag = vecnorm(dRad,2,1);            % 1 x N  radial distance
            rRk  = rMag / Rk;                    % 1 x N  strip r/R
    
            % radial & tangential unit vectors (guard the axis singularity)
            safe = rMag > 1e-9;
            rHat = zeros(3,N);
            rHat(:,safe) = dRad(:,safe) ./ rMag(safe);
            tHat = cross(repmat(axs,1,N), rHat, 1);   % tangential = axis x rHat
    
            % interpolate the induced-velocity profiles onto each strip's r/R
            % (0 outside the blade span -> strips beyond the tip get no wash)
            VaS = interp1(rRg, Va, rRk, 'linear', 0);
            VtS = interp1(rRg, Vt, rRk, 'linear', 0);
    
            % immersion: only strips inside the tube get wash
            outside = rRk > 1;
            VaS(outside) = 0;  VtS(outside) = 0;
    


            % strip-relative velocity added by the wake (see sign note below)

            vWash = vWash + cov.*( axs*VaS - spin * (tHat .* VtS) );
        end
    
        vLocal = vLocal + vWash;
    end



    % resolve chordwise and normal
    vChord = dot(vLocal,obj.Strips.chordHat,1);
    vNormal = dot(vLocal,obj.Strips.normalHat,1);

    alpha   = atan2(vNormal, vChord);          % local AoA in the chord-normal plane [rad]
    V2_2D   = vChord.^2 + vNormal.^2;          % exclude spanwise
    qStrip  = 0.5 * rho * V2_2D;               % per-strip dynamic pressure [Pa]

    % model inputs
    delta = -obj.Strips.ctrlIdx * Ctrl * obj.Controls.Direction;
    chord =  obj.Strips.chord;
    surfChord = obj.Strips.ctrlChordFrac.*chord;

    % select lifting model
    if ModelName == "simplePlate"
        [Ls,Ds,Ms] = obj.simplePlate(alpha, delta, chord, surfChord, qStrip);
    elseif ModelName == "ViternaCorrigan"
        [Ls,Ds,Ms] = obj.ViternaCorrigan(alpha, delta, chord, surfChord, qStrip);
    end

    maxAlpha = rad2deg(max(abs(alpha)));
    minAlpha = rad2deg(min(abs(alpha)));

    stripSpan = diff(obj.Strips.edges)*obj.Span;
    Ls = Ls .* stripSpan;
    Ds = Ds .* stripSpan;
    Ms = Ms .* stripSpan;

    localF = -Ds.*obj.Strips.chordHat - Ls.*obj.Strips.normalHat;
    F = sum(localF,2);
    M = sum(cross(r,localF,1),2) - sum(Ms.*obj.Strips.spanHat,2);
end
