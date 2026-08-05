function [M, cgV, Iv] = assembleInertia(varargin)
% assembleInertia  Combine StripModel parts into vehicle mass, CG, and inertia
% tensor (about the vehicle CG, body axes).
%
%   [M,cgV,Iv] = assembleInertia(wingR, wingL, fin)
%
% Each part must have had addMass() run first, populating .Inertias with:
%   .mass (1x1)   .cg (3x1)   .I (3x3, about own CG, body axes)

    N     = numel(varargin);
    m     = zeros(1,N);
    cg    = zeros(3,N);
    Ipart = cell(1,N);

    for i = 1:N
        S = varargin{i};
        m(i)     = S.PartMass;
        cg(:,i)  = S.PartCG;
        Ipart{i} = S.Inertias;
    end

    % ---- vehicle mass and CG ----
    M   = sum(m);
    cgV = (cg * m.') / M;                      % 3 x 1

    % ---- parallel-axis each part into the vehicle frame, then sum ----
    Iv = zeros(3);
    for i = 1:N
        d  = cg(:,i) - cgV;
        Iv = Iv + Ipart{i} + m(i) * ((d.'*d)*eye(3) - d*d.');
    end
end