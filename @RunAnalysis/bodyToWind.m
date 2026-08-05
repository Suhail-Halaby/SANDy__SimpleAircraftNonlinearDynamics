function Fw = bodyToWind(~, Fb, alpha, beta)
% Body [X;Y;Z] (X fwd, Y stbd, Z down) -> wind [D; Yw; L].
ca = cos(alpha); sa = sin(alpha);
cb = cos(beta);  sb = sin(beta);

% Body-to-wind DCM
T = [ ca*cb,   sb,   sa*cb ;
    -ca*sb,   cb,  -sa*sb ;
    -sa,      0,    ca    ];

Fwind = T * Fb;          % [along-wind ; side ; wind-down]
D  = -Fwind(1);          % drag opposes velocity
Yw =  Fwind(2);
L  = -Fwind(3);          % lift is up (-wind z)
Fw = [D; Yw; L];
end