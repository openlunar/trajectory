function setEarthMoon(ln)
% setEarthMoon(ln)
% 
% model for Earth-Moon CR3BP system
% 
% ln defines the lagrange point number as origin, ln=0 sets origin as CR3BP
% origin (system barycenter). ln=0 is the default if none is specified.

global mu BODY EARTH MOON SUN LN RL PRIM SEC AU RUNIT TUNIT VUNIT AUNIT

disp('Set EARTH-MOON System');
SUN.name        = 'SUN';
SUN.gm          = 132712440017.986999511718750; % [km^3/s^2]
SUN.radius      = 6.955e5;                      % [km] (approximate value)

G = 6.67408e-20;                                % [km^3/(kg*s^2)]

SEC.name        = 'MOON';
SEC.gm          = 4902.801;                     % [km^3/s^2]  (+-0.001)   http://ssd.jpl.nasa.gov/?sat_phys_par
SEC.radius      = 1737.5;                       % [km]        (+-0.1  )   http://ssd.jpl.nasa.gov/?sat_phys_par
SEC.sm          = 384400;                       % [km]        (+-?    )   http://ssd.jpl.nasa.gov/?sat_elem
SEC.ecc         = 0.0554;                       % [unitless]  (+-?    )   http://ssd.jpl.nasa.gov/?sat_elem
SEC.period      = 27.322 * 86400;               % [sec]       (+-?    )   http://ssd.jpl.nasa.gov/?sat_elem
SEC.rot         = 2 * pi / SEC.period;          % [rad/s] synchronous
SEC.color       = [0.6,0.6,0.6];                % Grey

PRIM.name       = 'EARTH';
PRIM.m          = 5.97237e24;                   % [kg]        (+-0.00028)
PRIM.gm         = G*PRIM.m;                     % [km^3/s^2]  (+-2.7  )   https://ssd.jpl.nasa.gov/?planet_phys_par
PRIM.gmbary     = PRIM.gm + SEC.gm;             % [km^3/s^2]  (+-2.7  )   https://ssd.jpl.nasa.gov/?planet_phys_par
PRIM.radius     = 6371.0084;                    % [km]        (+-4    )   http://ssd.jpl.nasa.gov/?planet_phys_par
PRIM.period     = 1.0000174 * 365.25 * 86400;   % [sec]       (+-?    )   http://ssd.jpl.nasa.gov/?planet_phys_par
PRIM.sm         = ((SUN.gm + PRIM.gmbary) * (PRIM.period / (2 * pi)) ^ 2) ^ (1 / 3); % km
PRIM.rot        = 1 / 0.99726968;               % [Rev/Day]   (+-?    )   http://ssd.jpl.nasa.gov/?planet_phys_par
PRIM.rotD       = PRIM.rot * 360 / 86400;       % [Deg/Sec]   (+-?    )
PRIM.color      = [0,0,1];                      % Blue

AU              = 149597927.000;                % km

mu              = SEC.gm / (SEC.gm + PRIM.gm);

BODY            = SEC;
MOON            = SEC;
EARTH           = PRIM;

RUNIT           = SEC.sm;
TUNIT           = SEC.period / (2 * pi);
VUNIT           = RUNIT / TUNIT;
AUNIT           = RUNIT / TUNIT ^ 2;

if nargin == 0
    ln          = 0;
    rL          = [0;0;0];
else
    if ~(ln == 0||ln == 1||ln == 2||ln == 3||ln == 4||ln == 5)
        error('lagrange point muste be integer 0<=ln<=5')
    elseif ln == 0
        rL      = zeros(3,1);
    else
        [Xeq, ~]= CR3BPEqPts(mu);
        rL      = Xeq(1:3, ln);
    end
end

% assign global variables
LN              = ln;
RL              = rL;

% Record of Revision
% Date          Programmer          Description of Changes
% 01/23/2020    Jared T. Blanchard  Code modified from Brian Anderson