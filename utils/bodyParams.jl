
mutable struct body
    name
    color
    μ # gravitational parameter  [km³/s²]
    m # mass                     [kg]
    a # semimajor axis           [km]
    ω # rotation rate            [rad/s]
    T # period                   [sec]
    r # radius                   [km]
    e # eccentricity             [NON]
end

earth = body("Earth",[0.,0.,1.],398600,5.97237e24,1.495978740473789e8 ,1.160576151685878e-5,3.155814910224000e7, 6371.0084,0)
moon = body("Moon",[.1,.1,.1],4901.801,7.34767309e22,384400,2.661666501955582e-6,2.3606208e6,1737.5,0.0554)
