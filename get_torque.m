% GET_TORQUE - Calculate double cage motor torque and stator current
% (without core loss component)
% Author: Julius Susanto
% 
% Usage: get_torque (slip,type,x)
%
% Where slip is the motor slip (pu)
%       x is a vector of motor equivalent parameters:
%           x = [Rs Xs Xm Rr1 Xr1 Rr2 Xr2]
%            Rs = stator resistance
%            Rs = stator reactance
%            Xm = magnetising reactance
%            Rr1 = rotor / inner cage resistance
%            Xr1 = rotor / inner cage reactance
%            Rr2 = outer cage resistance
%            Xr2 = outer cage reactance       
%             
% Outputs: motor torque (pu) as a real number and stator current (as a
% complex number and without core loss component)
%

function [torque ist] = get_torque(slip,x)

% Calculate admittances
Ys = 1/complex(x(1),x(2));
Ym = 1/complex(0,x(3));
Yr1 = 1/complex(x(4)/slip,x(5));
Yr2 = 1/complex(x(6)/slip,x(7));

% Calculate voltage and currents
u1 = Ys / (Ys + Ym + Yr1 + Yr2);
ir1 = abs (u1 * Yr1);
ir2 = abs (u1 * Yr2);

% Calculate torque and stator current
torque = x(4)/slip * ir1^2 + x(6)/slip * ir2^2;
ist = (1 - u1) * Ys;

end