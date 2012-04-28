% CALC_PQT - Calculates motor mechanical power, reactive power, breakdown
% torque and efficiency from equivalent circuit parameters (used for double
% cage model with core losses)
% 
% Usage: calc_pqt (slip,x)
%
% Where sf is the full load slip (pu)
%       x is a 8 x 1 vector of motor equivalent parameters:
%           x = [Rs Xs Xm Rr1 Xr1 Rr2 Xr2 Rc]
%            x(1) = Rs = stator resistance
%            x(2) = Xs = stator reactance
%            x(3) = Xm = magnetising reactance
%            x(4) = Rr1 = rotor / inner cage resistance
%            x(5) = Xr1 = rotor / inner cage reactance
%            x(6) = Rr2 = outer cage resistance
%            x(7) = Xr2 = outer cage reactance
%            x(8) = Rc = core resistance
%              
% Outputs: y is a vector [Pm Q Tb I_nl]
%

function [y] = calc_pqt(sf,x)

x = abs(x);

[T_fl i_s] = get_torque(sf,x);          % Calculate full-load torque and current
Pm = T_fl * (1 - sf);                   % Calculate mechanical power (at FL)
Sn = complex(1,0)*conj(i_s);
Q_fl = abs(imag(Sn));                   % Calculate reactive power input (at FL)
i_c = 1 / complex (x(8),0);             % Calculate core loss currents (at FL)
i_in = i_s + i_c;                       % Calculate total input current (at FL)
p_in = real(complex(1,0)*conj(i_in));   % Calculate input power (at FL)
eff_fl = Pm / p_in;                     % Calculate efficiency (at FL)

% Calculate breakdown torque with an interval search
T_b = 0;
for i=0.01:0.01:1
    T_i = get_torque(i,x);                                
    if T_i > T_b
        T_b = T_i;                      % Calculated breakdown torque
    end
end

[T_lr i_lr] = get_torque(1,x);
y = [Pm Q_fl T_b T_lr abs(i_lr + i_c) eff_fl];

end