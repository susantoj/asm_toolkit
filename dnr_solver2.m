% DNR_SOLVER2 - Damped Newton-Rhapson solver for double cage model with core losses
%              Solves for 6 circuit parameters [Xs Xm Rr1 Xr1 Rr2 Rc]
%              Rs and Xr2 are known
%              Includes change of variables
%              Includes adaptive step size (as per Pedra 2008)
%              Includes determinant check of jacobian matrix
%
% Usage: nr_solver2 (p, Rs, Xr2, lambda, max_iter)
%
% Where   p is a vector of motor performance parameters:
%         p = [sf eff pf Tb Tlr Ilr]
%           sf = full-load slip
%           eff = full-load efficiency
%           pf = full-load power factor
%           T_b = breakdown torque (as % of FL torque)
%           T_lr = locked rotor torque (as % of FL torque)
%           I_lr = locked rotor current
%		  Rs and Xr2 are known parameters
%         max_iter is the maximum number of iterations  
%
% OUTPUT:   x is a vector of motor equivalent parameters:
%           x = [Rs Xs Xm Rr1 Xr1 Rr2 Xr2 Rc]
%            x(1) = Rs = stator resistance
%            x(2) = Xs = stator reactance
%            x(3) = Xm = magnetising reactance
%            x(4) = Rr1 = rotor / inner cage resistance
%            x(5) = Xr1 = rotor / inner cage reactance
%            x(6) = Rr2 = outer cage resistance
%            x(7) = Xr2 = outer cage reactance
%            x(8) = Rc = core resistance
%           iter is the number of iterations
%           err is the squared error of the objective function
%           conv is a true/false flag indicating convergence

function [z iter err conv] = dnr_solver2(p, Rs, Xr2, lambda, max_iter)

% For TESTING ONLY
%p = [0.01 0.947 0.86 2.7 2.5 6.7];
%p = [0.0266 0.87 0.88 3.1 2.5 5.7];
%p = [0.0266 0.88 0.83 3.6 2.95 6];
%p = [0.015 0.96 0.86 2.1 1.7 6];


% Problem sets
%p = [0.027 0.915 0.81 3.300925926 2.199074074 7];
%p = [0.018 0.936 0.86 2.400149589 2.600037397 6.199494949];
%p = [0.04 0.84 0.8 3 2.698630137 3.043478261];

% DIgSI problem sets
% p = [0.0144444 0.935 0.79 2 1.2 5.3];

% Human-readable motor performance parameters
% And base value initialisation
sf = p(1);                          % Full-load slip (pu)
eff = p(2);                         % Full-load efficiency (pu)
pf = p(3);                          % Full-load power factor (pu)
T_fl = pf * eff / (1 - sf);         % Full-load torque (pu)
T_b = p(4) * T_fl;                  % Breakdown torque (pu)
T_lr = p(5) * T_fl;                 % Locked rotor torque (pu)
i_lr = p(6);                        % Locked rotor current (pu)
Pm_fl = pf * eff;                   % Mechanical power (at FL)
Q_fl = sin(acos(pf));               % Full-load reactive power (pu)

% Set initial conditions
z(3) = 1 / Q_fl;            %Xm
z(2) = 0.05 * z(3);         %Xs
z(4) = 1 / Pm_fl * sf;      %Rr1
z(1) = Rs;                  %Rs
z(5) = 1.2 * z(2);          %Xr1
z(6) = 5 * z(4);            %Rr2
z(7) = Xr2;                 %Xr2
z(8) = 12;

% Change of variables to constrained parameters (with initial values)
x(1) = z(4);
x(2) = z(6) - z(4);
x(3) = z(3);
x(4) = z(2);
x(5) = z(5) - z(7);
x(6) = z(8);

% Formulate solution
pqt = [Pm_fl Q_fl T_b T_lr i_lr eff];

% Set up NR algorithm parameters
h = 0.00001;
n = 0;
hn = 1;
hn_min = 0.0000001;
err_tol = 0.00001;
err = 1;
iter = 0;
conv = false;
gamma = 3;
beta = 3;

% Run NR algorithm
while (err > err_tol) && (iter < max_iter)
    
    % Evaluate objective function for current iteration
    y = (pqt - calc_pqt4(sf,z))./pqt;
    err0 = y*y';
    
    % Construct Jacobian matrix
    for i=1:6
        x(i) = x(i) + h;
        
        % Change of variables back to equivalent circuit parameters
        z(2) = x(4);
        z(3) = x(3);
        z(4) = x(1);
        z(5) = z(7) + x(5);
        z(6) = x(1) + x(2);
        z(8) = x(6);
        
        j(:,i) = ((pqt - calc_pqt4(sf,z))./pqt - y) / h;
        x(i) = x(i) - h;
    end
    
    % Check if jacobian matrix is singular and exit function if so
    if (det(j) == 0)
        return;
    end
    
    x_reset = x;
    y_reset = y;
    iter0 = iter;
    
    % Inner loop (descent direction check and step size adjustment)
    while (iter == iter0)
        % Calculate next iteration and update x
        delta_x = (j + lambda.*eye(6))\y';
        x = abs(x - hn * delta_x');

        % Change of variables back to equivalent circuit parameters
        z(2) = x(4);
        z(3) = x(3);
        z(4) = x(1);
        z(5) = z(7) + x(5);
        z(6) = x(1) + x(2);
        z(8) = x(6);
        
        % Calculate squared error terms
        y = (pqt - calc_pqt4(sf,z))./pqt;
        err = y*y';
        
        % Descent direction check and step size adjustment
        if (abs(err) >= abs(err0))
            n = n + 1;
            hn = 2^(-n);
            lambda = lambda * beta;
            x = x_reset;
            y = y_reset;
        else
            n = 0;
            lambda = lambda / gamma;
            iter = iter + 1;
        end
        
        % If descent direction isn't minimising, then there is no convergence
        if (hn < hn_min)
            return; 
        end
   end
end

if err < err_tol
    conv = true;
end

end
