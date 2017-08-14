% pyGAS from MATLAB -- illustrate use of the pyGAS functions from a MATLAB
% environment

clear
clc
close 'all'

% put location of pyGAS.py module on the python search path
% here, I assume that pyGAS.py is in the current directory
if count(py.sys.path,'') == 0  % <-- see if current directory is on path
    insert(py.sys.path,int32(0),''); %<-- if not; add it.
end

workFluid = py.pyGAS.He(); % get a handle to my pyGAS-derived object

r_p = 4;
eta_comp = 0.85;
eta_turb = 0.94;

T = nan(4,1);
T_s = nan(4,1);
P = nan(4,1);
s = nan(4,1);
s_s = nan(4,1);
h = nan(4,1);
h_s = nan(4,1);

% state point 1 - compressor inlet
T(1) = 278; % K
P(1) = 1000; % kPa
h(1) = workFluid.enthalpy(T(1));
s(1) = workFluid.entropy(T(1),P(1));

% compression 
P(2) = r_p*P(1);
T_s(2) = workFluid.isentropicWork(T(1),P(1),P(2));
h_s(2) = workFluid.enthalpy(T_s(2));
h(2) = h(1) - (h(1) - h_s(2))/eta_comp;

w_comp = h(1) - h(2);

% state point 2 - 3: isobaric heat addition to 972 K
T(3) = 972; % K - given
P(3) = P(2); % isentropic
h(3) = workFluid.enthalpy(T(3));
q_s = h(3) - h(2);
s(3) = workFluid.entropy(T(3),P(3));

% state point 3 - 4: expansion in a turbine
P(4) = P(1);
T_s(4) = workFluid.isentropicWork(T(3),P(3),P(4));
h_s(4) = workFluid.enthalpy(T_s(4));
h(4) = h(3) - eta_turb*(h(3) - h_s(4));

w_turbine = h(3) - h(4);

% state point 4 - 1: isobaric heat rejection
q_r = h(1) - h(4);

% 1st law balance
w_net = w_turbine + w_comp;
q_net = q_s + q_r;
eta_th = w_net/q_s;

fprintf('Net work: %g kJ/kg \n',w_net);
fprintf('Net heat: %g kJ/kg \n',q_net);
fprintf('Thermal efficiency: %4.1f percent \n',eta_th*100);



