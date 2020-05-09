close all
clear all
% setGraphicsDefaults()

%% Numerical Example 3 Buses Without Battery

profile on 

P_inj = sdpvar(3,1,'full','real'); % P injected vector for all buses
Q_inj = sdpvar(3,1,'full','real'); % Q injected vector for all buses
p_load = sdpvar(3,1,'full','real'); % P load vector for all buses
q_load = sdpvar(3,1,'full','real'); % Q load vector for all buses
p_bat = sdpvar(3,1,'full','real'); % P battery vector for all buses
q_bat = sdpvar(3,1,'full','real'); % Q battery vector for all buses
V_real = sdpvar(3,1,'full','real'); % real part of Voltage vector for all buses
DeltaV_real = sdpvar(3,1,'full','real'); % real part of DeltaV vector for all buses
V_imag = sdpvar(3,1,'full','real'); % imaginary part of Voltage vector for all buses
DeltaV_imag = sdpvar(3,1,'full','real'); % imaginary part of DeltaV vector for all buses
% P_inj_constraint = sdpvar(3,1,'full','real');
% Q_inj_constraint = sdpvar(3,1,'full','real');

V_totComplex = V_real+i*V_imag;
DeltaV_totComplex = DeltaV_real+i*DeltaV_imag;

%Define Admittance Matrix 
gkm = [0.1, 0.1, 0.1; 0.1, 0.1, 0.1; 0.1, 0.1, 0.1;];
bkm = 0.2*[1, 1, 1; 1, 1, 1; 1, 1, 1;];
gkm_shunt = 0.03*[1, 1, 1; 1, 1, 1; 1, 1, 1;];
bkm_shunt = 0.04*[1, 1, 1; 1, 1, 1; 1, 1, 1;]; 

Y(1,1) = gkm_shunt(1,2)+i*bkm_shunt(1,2)+gkm(1,2)+i*bkm(1,2);
Y(1,2) = gkm(1,2)+i*bkm(1,2);
Y(1,3) = 0; %zero because bus 1 not connected with bus 3
Y(2,1) = gkm(1,2)+i*bkm(1,2);
Y(2,2) = gkm_shunt(2,1)+i*bkm_shunt(2,1)+gkm(2,1)+i*bkm(2,1)+gkm_shunt(2,3)+i*bkm_shunt(2,3)+gkm(2,3)+i*bkm(2,3);
Y(2,3) = gkm(2,3)+i*bkm(2,3);
Y(3,1) = 0; %zero because bus 1 not connected with bus 3
Y(3,2) = gkm(3,2)+i*bkm(3,2);
Y(3,3) = gkm_shunt(3,2)+i*bkm_shunt(3,2)+gkm(3,2)+i*bkm(3,2);

G = real(Y);
B = imag(Y);

maxCurrent = 100; 
Imax = maxCurrent*ones(3,3); % assuming same Imax for all lines

VminValue = 0;
VmaxValue = 100;
Vmin = VminValue*ones(3,1);
Vmax = VmaxValue*ones(3,1);

%Define Constraints

Gamma_real = diag(G*V_real - B*V_imag);
Gamma_imag = diag(-(G*V_imag+B*V_real));
Ksi_real = diag(V_real)*G+diag(V_imag)*B;
Ksi_imag = diag(V_imag)*G-diag(V_real)*B;
alpha = G*V_real-B*V_imag;
beta = G*V_imag+B*V_real;
Pi_real = -diag(V_real)*alpha - diag(V_imag)*beta;
Pi_imag = -diag(V_imag)*alpha - diag(V_real)*beta;

%Power Constraints (both leq and geq so that only equal counts but maybe better way to define equality???)
P_contstraint1 = P_inj - ((Gamma_real+Ksi_real)*DeltaV_real+(-Gamma_imag+Ksi_imag)*DeltaV_imag-Pi_real+p_load) == 0; 
% P_contstraint2 = P_inj - (Gamma_real+Ksi_real)*DeltaV_real+(-Gamma_imag+Ksi_imag)*DeltaV_imag-Pi_real >= 0;
Q_contstraint1 = Q_inj - ((Gamma_imag+Ksi_imag)*DeltaV_real+(Gamma_real+Ksi_real)*DeltaV_imag-Pi_imag+q_load) == 0;
% Q_contstraint2 = Q_inj - (Gamma_imag+Ksi_imag)*DeltaV_real+(Gamma_real+Ksi_real)*DeltaV_imag-Pi_imag >= 0;

%Current Constraints
Current1_contstraint = abs((V_totComplex(1)-V_totComplex(2))*(gkm(1,2)+i*bkm(1,2))+V_totComplex(1)*(gkm_shunt(1,2)+i*bkm_shunt(1,2)))<=Imax(1,2); %constraint on I12
Current2_contstraint = abs((V_totComplex(2)-V_totComplex(1))*(gkm(1,2)+i*bkm(1,2))+V_totComplex(2)*(gkm_shunt(1,2)+i*bkm_shunt(1,2)))<=Imax(2,1); %constraint on I21
Current3_contstraint = abs((V_totComplex(2)-V_totComplex(3))*(gkm(2,3)+i*bkm(2,3))+V_totComplex(2)*(gkm_shunt(2,3)+i*bkm_shunt(2,3)))<=Imax(2,3); %constraint on I23
Current4_contstraint = abs((V_totComplex(3)-V_totComplex(2))*(gkm(2,3)+i*bkm(2,3))+V_totComplex(3)*(gkm_shunt(2,3)+i*bkm_shunt(2,3)))<=Imax(3,2); %constraint on I32

%Voltage Constraints
Voltage_constraintHigh = V_totComplex<=Vmax;
Voltage_constraintLow = V_totComplex>=Vmin;

%Gathering constraints
Constraints=[P_contstraint1;
              Q_contstraint1;
              Current1_contstraint;
              Current2_contstraint;
              Current3_contstraint;
              Current4_contstraint;
              Voltage_constraintLow];

          
%Objective Function
Objective = abs(P_inj(1))+abs(Q_inj(1))+P_inj(1);

%Define Options for solver
options = sdpsettings('verbose',1,'solver','IPOPT');

%Solving the problem
sol = optimize(Constraints,Objective,options);

try
    PinjectionsValue = value(P_inj);
    QinjectionsValue = value(Q_inj);
    ploadValue = value(Q_inj);
    qloadValue = value(Q_inj);
    Obj_value= value(Objective);
    fprintf('\nPinj: %f\n',PinjectionsValue)
    fprintf('Qinj: %f\n',QinjectionsValue)
    fprintf('Pload: %f\n',ploadValue)
    fprintf('Qload: %f\n',qloadValue)
    fprintf('Obj: %f\n',Obj_value)
catch
    disp('No such value defined to print')
end

profile viewer

%% Numerical Example 3 Buses Without Battery Over Time (24 Timesteps) (Essentially row wise addition of all times instead of 3x1 vectors)

profile on

TimeInterval = 24;

P_inj = sdpvar(3,TimeInterval,'full','real'); % P injected vector for all buses
Q_inj = sdpvar(3,TimeInterval,'full','real'); % Q injected vector for all buses
p_load = sdpvar(3,TimeInterval,'full','real'); % P load vector for all buses
q_load = sdpvar(3,TimeInterval,'full','real'); % Q load vector for all buses
p_bat = sdpvar(3,TimeInterval,'full','real'); % P battery vector for all buses
q_bat = sdpvar(3,TimeInterval,'full','real'); % Q battery vector for all buses
V_real = sdpvar(3,TimeInterval,'full','real'); % real part of Voltage vector for all buses
DeltaV_real = sdpvar(3,TimeInterval,'full','real'); % real part of DeltaV vector for all buses
V_imag = sdpvar(3,TimeInterval,'full','real'); % imaginary part of Voltage vector for all buses
DeltaV_imag = sdpvar(3,TimeInterval,'full','real'); % imaginary part of DeltaV vector for all buses
% P_inj_constraint = sdpvar(3,1,'full','real');
% Q_inj_constraint = sdpvar(3,1,'full','real');
S_dp = sdpvar(1,1,'full','complex');

V_totComplex = V_real+i*V_imag;
DeltaV_totComplex = DeltaV_real+i*DeltaV_imag;

%Define Admittance Matrix 
gkm = [0.1, 0.1, 0.1; 0.1, 0.1, 0.1; 0.1, 0.1, 0.1;];
bkm = 0.2*[1, 1, 1; 1, 1, 1; 1, 1, 1;];
gkm_shunt = 0.03*[1, 1, 1; 1, 1, 1; 1, 1, 1;];
bkm_shunt = 0.04*[1, 1, 1; 1, 1, 1; 1, 1, 1;];

Y(1,1) = gkm_shunt(1,2)+i*bkm_shunt(1,2)+gkm(1,2)+i*bkm(1,2);
Y(1,2) = gkm(1,2)+i*bkm(1,2);
Y(1,3) = 0; %zero because bus 1 not connected with bus 3
Y(2,1) = gkm(1,2)+i*bkm(1,2);
Y(2,2) = gkm_shunt(2,1)+i*bkm_shunt(2,1)+gkm(2,1)+i*bkm(2,1)+gkm_shunt(2,3)+i*bkm_shunt(2,3)+gkm(2,3)+i*bkm(2,3);
Y(2,3) = gkm(2,3)+i*bkm(2,3);
Y(3,1) = 0; %zero because bus 1 not connected with bus 3
Y(3,2) = gkm(3,2)+i*bkm(3,2);
Y(3,3) = gkm_shunt(3,2)+i*bkm_shunt(3,2)+gkm(3,2)+i*bkm(3,2);

G = real(Y);
B = imag(Y);

maxCurrent = 100; 
Imax = maxCurrent*ones(3,3); % assuming same Imax for all lines

VminValue = 0;
VmaxValue = 100;
Vmin = VminValue*ones(3,1);
Vmax = VmaxValue*ones(3,1);

%Define Constraints

Gamma_real = diag(G*sum(V_real,2) - B*sum(V_imag,2));
Gamma_imag = diag(-(G*sum(V_imag,2)+B*sum(V_real,2)));
Ksi_real = diag(sum(V_real,2))*G+diag(sum(V_imag,2))*B;
Ksi_imag = diag(sum(V_imag,2))*G-diag(sum(V_real,2))*B;
alpha = G*sum(V_real,2)-B*sum(V_imag,2);
beta = G*sum(V_imag,2)+B*sum(V_real,2);
Pi_real = -diag(sum(V_real,2))*alpha - diag(sum(V_imag,2))*beta;
Pi_imag = -diag(sum(V_imag,2))*alpha - diag(sum(V_real,2))*beta;

%Power Constraints 
P_contstraint1 = sum(P_inj,2) - ((Gamma_real+Ksi_real)*sum(DeltaV_real,2)+(-Gamma_imag+Ksi_imag)*sum(DeltaV_imag,2)-Pi_real+sum(p_load,2)) == 0; 
% P_contstraint2 = P_inj - (Gamma_real+Ksi_real)*DeltaV_real+(-Gamma_imag+Ksi_imag)*DeltaV_imag-Pi_real >= 0;
Q_contstraint1 = sum(Q_inj,2) - ((Gamma_imag+Ksi_imag)*sum(DeltaV_real,2)+(Gamma_real+Ksi_real)*sum(DeltaV_imag,2)-Pi_imag+sum(q_load,2)) == 0;
% Q_contstraint2 = Q_inj - (Gamma_imag+Ksi_imag)*DeltaV_real+(Gamma_real+Ksi_real)*DeltaV_imag-Pi_imag >= 0;

temp1 = sum(V_totComplex,2);

%Current Constraints
Current1_contstraint = abs((temp1(1)-temp1(2))*(gkm(1,2)+i*bkm(1,2))+temp1(1)*(gkm_shunt(1,2)+i*bkm_shunt(1,2)))<=Imax(1,2); %constraint on I12
Current2_contstraint = abs((temp1(2)-temp1(1))*(gkm(1,2)+i*bkm(1,2))+temp1(2)*(gkm_shunt(1,2)+i*bkm_shunt(1,2)))<=Imax(2,1); %constraint on I21
Current3_contstraint = abs((temp1(2)-temp1(3))*(gkm(2,3)+i*bkm(2,3))+temp1(2)*(gkm_shunt(2,3)+i*bkm_shunt(2,3)))<=Imax(2,3); %constraint on I23
Current4_contstraint = abs((temp1(3)-temp1(2))*(gkm(2,3)+i*bkm(2,3))+temp1(3)*(gkm_shunt(2,3)+i*bkm_shunt(2,3)))<=Imax(3,2); %constraint on I32

%Voltage Constraints
Voltage_constraintHigh = sum(V_totComplex,2)<=Vmax;
Voltage_constraintLow = sum(V_totComplex,2)>=Vmin;

%Gathering constraints
Constraints=[P_contstraint1;
              Q_contstraint1;
              Current1_contstraint;
              Current2_contstraint;
              Current3_contstraint;
              Current4_contstraint;
              Voltage_constraintLow];

          
Sinj = P_inj+i*Q_inj;     
          
%Objective Function
Objective = abs(sum(P_inj(1,:)))+abs(sum(Q_inj(1,:)))+sum(P_inj(1,:));

%Define Options for solver
options = sdpsettings('verbose',1,'solver','IPOPT');

%Solving the problem
sol = optimize(Constraints,Objective,options);

try
    Pinjections = value(P_inj);
    Qinjections = value(Q_inj);
    Obj_value= value(Objective);
    fprintf('\nPinj: %f\n',Pinjections)
    fprintf('Qinj: %f\n',Qinjections)
    fprintf('Obj: %f\n',Obj_value)
catch
    disp('No such value defined to print')
end

profile viewer

%% Numerical Example 3 Buses WITH Battery Over Time (12 Timesteps)

profile on

TimeInterval = 12;
DeltaT = 1;
eta = 0.2;
a_b = 0.5;

P_inj = sdpvar(3,TimeInterval,'full','real'); % P injected vector for all buses
Q_inj = sdpvar(3,TimeInterval,'full','real'); % Q injected vector for all buses
p_load = sdpvar(3,TimeInterval,'full','real'); % P load vector for all buses
q_load = sdpvar(3,TimeInterval,'full','real'); % Q load vector for all buses
p_bat = sdpvar(3,TimeInterval,'full','real'); % P battery vector for all buses
q_bat = sdpvar(3,TimeInterval,'full','real'); % Q battery vector for all buses
V_real = sdpvar(3,TimeInterval,'full','real'); % real part of Voltage vector for all buses
DeltaV_real = sdpvar(3,TimeInterval,'full','real'); % real part of DeltaV vector for all buses
V_imag = sdpvar(3,TimeInterval,'full','real'); % imaginary part of Voltage vector for all buses
DeltaV_imag = sdpvar(3,TimeInterval,'full','real'); % imaginary part of DeltaV vector for all buses
% P_inj_constraint = sdpvar(3,1,'full','real');
% Q_inj_constraint = sdpvar(3,1,'full','real');
S_dp = sdpvar(1,1,'full','complex');
p_bat = sdpvar(3,TimeInterval,'full','real'); % P battery vector for all buses
q_bat = sdpvar(3,TimeInterval,'full','real'); % Q battery vector for all buses
p_bat_plus = sdpvar(3,TimeInterval,'full','real'); % P battery charging vector for all buses
p_bat_minus = sdpvar(3,TimeInterval,'full','real'); % P battery discharging vector for all buses
SoE = sdpvar(3,TimeInterval,'full','real');
E_b = sdpvar(3,TimeInterval,'full','real');

V_totComplex = V_real+i*V_imag;
DeltaV_totComplex = DeltaV_real+i*DeltaV_imag;

%Define Admittance Matrix 
gkm = [0.1, 0.1, 0.1; 0.1, 0.1, 0.1; 0.1, 0.1, 0.1;];
bkm = 0.2*[1, 1, 1; 1, 1, 1; 1, 1, 1;];
gkm_shunt = 0.03*[1, 1, 1; 1, 1, 1; 1, 1, 1;];
bkm_shunt = 0.04*[1, 1, 1; 1, 1, 1; 1, 1, 1;];

Y(1,1) = gkm_shunt(1,2)+i*bkm_shunt(1,2)+gkm(1,2)+i*bkm(1,2);
Y(1,2) = gkm(1,2)+i*bkm(1,2);
Y(1,3) = 0; %zero because bus 1 not connected with bus 3
Y(2,1) = gkm(1,2)+i*bkm(1,2);
Y(2,2) = gkm_shunt(2,1)+i*bkm_shunt(2,1)+gkm(2,1)+i*bkm(2,1)+gkm_shunt(2,3)+i*bkm_shunt(2,3)+gkm(2,3)+i*bkm(2,3);
Y(2,3) = gkm(2,3)+i*bkm(2,3);
Y(3,1) = 0; %zero because bus 1 not connected with bus 3
Y(3,2) = gkm(3,2)+i*bkm(3,2);
Y(3,3) = gkm_shunt(3,2)+i*bkm_shunt(3,2)+gkm(3,2)+i*bkm(3,2);

G = real(Y);
B = imag(Y);

maxCurrent = 100; 
Imax = maxCurrent*ones(3,3); % assuming same Imax for all lines

VminValue = 0;
VmaxValue = 100;
Vmin = VminValue*ones(3,1);
Vmax = VmaxValue*ones(3,1);

SoEmaxValue = 1000;
SoEminValue = 0;
SoEmax = SoEmaxValue*ones(3,TimeInterval); 
SoEmin = SoEminValue*ones(3,TimeInterval);

SbValue = 1000;
Sb_rated = SbValue*ones(3,TimeInterval);

E_bminValue = 0;
E_bmaxValue = 1000;
E_bmin = E_bminValue*ones(3,TimeInterval); 
E_bmax = E_bmaxValue*ones(3,TimeInterval); 


%Define Constraints

Gamma_real = diag(G*sum(V_real,2) - B*sum(V_imag,2));
Gamma_imag = diag(-(G*sum(V_imag,2)+B*sum(V_real,2)));
Ksi_real = diag(sum(V_real,2))*G+diag(sum(V_imag,2))*B;
Ksi_imag = diag(sum(V_imag,2))*G-diag(sum(V_real,2))*B;
alpha = G*sum(V_real,2)-B*sum(V_imag,2);
beta = G*sum(V_imag,2)+B*sum(V_real,2);
Pi_real = -diag(sum(V_real,2))*alpha - diag(sum(V_imag,2))*beta;
Pi_imag = -diag(sum(V_imag,2))*alpha - diag(sum(V_real,2))*beta;

%Power Constraints 
P_contstraint1 = sum(P_inj,2) - ((Gamma_real+Ksi_real)*sum(DeltaV_real,2)+(-Gamma_imag+Ksi_imag)*sum(DeltaV_imag,2)-Pi_real+sum(p_load,2)) == 0; 
% P_contstraint2 = P_inj - (Gamma_real+Ksi_real)*DeltaV_real+(-Gamma_imag+Ksi_imag)*DeltaV_imag-Pi_real >= 0;
Q_contstraint1 = sum(Q_inj,2) - ((Gamma_imag+Ksi_imag)*sum(DeltaV_real,2)+(Gamma_real+Ksi_real)*sum(DeltaV_imag,2)-Pi_imag+sum(q_load,2)) == 0;
% Q_contstraint2 = Q_inj - (Gamma_imag+Ksi_imag)*DeltaV_real+(Gamma_real+Ksi_real)*DeltaV_imag-Pi_imag >= 0;

temp1 = sum(V_totComplex,2);

%Current Constraints
Current1_contstraint = abs((temp1(1)-temp1(2))*(gkm(1,2)+i*bkm(1,2))+temp1(1)*(gkm_shunt(1,2)+i*bkm_shunt(1,2)))<=Imax(1,2); %constraint on I12
Current2_contstraint = abs((temp1(2)-temp1(1))*(gkm(1,2)+i*bkm(1,2))+temp1(2)*(gkm_shunt(1,2)+i*bkm_shunt(1,2)))<=Imax(2,1); %constraint on I21
Current3_contstraint = abs((temp1(2)-temp1(3))*(gkm(2,3)+i*bkm(2,3))+temp1(2)*(gkm_shunt(2,3)+i*bkm_shunt(2,3)))<=Imax(2,3); %constraint on I23
Current4_contstraint = abs((temp1(3)-temp1(2))*(gkm(2,3)+i*bkm(2,3))+temp1(3)*(gkm_shunt(2,3)+i*bkm_shunt(2,3)))<=Imax(3,2); %constraint on I32

%Voltage Constraints
Voltage_constraintHigh = sum(V_totComplex,2)<=Vmax;
Voltage_constraintLow = sum(V_totComplex,2)>=Vmin;

%Battery (Written Very Inefficiently TODO: IMPROVE!!!)
SoE_constraint1 = SoE(:,2) == SoE(:,1) + (eta*p_bat_plus(:,1) - (1/eta)*p_bat_minus(:,1))*DeltaT;
SoE_constraint2 =SoE(:,3) == SoE(:,2) + (eta*p_bat_plus(:,2) - (1/eta)*p_bat_minus(:,2))*DeltaT;
SoE_constraint3 =SoE(:,4) == SoE(:,3) + (eta*p_bat_plus(:,3) - (1/eta)*p_bat_minus(:,3))*DeltaT;
SoE_constraint4 =SoE(:,5) == SoE(:,4) + (eta*p_bat_plus(:,4) - (1/eta)*p_bat_minus(:,4))*DeltaT;
SoE_constraint5 =SoE(:,6) == SoE(:,5) + (eta*p_bat_plus(:,5) - (1/eta)*p_bat_minus(:,5))*DeltaT;
SoE_constraint6 =SoE(:,7) == SoE(:,6) + (eta*p_bat_plus(:,6) - (1/eta)*p_bat_minus(:,6))*DeltaT;
SoE_constraint7 =SoE(:,8) == SoE(:,7) + (eta*p_bat_plus(:,7) - (1/eta)*p_bat_minus(:,7))*DeltaT;
SoE_constraint8 =SoE(:,9) == SoE(:,8) + (eta*p_bat_plus(:,8) - (1/eta)*p_bat_minus(:,8))*DeltaT;
SoE_constraint9 =SoE(:,10) == SoE(:,9) + (eta*p_bat_plus(:,9) - (1/eta)*p_bat_minus(:,9))*DeltaT;
SoE_constraint10 =SoE(:,11) == SoE(:,10) + (eta*p_bat_plus(:,10) - (1/eta)*p_bat_minus(:,10))*DeltaT;
SoE_constraint11 =SoE(:,12) == SoE(:,11) + (eta*p_bat_plus(:,11) - (1/eta)*p_bat_minus(:,11))*DeltaT;

SoEbound1max = SoE<= (1-a_b)*SoEmax; %collective notation (can only be done when all values same (?))
SoEbound1min = SoE>= a_b*SoEmin;
Eq13_constraint = p_bat_plus.^2 + p_bat_minus.^2 + q_bat.^2 <= Sb_rated.^2; %collective notation
Eq14_constraint = E_bmin - SoE <= E_b;
Eq15_constraint = 0 <= E_b;
Eq16_constraint = - E_bmax - SoE <= E_b;


%Gathering constraints
Constraints=[P_contstraint1;
              Q_contstraint1;
              Current1_contstraint;
              Current2_contstraint;
              Current3_contstraint;
              Current4_contstraint;
              Voltage_constraintLow;
              SoE_constraint1;
              SoE_constraint2;
              SoE_constraint3;
              SoE_constraint4;
              SoE_constraint5;
              SoE_constraint6;
              SoE_constraint7;
              SoE_constraint8;
              SoE_constraint9;
              SoE_constraint10;
              SoE_constraint11;
              SoEbound1max;
              SoEbound1min;
              Eq13_constraint;
              Eq14_constraint;
              Eq15_constraint;
              Eq16_constraint];
          
Sinj = P_inj+i*Q_inj;     
          
%Objective Function
Objective = abs(sum(P_inj(1,:)))+abs(sum(Q_inj(1,:)))+sum(P_inj(1,:)) + sum(sum(E_b));

%Define Options for solver
options = sdpsettings('verbose',1,'solver','IPOPT');

%Solving the problem
sol = optimize(Constraints,Objective,options);

try
    Pinjections = value(P_inj);
    Qinjections = value(Q_inj);
    Obj_value= value(Objective);
    fprintf('\nPinj: %f\n',Pinjections)
    fprintf('Qinj: %f\n',Qinjections)
    fprintf('Obj: %f\n',Obj_value)
catch
    disp('No such value defined to print')
end

profile viewer
