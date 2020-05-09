close all
clear all
% setGraphicsDefaults()

%% Numerical Example 3 Buses Without Battery (3x1 Vectors) (Not Realistic - Easier to test)

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
P_dp = sdpvar(1,1,'full','real'); %Real part of S^dp which is the last part (with weight w5) in eq 16
Q_dp = sdpvar(1,1,'full','real'); %Imag part of S^dp which is the last part (with weight w5) in eq 16

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
% Objective = abs(P_inj(1))+abs(Q_inj(1))+P_inj(1); %not including last part of obj (weight w5) 
Objective = abs(P_inj(1))+abs(Q_inj(1))+P_inj(1)+(P_inj(1)-P_dp)^2+(Q_inj(1)-Q_dp)^2; %including last part of obj (weight w5) 

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


%% Numerical Example 3 Buses WITH Battery Over Time (12 Timesteps)

profile on

TimeInterval = 12;
DeltaT = 1;
eta = 0.2;
a_b = 0.5;

p_load = [0 0 0 0 0 0 0 0 0 0 0 0; 1 5 3 14 18 26 35 24 22 25 32 12; 10 12 17 24 25 30 32 42 45 38 35 32;];
q_load = 0.1*[0 0 0 0 0 0 0 0 0 0 0 0; 1 5 3 14 18 26 35 24 22 25 32 12; 10 12 17 24 25 30 32 42 45 38 35 32;];

P_inj = sdpvar(3,TimeInterval,'full','real'); % P injected vector for all buses
Q_inj = sdpvar(3,TimeInterval,'full','real'); % Q injected vector for all buses
% p_load = sdpvar(3,TimeInterval,'full','real'); % P load vector for all buses
% q_load = sdpvar(3,TimeInterval,'full','real'); % Q load vector for all buses
% p_bat = sdpvar(3,TimeInterval,'full','real'); % P battery vector for all buses
% q_bat = sdpvar(3,TimeInterval,'full','real'); % Q battery vector for all buses
V_real = sdpvar(3,TimeInterval,'full','real'); % real part of Voltage vector for all buses
DeltaV_real = sdpvar(3,TimeInterval,'full','real'); % real part of DeltaV vector for all buses
V_imag = sdpvar(3,TimeInterval,'full','real'); % imaginary part of Voltage vector for all buses
DeltaV_imag = sdpvar(3,TimeInterval,'full','real'); % imaginary part of DeltaV vector for all buses
% P_inj_constraint = sdpvar(3,1,'full','real');
% Q_inj_constraint = sdpvar(3,1,'full','real');
% S_dp = sdpvar(1,1,'full','complex');
p_bat = sdpvar(3,TimeInterval,'full','real'); % P battery vector for all buses
q_bat = sdpvar(3,TimeInterval,'full','real'); % Q battery vector for all buses
p_bat_plus = sdpvar(3,TimeInterval,'full','real'); % P battery charging vector for all buses
p_bat_minus = sdpvar(3,TimeInterval,'full','real'); % P battery discharging vector for all buses
SoE = sdpvar(3,TimeInterval,'full','real');
E_b = sdpvar(3,TimeInterval,'full','real');
P_dp = sdpvar(1,TimeInterval,'full','real'); %Real part of S^dp which is the last part (with weight w5) in eq 16
Q_dp = sdpvar(1,TimeInterval,'full','real'); %Imag part of S^dp which is the last part (with weight w5) in eq 16

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

maxCurrent = 80; 
Imax = maxCurrent*ones(3,3); % assuming same Imax for all lines

VminValue = 0.9;
VmaxValue = 1.1;
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
ConstraintsNew = [];
%loop to create constraints for each time ("for every bus" doesnt need a loop since implemented as rows of V matrix)

for it = 1: TimeInterval
    Gamma_real = diag(G*V_real(:,it) - B*V_imag(:,it));
    Gamma_imag = diag(-(G*V_imag(:,it)+B*V_real(:,it)));
    Ksi_real = diag(V_real(:,it))*G+diag(V_imag(:,it))*B;
    Ksi_imag = diag(V_imag(:,it))*G-diag(V_real(:,it))*B;
    alpha = G*(V_real(:,it))-B*(V_imag(:,it));
    beta = G*(V_imag(:,it))+B*(V_real(:,it));
    Pi_real = -diag((V_real(:,it)))*alpha - diag((V_imag(:,it)))*beta;
    Pi_imag = -diag((V_imag(:,it)))*alpha - diag((V_real(:,it)))*beta;
    
    %Power Constraints
    P_contstraint1 = (P_inj(:,it)) - ((Gamma_real+Ksi_real)*(1-V_real(:,it))+(-Gamma_imag+Ksi_imag)*(1-V_imag(:,it))-Pi_real+(p_load(:,it))) == 0;
    % P_contstraint2 = P_inj - (Gamma_real+Ksi_real)*DeltaV_real+(-Gamma_imag+Ksi_imag)*DeltaV_imag-Pi_real >= 0;
    Q_contstraint1 = (Q_inj(:,it)) - ((Gamma_imag+Ksi_imag)*(1-V_real(:,it))+(Gamma_real+Ksi_real)*(1-V_imag(:,it))-Pi_imag+(q_load(:,it))) == 0;
    % Q_contstraint2 = Q_inj - (Gamma_imag+Ksi_imag)*DeltaV_real+(Gamma_real+Ksi_real)*DeltaV_imag-Pi_imag >= 0;
    
    VComplex = V_totComplex - V_real - i*V_imag == 0; %Assuming No DV in this example
    VatPCC = V_totComplex(1,:) - ones(1,12) - i*0 ==0; %Voltage at slack bus considered 1
    
    temp1 = V_totComplex(:,it);
    
    %Current Constraints
    Current1_contstraint = abs((temp1(1)-temp1(2))*(gkm(1,2)+i*bkm(1,2))+temp1(1)*(gkm_shunt(1,2)+i*bkm_shunt(1,2)))<=Imax(1,2); %constraint on I12
    Current2_contstraint = abs((temp1(2)-temp1(1))*(gkm(1,2)+i*bkm(1,2))+temp1(2)*(gkm_shunt(1,2)+i*bkm_shunt(1,2)))<=Imax(2,1); %constraint on I21
    Current3_contstraint = abs((temp1(2)-temp1(3))*(gkm(2,3)+i*bkm(2,3))+temp1(2)*(gkm_shunt(2,3)+i*bkm_shunt(2,3)))<=Imax(2,3); %constraint on I23
    Current4_contstraint = abs((temp1(3)-temp1(2))*(gkm(2,3)+i*bkm(2,3))+temp1(3)*(gkm_shunt(2,3)+i*bkm_shunt(2,3)))<=Imax(3,2); %constraint on I32
    
    %Voltage Constraints
    Voltage_constraintHigh = (temp1)<=Vmax;
    Voltage_constraintLow = (temp1)>=Vmin;
    
    ConstraintsNew = [ConstraintsNew; P_contstraint1; Q_contstraint1; Current1_contstraint; Current2_contstraint; 
        Current3_contstraint; Current4_contstraint; Voltage_constraintHigh; Voltage_constraintLow; VComplex; VatPCC]; 

end

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

NoPBatInjPCC = p_bat(3,:) == 0; %no P battery injection at first bus (PCC)
NoQBatInjPCC = q_bat(3,:) == 0; %no Q battery injection at first bus (PCC)
NoPBatInjPCCplus = p_bat_plus(3,:) == 0; 
NoPBatInjPCCminus = p_bat_minus(3,:) == 0; 
% NoPBatInjPCC3 = p_bat(3,:) == 0; %no P battery injection at bus 3
% NoQBatInjPCC3 = q_bat(3,:) == 0; %no Q battery injection at bus 3

%loop to create constraints for each time ("for every bus" doesnt need a loop since implemented as rows of V matrix)

for it = 1: TimeInterval
    
    SoEbound1max = SoE(:,it)<= (1-a_b)*SoEmax(:,it); 
    SoEbound1min = SoE(:,it)>= a_b*SoEmin(:,it);
    Eq13_constraint = p_bat_plus(:,it).^2 + p_bat_minus(:,it).^2 + q_bat(:,it).^2 <= Sb_rated(:,it).^2; 
    Eq14_constraint = E_bmin(:,it) - SoE(:,it) <= E_b(:,it);
    Eq15_constraint = 0 <= E_b(:,it);
    Eq16_constraint = - E_bmax(:,it) - SoE(:,it) <= E_b(:,it);
    
    ConstraintsNew = [ConstraintsNew; SoEbound1max; SoEbound1min; Eq13_constraint; Eq14_constraint; Eq15_constraint; Eq16_constraint; NoPBatInjPCC; NoQBatInjPCC;
        NoPBatInjPCCplus; NoPBatInjPCCminus];
    
end


%Gathering constraints
% Constraints=[P_contstraint1;
%               Q_contstraint1;
%               Current1_contstraint;
%               Current2_contstraint;
%               Current3_contstraint;
%               Current4_contstraint;
%               Voltage_constraintLow;
%               SoE_constraint1;
%               SoE_constraint2;
%               SoE_constraint3;
%               SoE_constraint4;
%               SoE_constraint5;
%               SoE_constraint6;
%               SoE_constraint7;
%               SoE_constraint8;
%               SoE_constraint9;
%               SoE_constraint10;
%               SoE_constraint11;
%               SoEbound1max;
%               SoEbound1min;
%               Eq13_constraint;
%               Eq14_constraint;
%               Eq15_constraint;
%               Eq16_constraint];

ConstraintsNew = [ConstraintsNew;
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
              SoE_constraint11];
          
% Sinj = P_inj+i*Q_inj;     
          
%Define objective function weights (same indexes refer to same parts as in CoDistFlow paper)
w1 = 1; 
w2 = 1;
w3 = 1;
w4 = 1;
w5 = 1;

%Objective Function
% Objective = abs(sum(P_inj(1,:)))+abs(sum(Q_inj(1,:)))+sum(P_inj(1,:)) + sum(sum(E_b));
Objective = w3*abs(sum(P_inj(1,:))) + w2*abs(sum(Q_inj(1,:))) + w4*sum(P_inj(1,:)) + w1*sum(sum(E_b)) + w5*(sum(P_inj(1,:))-sum(P_dp(1,:)))^2 + w5*(sum(Q_inj(1,:))-sum(Q_dp(1,:)))^2;

%Define Options for solver
options = sdpsettings('verbose',1,'solver','IPOPT');

%Solving the problem
sol = optimize(ConstraintsNew,Objective,options);

try
%     Pinjections = value(P_inj);
%     Qinjections = value(Q_inj);
%     Obj_value= value(Objective);
%     fprintf('\nPinj: %f\n',Pinjections)
%     fprintf('Qinj: %f\n',Qinjections)
%     fprintf('Obj: %f\n',Obj_value)
catch
    disp('No such value defined to print')
end

profile viewer

%% Plotting 

timeTotal = linspace(1,12,12);

%Plot p_load at each bus
figure 
for i = 1:size(p_load,1)
    plot(timeTotal, p_load(i,:))
    hold on
    grid on
end
% title('')
legend('PCC','Bus 2','Bus 3')
xlabel('Time [h]')
ylabel('P_{load} [W]')

figure 
for i = 1:size(p_load,1)
    plot(timeTotal, q_load(i,:))
    hold on
    grid on
end
% title('')
legend('PCC','Bus 2','Bus 3')
xlabel('Time [h]')
ylabel('Q_{load} [W]')

%Plot V complex at each bus
figure 
for i = 1:size(p_load,1)
    plot(timeTotal, value(V_totComplex(i,:)))
    hold on
    grid on
end
% title('')
legend('Bus 1','Bus 2','Bus 3')
xlabel('Time [h]')
ylabel('Voltage [V]')

%Plot P_inj at each bus
figure 
for i = 1:size(p_load,1)
    plot(timeTotal, value(P_inj(i,:)))
    hold on
    grid on
end
% title('')
legend('PCC','Bus 2','Bus 3')
xlabel('Time [h]')
ylabel('P [V]')

%Plot P_bat_plus at each bus
figure 
for i = 1:size(p_load,1)
    plot(timeTotal, value(p_bat_plus(i,:)))
    hold on
    grid on
end
% title('')
legend('PCC','Bus 2','Bus 3')
xlabel('Time [h]')
ylabel('P_{bat plus} [V]')

%Plot P_bat_minus at each bus
figure 
for i = 1:size(p_load,1)
    plot(timeTotal, value(p_bat_minus(i,:)))
    hold on
    grid on
end
% title('')
legend('PCC','Bus 2','Bus 3')
xlabel('Time [h]')
ylabel('P_{bat minus} [V]')
