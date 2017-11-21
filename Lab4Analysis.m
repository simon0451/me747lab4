clear all; close all;

%% Part 1 a)
bucket = 15.9;
weights = [0,20,40,60,80,100,120,140,160,180,200]+bucket; % grams
weightsOzf = weights/1000*9.81*0.224809*16; % convert to ozf

incV = [475,501,527,554,582,607,633,658,684,712,737]/1000; % [V]
decV = [473,503,528,552,580,607,633,659,686,711,738]/1000; % [V]

% plot eo vs weight
figure(1);
hold on;
plot(weightsOzf,incV,'k--');
plot(weightsOzf,decV,'b');
xlabel('Weights (Oz_f)','FontSize',12);
ylabel('Voltage Output e_o (V)','FontSize',12);
legend('Location','best','Voltage from increasing deflection',...
    'Voltage from decreasing deflection');

% The hysterisis is low so the average is used for the sensitivity
avgV = (incV + decV)/2;
weightVFit = polyfit(weights/1000*9.81,avgV,1);
sensWeight = weightVFit(1); % 0.0372 V/N
sensWeightEng = sensWeight/3.5969431019354;

%% Part 1 b)
eo_vs_time = importdata('data 1_4.lvm','\t',28);
eo_time = eo_vs_time.data(:,1);
eo_volt = eo_vs_time.data(:,2);

for i = 1:length(eo_volt)
    if (eo_volt(i) < 0.7645)
        startIndex = i-3;
        break;
    end
end
eo_time = eo_time(startIndex-100:3500) - eo_time(startIndex);

eo_volt = eo_volt(startIndex-100:3500);
endMean = mean(eo_volt(1300:end));
% eo_volt = eo_volt - endMean;

figure(2);
plot(eo_time,eo_volt);
xlabel('Time (s)','FontSize',12);
ylabel('Voltage (V)','FontSize',12);

%% Part 1 c) - e)
[~,minV] = peakdet(eo_volt,0.05);
osPeak = minV(2); % peak voltage
tp = eo_time(minV(1)); % peak time
ospct = (endMean - osPeak)/endMean; % overshoot percentage

% damping ratio from the graph given overshoot percentage
zeta = 0.55;

% From peak time to omegaN
omegaN = pi/(tp*sqrt(1 - zeta^2)); % [rad/s]

% change in voltage after moving to the side
delta_eo = 60/1000; % [V]
% equivalent weight from sensitivity
delta_weight = delta_eo/sensWeight; % N
% from weight to mass: N -> kg
delta_mass_eff = delta_weight/9.81;
delta_mass_eff_eng = delta_mass_eff*35.27392; % [kg] -> [ozm]

% from weight_eff and omegaN to spring constant
K = omegaN^2*delta_mass_eff; % [kg/s^2] = [N/m]
KEng = K*3.5969431019354/39.3701; % [N/m] -> [Ozf/in]

% From spring constant to B using omegaN and zeta:
B = 2*zeta/omegaN*K;
BEng = 2*zeta/omegaN*KEng;

% sensitivity: V/N -> V/m/s^2 = V/(kg-m/s^2) * kg
sensAcc = sensWeight*delta_mass_eff;
sensAccEng = sensAcc/39.3701; % V/m/s^2 -> V/ in/s^2

% maximum acceleration
maxAcc = 778/1000/sensAcc; % [m/s^2] ~ 12g
maxAccEng = maxAcc*39.3701; % [m/s^2] -> [in/s^2]

% TF: MREs/2RoK / (M/K s^2 + B/K s + 1): MREs/2RoK -> steady-state gain
% DC potentiometer accelerometer, bandwidth 100 Hz, sensitivity = gain
sys = tf(sensAccEng,[delta_mass_eff_eng/KEng BEng/KEng 1]);
figure(3);
bode(sys);

%% Part 2.1 a)
mass = [0.2,0.4,0.6,0.8,1,1.2,1.4,1.6]; % [kg]
eo = [21.1,42.8,63.4,83.6,104.7,125.68,148.56,173.35]; % [mV]

weight = mass*9.81; % [N]
weightEng = weight*0.224809; % [N] -> [lbf]

calFit = polyfit(weight,eo,1);
calSens = calFit(1); % mV/N
calSensEng = calSens/0.224809; % [mV/N] -> [mV/lbf]
calV = calSens*weight + calFit(2);

maxdiffPct = max(eo - calV)/calV(end)*100;

figure(4);
plot(weightEng,eo);
xlabel('Mass (kg)','FontSize',12);
ylabel('Voltage (mV)','FontSize',12);

%% Part 2.1 b)
decayData = importdata('data 2_1_3.lvm','\t',29);
decayTime = decayData.data(:,1);
decayVolt = decayData.data(:,2);

decaySmooth = smooth(smooth(decayVolt));
decaySmooth = decaySmooth - decaySmooth(end);

startVoltage = 0.1222; % pure speculation
startTime = 0.06422;

tauVoltage = (decayVolt(end)-startVoltage)*0.632 + startVoltage;

for i = 1:length(decayTime)
    if (decaySmooth(i) < tauVoltage && decayTime(i) > startTime)
        tauPoint = i;
        break;
    end
end

tauTime = decayTime(tauPoint);

figure(5);
hold on;
plot(decayTime,decaySmooth);
plot([0 tauTime],[tauVoltage tauVoltage],'r--');
plot([tauTime tauTime],[-0.02 tauVoltage],'r--');
xlim([0 35]);
ylim([0 0.14]);
xlabel('Time (s)','FontSize',12);
ylabel('Voltage (V)','FontSize',12);

%% Part 2.1 c)
PCBfrequencyData = importdata('data 2_1_4.lvm','\t',29);
wnTime = PCBfrequencyData.data(:,1);
wnVolt = PCBfrequencyData.data(:,2);

[maxV,~] = peakdet(wnVolt(1000:3500),0.03);

timeDiffPCB = diff(wnTime(maxV(:,1)+999));
timeDiffPCBAvg = mean(timeDiffPCB);
PCBWn = 1/timeDiffPCBAvg*2*pi; % [rad/s]

figure(6);
hold on;
plot(wnTime,wnVolt);
plot(wnTime(maxV(:,1)+999),wnVolt(maxV(:,1)+999),'O');
xlabel('Time (s)','FontSize',12);
ylabel('Voltage (V)','FontSize',12);

%% Part 2.2 a) and b)
foamMassData = importdata('data 2_2_1.lvm','\t',29);
foamTime = foamMassData.data(3000:end,1);
foamVolt = foamMassData.data(3000:end,2);

foamTime = foamTime - 0.613;

foamVolt = foamVolt - mean(foamVolt(1:100));

[maxV,~] = peakdet(foamVolt,0.002);

TdFoam = diff(foamTime(maxV(:,1)));
TdFoam = TdFoam(1);

n = 3;
sigmaFoam = (1/(n-1))*log(maxV(2,2)/maxV(4,2));
zetaFoam = sigmaFoam/sqrt(4*pi^2 + sigmaFoam^2);

wdFoam = 2*pi/TdFoam;
wnFoam = wdFoam/sqrt(1 - zetaFoam^2);

% 2 lbf plate = 8.89644 N = 8.89644/9.81 kg
mFoam = 8.89644/9.81; % [kg]
mFoamEng = 2/32.2; % [lbm]

KFoam = wnFoam^2*mFoam; % [kg/s^2] = N/m
KFoamEng = KFoam*0.2248089431/39.3701; % [lbf/in]

BFoam = 2*zetaFoam/wnFoam*KFoam;
BFoamEng = 2*zetaFoam/wnFoam*KFoamEng; % [lbf/ in/s]

foamSys = tf(0.007, [1/(wnFoam^2) 2*zetaFoam/wnFoam 1]);
figure(7);
impulse(foamSys);
hold on;
plot(foamTime,foamVolt*1000/calSensEng,'r');

%% Part 3 a)
lvtAcceleration = importdata('data 2_3_4 - accel.lvm','\t',32);
%lvtForce = importdata('data 2_3_3 - force.lvm','\t',32);
lvtA = lvtAcceleration.data;
%lvtF = lvtForce.data;

lvtA = lvtA(1:3000,:);

% reset the start time
for i = 1:length(lvtA(:,1))
    if (lvtA(i,1) > -0.05029) % number pulled from the plot
        lvtStart = i;
        break;
    end
end

lvtA(:,[1 3]) = lvtA(:,[1 3]) - lvtA(lvtStart,1);

% Find the time when the body experience full acceleration
accStart = 0;
accEnd = 0;
for i = lvtStart:length(lvtA(:,1))
    if (lvtA(i,1) > 0.029 && accEnd == 0)
        accEnd = i;
        break;
    elseif (lvtA(i,1) > 0.0062 && accStart == 0) % number pulled from the plot
        accStart = i;
    end
end

avgAccV = mean(lvtA(accStart:accEnd,2)); % [V] This number should be gravity
accSens = avgAccV/-9.81; % V/ m/s^2
accSensEng = accSens/39.3701; % [V/ m/s^2] -> [V/ in/s^2]

% total change in velocity should equal to the average acceleration
% (gravity) times the time span, use the period that I am confident that
% the data is good
dTime = lvtA(accEnd,1) - lvtA(accStart,1);
dVel = -9.81*dTime; % [m/s]
dVelEng = dVel*39.3701; % [in/s]
dVelV = lvtA(accEnd,4) - lvtA(accStart,4); % channel 1 is LVT
lvtSens = dVelV/dVel; % [V/ m/s]
lvtSensEng = dVelV/dVelEng; % [V/ in/s]

figure(8);
hold on;
plot(lvtA(:,1),lvtA(:,2),'b--');
plot(lvtA(:,3),lvtA(:,4),'r');
xlabel('Time (s)','FontSize',12);
ylabel('LVT and Accelerometer Output (V)','FontSize',12);

%% Part 3 b)
intAcc = cumtrapz(lvtA(:,1),lvtA(:,2)/accSensEng);

figure(9);
hold on;
plot(lvtA(:,1),intAcc,'--');
plot(lvtA(:,3),lvtA(:,4)/lvtSensEng);
xlabel('Time (s)','FontSize',12);
ylabel('Velocity (in/s)','FontSize',12);

%% Part 3 c)
figure(10);
hold on;
plot(lvtA(:,1),lvtA(:,2),'b');
plot(lvtA(:,3),lvtA(:,4),'r');
plot([0 0],[-3 2],'k--'); % beginning of drop
plot([0.07 0.07],[-3 2],'r--'); % first contact with foam
plot([0.0909 0.0909],[-3 2],'b--'); % maximum displacement in foam
plot([0.11 0.11],[-3 2],'g--'); % leaving the foam
plot([0.2494 0.2494],[-3 2],'--','Color',[0.5 0.5 0.5]); % Permanent on foam
xlabel('Time (s)','FontSize',12);
ylabel('LVT and Accelerometer Output (V)','FontSize',12);

%% Part 3 d)
maxV = max(abs(lvtA(:,4)/lvtSensEng)); % [in/s] going downwards

%% Part 3 e)
% lvtF = lvtF(1:3000,:);
fSensEng = 0.491; % [V/lbf]
fSens = fSensEng/4.44822; % [V/N]

[~,ind] = max(abs(lvtF(:,4)));
maxF = abs(lvtF(ind,2)-lvtF(1,2))/fSens; % [N] maximum force at maximum velocity
ssF = abs(lvtF(end,2)-lvtF(1,2))/fSens; % [N] force at steady state
totalM = ssF/9.81; % [kg]
totalMEng = totalM*2.20462; % [lbm]

figure(11);
hold on;
plot(lvtF(:,1),lvtF(:,2),'b--');
plot(lvtF(:,3),lvtF(:,4),'r');
xlabel('Time (s)','FontSize',12);
ylabel('LVT and Force Sensor Output (V)','FontSize',12);