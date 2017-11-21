%Simon Popecki
%ME 747
%13 Novermber 2017
%Lab 4
clear all;
close all;


%% Part 1 Potentiometer Accelerometer
load lab4data.mat
mass_nobucket = [20 40 60 80 100 120 140 160 180 200]; %grams
bucketWeight = 15.9; %grams
testmass = mass_nobucket+bucketWeight; %testmass (g) is the total mass on the sensor
testmassoz = testmass*0.035274; %converting from grams to ounces because Prof. Fussell likes imperial
eoinc = [501 527 554 582 607 633 658 684 712 737]; %mV, output voltage with weight increasing
eodec = [503 528 552 580 607 633  659 686 711 738]; %mV, output voltage with weight decreasing
eoavg = (eoinc+eodec)/2; %averageing the output voltages of increasing and decreasing values
bestfit = polyfit(testmassoz,eoavg,1); %making a best fit line to find sensitivity

sensitivity = bestfit(1)/1000; %the sensitivity of the scale in volts per ounce

% Plotting eo vs. weight
figure(1)
plot(testmassoz,eoinc,testmassoz,eodec)
title('Potentiometer Accelerometer Proof Weight vs. Output Voltage')
xlabel('Proof Weight (Oz)')
ylabel('Output Voltage (mV)')
legend('Increasing Weight','Decreasing Weight','location','northwest')
grid on

% Plotting eo vs. time
figure(2)
plot(time,eout)
title('Accelerometer Voltage Response')
xlabel('Time (s)')
ylabel('Output Voltage (V)')
grid on

bottompeak = .39; %V, from manual inspection of figure 2 - the lowest peak
finalvalue = .4528; %V, from manual inspection of figure 2 - the ending value
percentovershoot = finalvalue/bottompeak; %finding the percent overshoot for use with the table on canvas
zeta = .52; %from the table, estimated value
peaktime = .029; %s, the time to the peak - obtained visually from plot
omegan = pi/(peaktime*(sqrt(1-(zeta^2)))); %rad/s, natural frequency
w = ((.454-.392)/sensitivity); %Oz, equivalent weight
m = w*0.00194256; %slugs, effective mass
k = (16/12)*(omegan^2)*m; %oz/in, the spring constant - !highly dependant on peak time (16/12 for converison from slugs to Oz/in)
b = 2*zeta/omegan*k; %damping ratio, unknown imperial units with ounces and inches

accelsensitivity = ((m*16/12)*sensitivity); %((Oz*s*s)/ft)*(V/Oz) -> V/(in/(s*s)) divided by 12 to convert from feet to inches, 16 to convert to Oz
maxaccel = (.7781/accelsensitivity); %in/(s*s), maximum sensor acceleration (is about 12 and a half Gs) (.7781 is resting voltage)

system = tf(accelsensitivity,[(m*16/12)/k b/k 1]);
figure(3);
bode(system);
title('Output Voltage vs. Acceleration')
grid on

%% Part 2.1 Piezoelectric Force Sensor

%calibration data
%everything is called -two now because it's part two
masstwo = [.2 .4 .6 .8 1 1.2 1.4 1.6]; %kg
weighttwo = (masstwo*2.20462)*16; %converting from kilograms to Earth Oz
eouttwo = [21.1 42.8 63.4 83.6 104.7 125.68 148.56 173.35]; %mV

figure(4)
plot(weighttwo,eouttwo,'-o')
title('Voltage Response of the Piezoelectric Force Sensor')
xlabel('Proof Weight (Oz)')
ylabel('Output Voltage (mV)')
grid on

bestfittwo = polyfit(weighttwo,eouttwo,1); %creating a linear best fit line to determine sensitivity
sensitivitytwo = bestfittwo(1); %mV/Oz, the slope is the sensitivity
sensitivitytwo = sensitivitytwo/1000; %V/Oz, converting from millivolts to volts per ounce

calibrationline = sensitivitytwo*weighttwo+(bestfittwo(2)/1000);
maxdifference = max((eouttwo/1000)-calibrationline)/(calibrationline(end)); %maximum difference determines the error
fserror = maxdifference*100; %full scale error percent

load lab4data.mat

smdvolt = smooth(decayvoltage);
smdvolt = smooth(smdvolt); % !warning! more advanced smoothing techniques require too much computational power

vstart = .11; %from visual inspection of figure 5
tstart = .06291; %s, start time corresponding to vstart

vtau = .632*(smdvolt(end)-vstart)+vstart; %finding the time constant with the 63.2% method

for i = length(decaytime):-1:1
    if smdvolt(i) > vtau %finding the point where voltage exceeds vtau, working in backwards direction
        tauindex = i;
        break;
    end
end

ttau = decaytime(tauindex); %s, finding the point in time of tau

figure(5)
plot(decaytime,smdvolt,ttau,vtau,'o')
title('Piezo Decay Test Plot')
xlabel('Time (s)')
ylabel('Output Voltage (V)')
legend('Voltage Response','First Time Constant')
grid on

figure(6)
plot(omegat,omegav)


deltat = .001277-.0006732; %time difference between points, from inspection of figure 6 - second peaks
dampednaturalfrequency = 1/deltat*2*pi; %rad/s

%% Part 2.2 Piezoelectric Force Sensor

load lab4data.mat
avgv = mean(dvoltage(1:200)); %finding an average voltage to adjust for shift
dvoltage = dvoltage-avgv; %subtracting the difference from the voltage array
timeint = dtime-.6134; %setting the time such that the peturbation starts at t = 0 s

[~,~,maxvoltage,~] = peakdet(dvoltage,0.002);
Td = diff(timeint(maxvoltage));
Tdselect = Td(1); %picking first value

%log decrement method
n = 2; %two peak separation
sigman = (1/(n-1))*log(.0144/.003298); %peak numbers picked from visual inspection
zetan = sigman/sqrt((4*pi^2)+(sigman^2));

wdn = 2*pi/Tdselect; %rad/s, finding damped natural frequency
wnn = wdn/sqrt(1-(zetan^2)); %rad/s, converting to undamped with zetan

mn = 2/32.17; %slugs, mass of the plate
kn = ((wnn^2)*mn)/12; %lb/in, spring constant
K = 1/(kn*16); %gain K
bn = 2*zetan/wnn*kn; %damping coefficient, should be lb*s/in

system2 = tf(K, [1/(wnn^2) 2*zetan/wnn 1]);

figure(7)
impulse(system2)
hold on
plot(timeint,dvoltage/1000/sensitivitytwo,'g')
legend('Simulated System','Experimental Parameters')

figure(8)
plot(timeint,dvoltage)
title('test')
xlabel('Time (s)')
ylabel('Voltage (V)')
grid on

%% Part 3 Vibration Analysis
load lab4data.mat
acceltime = LVTaccel(:,1);
acceltime = acceltime+.04929;
accel = LVTaccel(:,2);
LVT = LVTaccel(:,4); %measures velocity

gvolts = -.1086; %V, from gravity, the acceleration during freefall, picked from a stable point on the plot
accelerometersensitivity = gvolts/386.09; %V/(in/s^2), the sensitivity of the acceleromter

%using method suggested by lab partner
startindex = 0;
endindex = 0;
for i = 1:length(acceltime)
    if (acceltime(i)>0.029 && endindex==0)
        endindex = i;
        break;
    elseif (acceltime(i)>0.0064 && startindex==0)
        startindex = i;
    end
end

timeint = acceltime(endindex)-acceltime(startindex);
dvelocity = timeint*386.09; %in/s
LVTinterval = LVT(endindex)-LVT(startindex);
LVTsensitivity = abs(LVTinterval/dvelocity); %LVT sensitivity V/(in/s)

figure(9)
plot(acceltime,accel,acceltime,LVT)%,acceltime(startindex),accel(startindex),'o',acceltime(endindex),accel(endindex),'o')
legend('Accelerometer','Velocity (LVT)')
title('LVT and Accelerometer Voltage Response')
ylabel('Voltage (V)')
xlabel('Time (s)')
grid on

%part b
integratedacceleration = cumtrapz(acceltime,(accel/accelerometersensitivity));

figure(10)
plot(acceltime,integratedacceleration,acceltime,LVT/LVTsensitivity)
title('Unadjusted Integrated Acceleration vs. LVT Response')
ylabel('Velocity (in/s)')
xlabel('Time (s)')
grid on

%normalizationslope = (integratedacceleration(end)-integratedacceleration(1))/(acceltime(end)-acceltime(1)); %slope of the skewed curve
%Normalizing
deltaY = integratedacceleration(end)-integratedacceleration(1); %finding the total vertical skew
biasRange = linspace(0,deltaY,length(integratedacceleration)); %creating an array of deltaYs
biasRange = biasRange'; %turning row vector into column vector
normia = integratedacceleration-biasRange; %subtracting bias from the original data

figure(11)
plot(acceltime,normia,acceltime,LVT/LVTsensitivity)
title('Adjusted Integrated Acceleration vs. LVT Response')
ylabel('Velocity (in/s)')
xlabel('Time (s)')
grid on

%part c

figure(12)
plot(acceltime,accel,acceltime,LVT,[0 0],[-3 3],'--',[0.068 0.068],[-3 3],'--',[0.087 0.087],[-3 3],'--',[0.12 0.12],[-3 3],'--',[0.57 0.57],[-3 3],'--') 
title('Annotated Response')
xlabel('Time (s)')
ylabel('LVT and Accelerometer Output (V)')
legend('Accelerometer','Velocity (LVT)','Beginning of Drop','First Contact','Max Depth','Bouncing Off','Resting on Foam')

%part d
%calculate maximum core velocity
corevelocity = abs(LVT/LVTsensitivity);
maxcorev = max(corevelocity); %in/s

%part e

forcesensitivity = 0.491; %V/lb

[~,ind] = max(abs(ch1)); %finding the index of maximum value on ch1 (force, but in volts) 
fmax = ((abs(ch0(ind)-ch0(1)))/forcesensitivity); %lbs, maximum force in array

fsteady = (abs(ch0(end)-ch0(1)))/forcesensitivity; %lbs, force at rest/steady state
masspd = fsteady/32.17; %slugs, mass of the core

figure(13)
plot(ftime,ch0,ftime,ch1)
title('LVT and Force Sensor Response')
xlabel('Time (s)')
ylabel('LVT and Force Sensor Response (V)')
grid on
legend('Force Sensor','LVT')











