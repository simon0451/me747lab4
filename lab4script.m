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
    if smdvolt(i) > vtau %finding the point where voltage exceeds vtau, working backwards
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







