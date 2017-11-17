%Simon Popecki
%ME 747
%13 Novermber 2017
%Lab 4
clear all;
close all;
load lab4data.mat

%% Part 1 Potentiometer Accelerometer
mass_nobucket = [20 40 60 80 100 120 140 160 180 200]; %grams
bucketWeight = 15.9; %grams
testmass = mass_nobucket+bucketWeight; %testmass (g) is the total mass on the sensor
testmassoz = testmass*0.035274; %converting from grams to ounces because Fussell likes imperial
eoinc = [501 527 554 582 607 633 658 684 712 737]; %mV, output voltage with weight increasing
eodec = [503 528 552 580 607 633  659 686 711 738]; %mV, output voltage with weight decreasing
eoavg = (eoinc+eodec)/2; %averageing the output voltages of increasing and decreasing values
bestfit = polyfit(testmassoz,eoavg,1); %making a best fit line to find sensitivity

sensitivity = bestfit(1)/1000; %the sensitivity of the scale in volts per ounce

% Plotting eo vs. weight
figure (1)
plot(testmassoz,eoinc,testmassoz,eodec)
title('Potentiometer Accelerometer Proof Weight vs. Output Voltage')
xlabel('Proof Weight (Oz)')
ylabel('Output Voltage (mV)')
legend('Increasing Weight','Decreasing Weight','location','northwest')
grid on

% Plotting eo vs. time
figure (2)
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
m = ((.454-.392)/sensitivity)*0.00194256; %slugs, effective mass
k = (16/12)*(omegan^2)*m; %oz/in, the spring constant - !highly dependant on peak time (16/12 for converison from slugs to Oz/in)













