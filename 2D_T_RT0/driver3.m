function [T] = driver3(trials)
%
%test program that takes in number of trials
%and outputs time
Time1 = zeros(trials,1);
Time2 = zeros(trials,1);
%Time3 = zeros(trials,1);
count = 0;
nxy = 1;
Elements = zeros(trials,1);
DegreesOfFreedom = zeros(trials,1);
while (count < trials) 
    count = count+1;
    [Times,DegreesOfFreedom(count)] = driver2(nxy,nxy);
    Time1(count) = Times(1);
    Time2(count) = Times(2);
    Elements(count) = 2*nxy * nxy;
    nxy = nxy * 2;
end
T = table(Elements, DegreesOfFreedom, Time1,Time2);
return %end of function