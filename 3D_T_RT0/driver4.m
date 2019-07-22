function [T] = driver4(trials)

%test program that takes in number of trials
%and outputs time
%Includes non-vectorized
Time1 = zeros(trials,1); %Standard 1
Time2 = zeros(trials,1); %Standard assembly
Time3 = zeros(trials,1); %Vectorized 1
Time4 = zeros(trials,1); %Vectorized assembly
count = 0;
nxyz = 4;
Elements = zeros(trials,1);
DegreesOfFreedom = zeros(trials,1);
while (count < trials) 
    count = count+1;
    [Times,DegreesOfFreedom(count)] = driver(nxyz,nxyz,nxyz);
    Time1(count) = Times(1);
    Time2(count) = Times(2);
    Time3(count) = Times(3);
    Time4(count) = Times(4);
    Elements(count) = 6*nxyz * nxyz*nxyz;
    nxyz = nxyz * 2;
end
T = table(Elements, DegreesOfFreedom, Time1, Time2, Time3, Time4);
return %end of function