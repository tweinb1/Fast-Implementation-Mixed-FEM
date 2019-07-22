function [T] = driver4(trials)
% Test program that takes in number of trials
% and outputs time.
% The starting number of elements is 16, 
% and is quadrupled each time.
% This returns a table providing data used in the paper, 
% and runs both vectorized and non-vectorized

Time1 = zeros(trials,1); %Standard 1
Time2 = zeros(trials,1); %Standard assembly
Time3 = zeros(trials,1); %Vectorized 1
Time4 = zeros(trials,1); %Vectorized assembly
count = 0;
nxy = 4;
Elements = zeros(trials,1);
DegreesOfFreedom = zeros(trials,1);
while (count < trials) 
    count = count+1;
    [Times,DegreesOfFreedom(count)] = driver(nxy,nxy);
    Time1(count) = Times(1);
    Time2(count) = Times(2);
    Time3(count) = Times(3);
    Time4(count) = Times(4);
    Elements(count) = 2*nxy * nxy;
    nxy = nxy * 2;
end
T = table(Elements, DegreesOfFreedom, Time1, Time2, Time3, Time4);
return %end of function