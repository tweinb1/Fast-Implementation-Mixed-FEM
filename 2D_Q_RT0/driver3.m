function [T] = driver3(trials)
% Test program that takes in number of trials
% and outputs time.
% The starting number of elements is 16, 
% and is quadrupled each time.
% This returns a table providing data used in the paper, 
% but only runs the vectorized code.

Stiffness = zeros(trials,1);
Assembly = zeros(trials,1);
count = 0;
nxy = 4;
Elements = zeros(trials,1);
DegreesOfFreedom = zeros(trials,1);
while (count < trials) 
    count = count+1;
    [Times,DegreesOfFreedom(count)] = driver2(nxy,nxy);
    Stiffness(count) = Times(1);
    Assembly(count) = Times(2);
    Elements(count) = nxy * nxy;
    nxy = nxy * 2;
end
T = table(Elements, DegreesOfFreedom, Stiffness, Assembly);
return %end of function