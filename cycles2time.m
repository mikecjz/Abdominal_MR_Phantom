function [respCat,numInCat] = cycles2time(cyclescat)
%This function expands cycle catagories by number of time points in each
%resporatory cycles.

%Input:
% cyclescat [1 x n]: array of resporatory cycles types [1,2,3]
% 
%Output:
% respCat   [1 x totalTP]: array of expanded resporatory cycles types.
% totalTP: total timeoints
%
% numInCat  [1 x totalTP]: order of timepoints in their current cycle.

respCat = [];
numInCat = [];

for i = 1:length(cyclescat)
   switch  cyclescat(i)
       case 1 % normal
           newTPs = cyclescat(i)*ones(1,70);
           respCat = [respCat,newTPs];
           order = 1:70;
           numInCat = [numInCat,order];
       case 2 % long
           newTPs = cyclescat(i)*ones(1,105);
           respCat = [respCat,newTPs];
           order = 1:105;
           numInCat = [numInCat,order];
       case 3 % short
           newTPs = cyclescat(i)*ones(1,53);
           respCat = [respCat,newTPs];
           order = 1:53;
           numInCat = [numInCat,order];
       otherwise
           error('Invalid elements among imputs. valid eliments: {1,2,3}')
   end
           
end

    
end