function cat = generatecycletype(numcycle)
% input:
% numcycle   [1 x 1] : number breathing cycles
% 
% output:
% cat        [1 x numcycle] : cycle catagory(1: normal 2: long 3: short)


rng(11);
randNum = rand(1,numcycle)*15;
cat = zeros(size(randNum));
cat(randNum<=1) = 2;
cat(randNum>1&randNum<=3) = 3;
cat(randNum>3) = 1;

end