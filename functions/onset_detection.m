function onset = onset_detection(x, type)

if nargin<2
    type = 1;
end

x_shift = [0 x(1:end-1)];
if nargin<2
    onset = (x>0 & x_shift==0);   
else
    onset = (x==type & x_shift==0);
end
