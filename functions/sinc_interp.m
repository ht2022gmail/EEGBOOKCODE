%% Whittaker–Shannon interpolation Algorithm
% https://jp.mathworks.com/matlabcentral/answers/2038371-sinc-interpolation-without-using-predefined-function
function y = sinc_interp(t, T, ts, xn)
    y      = zeros(length(ts), 1);
    for n  = 1:length(ts)   
        y  = y + xn(n)*sinc((t - (n-1)*T)/T);
    end
end