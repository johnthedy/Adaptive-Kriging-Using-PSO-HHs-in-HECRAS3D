function [mu,sigma,nRV,dist]=problem(probtype)
% Mean, standard deviation, and random variable distribution are defined in
% this matlab file

% In this example case, there are four ramdon variables which is rainfall
% intensity, manning coefficient, skewness and kurtosis of pearson pdf.
% skewness and kurtosis of pearson pdf is used to defined rainfall
% distribution along raining duration.

if probtype==1
    mu=[419.9 0.05 0.30868 2.1594];sigma=0.1*mu;nRV=4;dist={'Normal' 'Normal' 'Normal' 'Normal'};
end
end