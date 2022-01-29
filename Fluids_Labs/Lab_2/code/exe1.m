clc
close all
clear all

time_stab=zeros(1,16);
input=[0.00,-0.314,-0.411,-0.695,-0.891,-1.186,-1.676,-2.656,...
       -5.598,-10.501,0.126,1.097,2.078,4.039,8.942,13.845]';
files = dir("calibr*.mat");


%% COMPUTE OUTPUT BY EXPERIMENTAL DATA
output=zeros(16,1);
min_delta= zeros(1,16);

for k=1:16
    
load(files(k).name)
n=length(V);
a=1:n;
res=cumsum(V)./(1:n);
% figure();
% subplot(2,1,1)
% plot(a/100,res);
% xlabel('time [s]');
% ylabel('cumulative voltage [v]');
% hold on
% z=[1,n/100];
% plot(z,[res(end),res(end)]);
variance= zeros(1,n);
find=0;
temp=k;
for i=1:n
variance(i)= var(res(i:end));
if (variance(i)/variance(1)<1.5e-2 && find==0)
    time_stab(k)= i/100;
    fprintf('Estimate of time acquisition for calibr %.1d with stopping criteria relative variance= %.2d is:\n',temp,1.5e-2);
    fprintf('%.2d seconds.\n',time_stab(k));
    find=1;
end  

min_delta(i)= abs( min( V(1:end-1)-V(2:end)));
end
% subplot(2,1,2)
% plot(a,variance, 'Linewidth', 1)

output(temp)=res(end);
end

%% LINEAR REGRESSION
X=[ones(length(input),1) input];
b1 = X\output;
my_outp = X*b1;
scatter(input,output)
hold on
plot(input,my_outp)
xlabel('force [N]')
ylabel('voltage [V]')
title('linear regression')

fprintf('----------------------\n\n')
fprintf('intercept: %.2d\n', b1(1));
fprintf('slope: %.2d\n\n', b1(2));

%% FIND UNCERTAINTY

% method 1
abs_diff= abs(((output-b1(1))/b1(2))-input);
U1= max(abs_diff);
fprintf('----------------------\n\n')
fprintf('for method 1 uncertainty U1: %.2d\n\n', U1);

% method 2
U2=2*sqrt(sum(((output-b1(1))/b1(2)-input).^2)/14);    
fprintf('----------------------\n\n')
fprintf('for method 2 uncertainty U2: %.2d\n',U2); 

%% RESOLUTION
resolution_volt = min(min_delta);
resolution_force = abs( resolution_volt/b1(2));

fprintf('----------------------\n\n')
fprintf('resolution in Volt: %.2d\n',resolution_volt); 
fprintf('resolution in Newton: %.2d\n',resolution_force); 