%% Simulation of dixon2003
clear; close all; clc;
%% Time interval and simulation time
Step = 0.0005;T_end = 5;
t = 0:Step:T_end;
%% Variable
x1=cell(1,size(t,2));
x2=cell(1,size(t,2));
x3=cell(1,size(t,2));
y1=cell(1,size(t,2));
y2=cell(1,size(t,2));
fm_1=cell(1,size(t,2));
fm_2=cell(1,size(t,2));
ym_1=cell(1,size(t,2));
ym_2=cell(1,size(t,2));
%% Parameter
k_s1=30;
alp_1=8;
gam_1=6;
k_s2=40;
alp_2=5;
gam_2=1.5;
A=[-0.2 0.4 -0.6;0.1 -0.2 0.3;0.3 -0.4 0.4];
B=[0.5;0.25;0.3];
%% Initial value
x1{1}=1;
x2{1}=1.5;
x3{1}=2.5;
fm_1{1}=1;
fm_2{1}=1;
ym_1{1}=x1{1}/x3{1};
ym_2{1}=x2{1}/x3{1};
%% Simulation
for i=1:size(t,2)
    y1{i}=x1{i}/x3{i};
    y2{i}=x2{i}/x3{i};
    e1=y1{i}-ym_1{i};
    e2=y2{i}-ym_2{i};
    dx=A*[x1{i};x2{i};x3{i}]+B;
    if i==size(t,2)
        break
    end
    %% Update state
    fm_1{i+1}=fm_1{i}+Step*(-(k_s1+alp_1)*fm_1{i}+gam_1*satlins(e1/0.001)+alp_1*k_s1*e1);
    fm_2{i+1}=fm_2{i}+Step*(-(k_s2+alp_2)*fm_2{i}+gam_2*satlins(e2/0.001)+alp_2*k_s2*e2);
    ym_1{i+1}=ym_1{i}+Step*(A(1,3)+(A(1,1)-A(3,3))*y1{i}+A(1,2)*y2{i}-A(3,1)*y1{i}^2-A(3,2)*y1{i}*y2{i}+fm_1{i});
    ym_2{i+1}=ym_2{i}+Step*(A(2,3)+A(2,1)*y1{i}+(A(2,2)-A(3,3))*y2{i}-A(3,2)*y2{i}^2-A(3,1)*y1{i}*y2{i}+fm_2{i});
    x1{i+1}=x1{i}+Step*dx(1);
    x2{i+1}=x2{i}+Step*dx(2);
    x3{i+1}=x3{i}+Step*dx(3);
end

y1_arr=cell2mat(y1);
y2_arr=cell2mat(y2);
ym_1_arr=cell2mat(ym_1);
ym_2_arr=cell2mat(ym_2);
fm_1_arr=cell2mat(fm_1);
fm_2_arr=cell2mat(fm_2);
x3_arr=cell2mat(x3);
f1_arr=(B(1)-B(3)*y1_arr)./x3_arr;
f2_arr=(B(2)-B(3)*y2_arr)./x3_arr;
figure(1);
plot(t,f1_arr-fm_1_arr);
figure(2);
plot(t,f2_arr-fm_2_arr);
figure(3);
xm_3_arr=1./sqrt((fm_1_arr.^2+fm_2_arr.^2)./((B(1)-B(3)*y1_arr).^2+(B(2)-B(3)*y2_arr).^2));
plot(t,x3_arr-xm_3_arr);