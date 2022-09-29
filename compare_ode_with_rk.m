clc;clear;close all;

%% 主函数
% 时间变量
t0 = 0;
tf = 100000;
% 初值
y10 = 1;
y20 = 1;
y30 = 3;
y = [y10,y20,y30];
% RK 算法步长
h = 0.25;
% 对该微分方程组用ode45和自编的龙格库塔函数进行比s较，调用如下：
% ode45()函数
tic;
[T,F] = ode45(@fun,[t0 tf],y);
time_record_ode = toc;
toc;

% 自编RK函数
tic;
[T1,F1]=runge_kutta1(@fun,y,h,t0,tf);
time_record_rk = toc;
toc;

%% 画图
figure('color',[1 1 1],'position',[600,200,500*1.2,416*1.2]);
subplot(121);
plot(T,F,'LineWidth',1.5);
title(['ode45','($t_f$=',num2str(tf),')'],'Interpreter','Latex');
set(gca,'FontSize',15,'FontName','Times New Roman','LineWidth',1.5);

subplot(122);
plot(T1,F1,'LineWidth',1.5);
title(['RK4','($t_f$=',num2str(tf),'),($h$=',num2str(h),')'],'Interpreter','Latex');
set(gca,'FontSize',15,'FontName','Times New Roman','LineWidth',1.5);

% 保存数据
str = ['data/tf_',num2str(tf),'_h_',num2str(h),'.mat'];
str_fig = ['img/tf_',num2str(tf),'_h_',num2str(h),'.jpg'];
save(str);
saveas(gcf,str_fig)
%% 子函数部分
% 微分方程
function dy = fun(t,y)
dy = zeros(3,1);%初始化列向量
dy(1) = y(2) * y(3);
dy(2) = -y(1) + y(3);
dy(3) = -0.51 * y(1) * y(2);
end

%--- RK 算法 ---%
% ufunc - 微分方程组的函数名
% y0 - 初始值
% h - 步长
% a - 初始时间
% b - 末端时间
function [x,y]=runge_kutta1(ufunc,y0,h,a,b)
% 步数
n = floor((b-a)/h);
% 时间起点
x = zeros(n,1);
x(1) = a;
% 初始值
y(:,1)=y0;

% 算法开始
for iter=1:n
% 时间变量
x(iter+1) = x(iter)+h;
% 开始迭代
k1 = ufunc(x(iter),       y(:,iter));
k2 = ufunc(x(iter) + h/2, y(:,iter) + h*k1 / 2);
k3 = ufunc(x(iter) + h/2, y(:,iter) + h*k2 / 2);
k4 = ufunc(x(iter) + h,   y(:,iter) + h*k3 );
% 得到结果
y(:,iter+1) = y(:,iter) + h * (k1 + 2*k2 + 2*k3 + k4) / 6;
end

end