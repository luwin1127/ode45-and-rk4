clc;clear;close all;

%% ������
% ʱ�����
t0 = 0;
tf = 100000;
% ��ֵ
y10 = 1;
y20 = 1;
y30 = 3;
y = [y10,y20,y30];
% RK �㷨����
h = 0.25;
% �Ը�΢�ַ�������ode45���Ա����������������б�s�ϣ��������£�
% ode45()����
tic;
[T,F] = ode45(@fun,[t0 tf],y);
time_record_ode = toc;
toc;

% �Ա�RK����
tic;
[T1,F1]=runge_kutta1(@fun,y,h,t0,tf);
time_record_rk = toc;
toc;

%% ��ͼ
figure('color',[1 1 1],'position',[600,200,500*1.2,416*1.2]);
subplot(121);
plot(T,F,'LineWidth',1.5);
title(['ode45','($t_f$=',num2str(tf),')'],'Interpreter','Latex');
set(gca,'FontSize',15,'FontName','Times New Roman','LineWidth',1.5);

subplot(122);
plot(T1,F1,'LineWidth',1.5);
title(['RK4','($t_f$=',num2str(tf),'),($h$=',num2str(h),')'],'Interpreter','Latex');
set(gca,'FontSize',15,'FontName','Times New Roman','LineWidth',1.5);

% ��������
str = ['data/tf_',num2str(tf),'_h_',num2str(h),'.mat'];
str_fig = ['img/tf_',num2str(tf),'_h_',num2str(h),'.jpg'];
save(str);
saveas(gcf,str_fig)
%% �Ӻ�������
% ΢�ַ���
function dy = fun(t,y)
dy = zeros(3,1);%��ʼ��������
dy(1) = y(2) * y(3);
dy(2) = -y(1) + y(3);
dy(3) = -0.51 * y(1) * y(2);
end

%--- RK �㷨 ---%
% ufunc - ΢�ַ�����ĺ�����
% y0 - ��ʼֵ
% h - ����
% a - ��ʼʱ��
% b - ĩ��ʱ��
function [x,y]=runge_kutta1(ufunc,y0,h,a,b)
% ����
n = floor((b-a)/h);
% ʱ�����
x = zeros(n,1);
x(1) = a;
% ��ʼֵ
y(:,1)=y0;

% �㷨��ʼ
for iter=1:n
% ʱ�����
x(iter+1) = x(iter)+h;
% ��ʼ����
k1 = ufunc(x(iter),       y(:,iter));
k2 = ufunc(x(iter) + h/2, y(:,iter) + h*k1 / 2);
k3 = ufunc(x(iter) + h/2, y(:,iter) + h*k2 / 2);
k4 = ufunc(x(iter) + h,   y(:,iter) + h*k3 );
% �õ����
y(:,iter+1) = y(:,iter) + h * (k1 + 2*k2 + 2*k3 + k4) / 6;
end

end