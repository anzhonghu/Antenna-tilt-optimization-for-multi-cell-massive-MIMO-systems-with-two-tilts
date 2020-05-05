%固定功率和天线数，当小区半径改变时，储存对应每个半径，容量取最大值时的两个倾角，图示为对应每个最大容量时的倾角随小区半径的变化
clear all;
close all;
h = 30;
r2 = 10;
M1 = 50;
M2 = 50;
M = 100;
K = 10;
N=20;
SLLaz = 25;
SLLel = 20;
SLLtot = 25;
alphamax = 0;
phi3dB = (70*pi)/180;
theta3dB = (7*pi)/180;
L=3;
Ni = 1e2;
P=10e+9;
Qs = (1:10);
Qss=20*Qs;
R1=zeros(N,N);
R2=zeros(1,10);
R3=zeros(1,10);

for p=1:length(Qs)
    r1=20*p;
    th1 = pi/2 - atan(r1/h);
    th2 = pi/2 - atan(r2/h);
    dthe = th2 - th1;
    for j = 1:N
        theta1 = ((j-1)/N)*(pi/2);
        for q = 1:N
            if q<j
                continue;
            end
            theta2 = ((q-1)/N)*(pi/2);
            sum1=0;
            for m = 1:M
                phi = ((m-1)/M)*(2*pi/3)-(pi/3);
                for n = 1:M
                    theta3 = th1 + ((n-1)/M)*(dthe);
                    d=h/sin(theta3);
                    alphatheta1 = -min(min(12*((phi)/phi3dB)^2, SLLaz)...
                        + min(12*((theta3 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                    alphatheta2 = -min(min(12*((phi)/phi3dB)^2, SLLaz)...
                        + min(12*((theta3 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                    alphatheta1=10^(alphatheta1/10)*d^(-1.5);
                    alphatheta2=10^(alphatheta2/10)*d^(-1.5);
                    sum1 = sum1 + K*((1 /(M1*(alphatheta1^2) + M2*(alphatheta2^2)))*(3/(2*pi))...
                        *((2*(h^2))/(r1^2-r2^2))*cot(theta3)*(1+(cot(theta3))^2)*((2*pi)/(3*M))*(dthe/M));
                end
            end
            R=K*log2(1+(P/sum1));
            R1(j,q)=R;
        end
    end
    [R1a, R1b] = find(R1==max(max(R1)));
    R2(1,p)=((R1a-1)/N)*(pi/2);
    R3(1,p)=((R1b-1)/N)*(pi/2);
end

plot(Qss, R2,'rs-'),xlabel('Cell radius r1[m]'),ylabel('Tilt angle[rad]');
hold on;
plot(Qss, R3,'gd-'),xlabel('Cell radius r1[m]'),ylabel('Tilt angle[rad]');
grid on;
legend('\theta_1','\theta_2','Location','NorthEast');