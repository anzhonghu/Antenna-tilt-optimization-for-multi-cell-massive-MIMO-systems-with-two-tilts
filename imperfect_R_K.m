%固定功率和倾角，容量随用户数K变化曲线
clear all;
close all;
h = 30;
r1 = 120;
r2 = 10;
M=100;
M1=50;
M2=50;
SLLaz = 25;
SLLel = 20;
SLLtot = 25;
alphamax = 0;
phi3dB = (70*pi)/180;
theta3dB = (7*pi)/180;
L=3;
theta1=0.3142;
theta2=0.55;
theta3=0.34;     %倾角指向中心
theta4=0.34;
theta5=0.2808;     %倾角指向两个扇区
theta6=0.4571;
sumrate = zeros(1,5);
sumrate1 = zeros(1,5);
sumrate_1 = zeros(1,5);
sumrate_2 = zeros(1,5);
R = zeros(1,5);
th1 = pi/2 - atan(r1/h);
th2 = pi/2 - atan(r2/h);
dthe = th2 - th1;
Qs = (1:5);
Qss=5*Qs;
Ni = 1e2;
P=10e+9;

for i = 1:Ni
    rho11=0;
    rho22=0;
    rho33=0;
    for q=1:length(Qs)
        K=5*q;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%以第1个小区为基准
        k1 = 0;
        A_1 = zeros (K,2);
        H_11 = zeros(M,K);      %第1个小区基站到其用户的信道Hjj
        H_21 = zeros(M,K);      %第2个小区基站到第一个小区用户的信道Hj'j
        H_31 = zeros(M,K);      %第3个小区基站到第一个小区用户的信道Hj'j
        H_11_1 = zeros(M,K);
        H_11_2 = zeros(M,K);
        H_21_1 = zeros(M,K);
        H_21_2 = zeros(M,K);
        H_31_1 = zeros(M,K);
        H_31_2 = zeros(M,K);
        
        alphatheta1_1=zeros(1,K);
        alphatheta2_1=zeros(1,K);
        alphatheta1_2=zeros(1,K);
        alphatheta2_2=zeros(1,K);
        alphatheta1_3=zeros(1,K);
        alphatheta2_3=zeros(1,K);
        alphas = zeros(K,8);
        angles = zeros(Ni*K,18);
        while k1<K                                                             %第一个小区用户位置
            x1 = [-(sqrt(3)/2)*r1 (sqrt(3)/2)*r1];
            y1 = [0 -r1];
            x1 = (max(x1)-min(x1))*rand(1,1)+min(x1);
            y1 = (max(y1)-min(y1))*rand(1,1)+min(y1);
            
            if (r2<=sqrt(x1^2+y1^2)) && (sqrt(x1^2+y1^2)<=r1) && atan((abs(x1))/(abs(y1)))<=pi/3
                A_1(k1+1,1) = x1;
                A_1(k1+1,2) = y1;
                phi_11 = atan(x1/(abs(y1)));                                   %第1个小区基站到其用户的信道Hjj
                theta3_11 = (pi/2)-atan((sqrt(x1^2+y1^2))/h);
                d_11=sqrt(h^2+x1^2+y1^2);
                alphatheta1_11 = -min(min(12*(phi_11/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_11 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_11 = -min(min(12*(phi_11/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_11 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta1_11_1 = -min(min(12*(phi_11/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_11 - theta3)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_11_1 = -min(min(12*(phi_11/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_11 - theta4)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta1_11_2 = -min(min(12*(phi_11/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_11 - theta5)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_11_2 = -min(min(12*(phi_11/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_11 - theta6)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                
                alphatheta1_11=10^(alphatheta1_11/10)*d_11^(-1.5);
                alphatheta2_11=10^(alphatheta2_11/10)*d_11^(-1.5);
                alphatheta1_11_1=10^(alphatheta1_11_1/10)*d_11^(-1.5);
                alphatheta2_11_1=10^(alphatheta2_11_1/10)*d_11^(-1.5);
                alphatheta1_11_2=10^(alphatheta1_11_2/10)*d_11^(-1.5);
                alphatheta2_11_2=10^(alphatheta2_11_2/10)*d_11^(-1.5);
                
                alphatheta1_1(:,k1+1)=alphatheta1_11;
                alphatheta2_1(:,k1+1)=alphatheta2_11;
                h1_11=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                h2_11=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                h_11=[alphatheta1_11*h1_11;alphatheta2_11*h2_11];
                h_11_1=[alphatheta1_11_1*h1_11;alphatheta2_11_1*h2_11];
                h_11_2=[alphatheta1_11_2*h1_11;alphatheta2_11_2*h2_11];
                
                H_11(:,k1+1)=h_11;
                H_11_1(:,k1+1)=h_11_1;
                H_11_2(:,k1+1)=h_11_2;
                angles((i-1)*K+k1+1, 1) = phi_11;
                angles((i-1)*K+k1+1, 2) = theta3_11;
                
                phi_21 = -((pi/3)-atan(abs(x1-r1)/(y1+(sqrt(3))*r1)));         %第2个小区基站到第一个小区用户的信道Hj'j
                theta3_21 = (pi/2)-atan((sqrt((x1-r1)^2+(y1+(sqrt(3))*r1)^2))/h);
                d_21=sqrt(h^2+(x1-r1)^2+(y1+(sqrt(3))*r1)^2);
                alphatheta1_21 = -min(min(12*(phi_21/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_21 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_21 = -min(min(12*(phi_21/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_21 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta1_21_1 = -min(min(12*(phi_21/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_21 - theta3)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_21_1 = -min(min(12*(phi_21/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_21 - theta4)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta1_21_2 = -min(min(12*(phi_21/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_21 - theta5)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_21_2 = -min(min(12*(phi_21/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_21 - theta6)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                
                alphatheta1_21=10^(alphatheta1_21/10)*d_21^(-1.5);
                alphatheta2_21=10^(alphatheta2_21/10)*d_21^(-1.5);
                alphatheta1_21_1=10^(alphatheta1_21_1/10)*d_21^(-1.5);
                alphatheta2_21_1=10^(alphatheta2_21_1/10)*d_21^(-1.5);
                alphatheta1_21_2=10^(alphatheta1_21_2/10)*d_21^(-1.5);
                alphatheta2_21_2=10^(alphatheta2_21_2/10)*d_21^(-1.5);
                
                alphatheta1_2(:,k1+1)=alphatheta1_21;
                alphatheta2_2(:,k1+1)=alphatheta2_21;
                h1_21=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                h2_21=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                h_21=[alphatheta1_21*h1_21;alphatheta2_21*h2_21];
                h_21_1=[alphatheta1_21_1*h1_21;alphatheta2_21_1*h2_21];
                h_21_2=[alphatheta1_21_2*h1_21;alphatheta2_21_2*h2_21];
                
                H_21(:,k1+1)=h_21;
                H_21_1(:,k1+1)=h_21_1;
                H_21_2(:,k1+1)=h_21_2;
                
                angles((i-1)*K+k1+1, 3) = phi_21;
                angles((i-1)*K+k1+1, 4) = theta3_21;
                alphas(k1+1,1) = alphatheta1_21;
                alphas(k1+1,2) = alphatheta2_21;
                
                phi_31 =atan((y1+(sqrt(3))*r1)/(x1+r1))-pi/6;                  %第3个小区基站到第一个小区用户的信道Hj'j
                theta3_31 = (pi/2)-atan((sqrt((x1+r1)^2+(y1+(sqrt(3))*r1)^2))/h);
                d_31=sqrt(h^2+(x1+r1)^2+(y1+(sqrt(3))*r1)^2);
                alphatheta1_31 = -min(min(12*(phi_31/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_31 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_31 = -min(min(12*(phi_31/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_31 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta1_31_1 = -min(min(12*(phi_31/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_31 - theta3)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_31_1 = -min(min(12*(phi_31/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_31 - theta4)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta1_31_2 = -min(min(12*(phi_31/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_31 - theta5)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_31_2 = -min(min(12*(phi_31/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_31 - theta6)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                
                alphatheta1_31=10^(alphatheta1_31/10)*d_31^(-1.5);
                alphatheta2_31=10^(alphatheta2_31/10)*d_31^(-1.5);
                alphatheta1_31_1=10^(alphatheta1_31_1/10)*d_31^(-1.5);
                alphatheta2_31_1=10^(alphatheta2_31_1/10)*d_31^(-1.5);
                alphatheta1_31_2=10^(alphatheta1_31_2/10)*d_31^(-1.5);
                alphatheta2_31_2=10^(alphatheta2_31_2/10)*d_31^(-1.5);
                
                alphatheta1_3(:,k1+1)=alphatheta1_31;
                alphatheta2_3(:,k1+1)=alphatheta2_31;
                h1_31=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                h2_31=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                h_31=[alphatheta1_31*h1_31;alphatheta2_31*h2_31];
                h_31_1=[alphatheta1_31_1*h1_31;alphatheta2_31_1*h2_31];
                h_31_2=[alphatheta1_31_2*h1_31;alphatheta2_31_2*h2_31];
                
                H_31(:,k1+1)=h_31;
                H_31_1(:,k1+1)=h_31_1;
                H_31_2(:,k1+1)=h_31_2;
                
                angles((i-1)*K+k1+1, 5) = phi_31;
                angles((i-1)*K+k1+1, 6) = theta3_31;
                alphas(k1+1,3) = alphatheta1_31;
                alphas(k1+1,4) = alphatheta2_31;
                
                k1=k1+1;
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%以第2个小区为基准
        k2 = 0;
        A_2 = zeros(K, 2);
        H_22 = zeros(M,K);      %第2个小区基站到其自己用户的信道Hjj
        H_12 = zeros(M,K);      %第1个小区基站到第2个小区用户的信道Hj'j
        H_32 = zeros(M,K);      %第3个小区基站到第2个小区用户的信道Hj'j
        H_22_1 = zeros(M,K);
        H_22_2 = zeros(M,K);
        H_12_1 = zeros(M,K);
        H_12_2 = zeros(M,K);
        H_32_1 = zeros(M,K);
        H_32_2 = zeros(M,K);
        
        alphatheta1_4=zeros(1,K);
        alphatheta2_4=zeros(1,K);
        alphatheta1_5=zeros(1,K);
        alphatheta2_5=zeros(1,K);
        alphatheta1_6=zeros(1,K);
        alphatheta2_6=zeros(1,K);
        while k2 < K                                                           %第2个小区用户位置
            x2 = [-r1 0];
            y2 = [-r1/2 r1];
            x2 = (max(x2)-min(x2))*rand(1,1)+min(x2);
            y2 = (max(y2)-min(y2))*rand(1,1)+min(y2);
            
            if (r2<=sqrt(x2^2+y2^2)) && (sqrt(x2^2+y2^2)<=r1) && ((y2<0 && atan(abs(y2/x2))<=pi/6) || y2>0)
                
                A_2(k2+1,1)=x2;
                A_2(k2+1,2)=y2;
                
                phi_22 = atan(y2/x2)+pi/6;                                     %第2个小区基站到其用户的信道Hj'j'
                theta3_22 = (pi/2)-atan((sqrt(x2^2+y2^2))/h);
                d_22=sqrt(h^2+x2^2+y2^2);
                alphatheta1_22= -min(min(12*(phi_22/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_22 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_22 = -min(min(12*(phi_22/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_22 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta1_22_1= -min(min(12*(phi_22/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_22 - theta3)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_22_1 = -min(min(12*(phi_22/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_22 - theta4)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta1_22_2= -min(min(12*(phi_22/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_22 - theta5)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_22_2 = -min(min(12*(phi_22/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_22 - theta6)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                
                alphatheta1_22=10^(alphatheta1_22/10)*d_22^(-1.5);
                alphatheta2_22=10^(alphatheta2_22/10)*d_22^(-1.5);
                alphatheta1_22_1=10^(alphatheta1_22_1/10)*d_22^(-1.5);
                alphatheta2_22_1=10^(alphatheta2_22_1/10)*d_22^(-1.5);
                alphatheta1_22_2=10^(alphatheta1_22_2/10)*d_22^(-1.5);
                alphatheta2_22_2=10^(alphatheta2_22_2/10)*d_22^(-1.5);
                
                alphatheta1_4(:,k2+1)=alphatheta1_22;
                alphatheta2_4(:,k2+1)=alphatheta2_22;
                h1_22=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                h2_22=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                
                h_22=[alphatheta1_22*h1_22;alphatheta2_22*h2_22];
                h_22_1=[alphatheta1_22_1*h1_22;alphatheta2_22_1*h2_22];
                h_22_2=[alphatheta1_22_2*h1_22;alphatheta2_22_2*h2_22];
                
                H_22(:,k2+1)=h_22;
                H_22_1(:,k2+1)=h_22_1;
                H_22_2(:,k2+1)=h_22_2;
                
                angles((i-1)*K+k2+1, 7) = phi_22;
                angles((i-1)*K+k2+1, 8) = theta3_22;
                alphas(k2+1,5) = alphatheta1_22;
                alphas(k2+1,6) = alphatheta2_22;
                
                
                phi_12 = atan((x2+r1)/(abs(y2-r1*sqrt(3))));                   %第1个小区基站到第2个小区用户的信道Hjj'
                theta3_12 = (pi/2)-atan((sqrt((x2+r1)^2+(y2-r1*sqrt(3))^2))/h);
                d_12=sqrt(h^2+(x2+r1)^2+(y2-r1*sqrt(3))^2);
                alphatheta1_12 = -min(min(12*(phi_12/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_12 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_12 = -min(min(12*(phi_12/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_12 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta1_12_1 = -min(min(12*(phi_12/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_12 - theta3)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_12_1 = -min(min(12*(phi_12/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_12 - theta4)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta1_12_2 = -min(min(12*(phi_12/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_12 - theta5)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_12_2 = -min(min(12*(phi_12/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_12 - theta6)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                
                alphatheta1_12=10^(alphatheta1_12/10)*d_12^(-1.5);
                alphatheta2_12=10^(alphatheta2_12/10)*d_12^(-1.5);
                alphatheta1_12_1=10^(alphatheta1_12_1/10)*d_12^(-1.5);
                alphatheta2_12_1=10^(alphatheta2_12_1/10)*d_12^(-1.5);
                alphatheta1_12_2=10^(alphatheta1_12_2/10)*d_12^(-1.5);
                alphatheta2_12_2=10^(alphatheta2_12_2/10)*d_12^(-1.5);
                
                alphatheta1_5(:,k2+1)=alphatheta1_12;
                alphatheta2_5(:,k2+1)=alphatheta2_12;
                h1_12=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                h2_12=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                h_12=[alphatheta1_12*h1_12;alphatheta2_12*h2_12];
                h_12_1=[alphatheta1_12_1*h1_12;alphatheta2_12_1*h2_12];
                h_12_2=[alphatheta1_12_2*h1_12;alphatheta2_12_2*h2_12];
                
                H_12(:,k2+1)=h_12;
                H_12_1(:,k2+1)=h_12_1;
                H_12_2(:,k2+1)=h_12_2;
                angles((i-1)*K+k2+1, 9) = phi_12;
                angles((i-1)*K+k2+1, 10) = theta3_12;
                
                phi_32 = atan(y2/(x2+2*r1))-pi/6;                              %第3个小区基站到第2个小区用户的信道Hj''j'
                theta3_32 = (pi/2)-atan((sqrt((x2+2*r1)^2+y2^2))/h);
                d_32=sqrt(h^2+(x2+2*r1)^2+y2^2);
                alphatheta1_32= -min(min(12*(phi_32/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_32 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_32 = -min(min(12*(phi_32/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_32 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta1_32_1= -min(min(12*(phi_32/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_32 - theta3)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_32_1 = -min(min(12*(phi_32/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_32 - theta4)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta1_32_2= -min(min(12*(phi_32/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_32 - theta5)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_32_2 = -min(min(12*(phi_32/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_32 - theta6)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                
                alphatheta1_32=10^(alphatheta1_32/10)*d_32^(-1.5);
                alphatheta2_32=10^(alphatheta2_32/10)*d_32^(-1.5);
                alphatheta1_32_1=10^(alphatheta1_32_1/10)*d_32^(-1.5);
                alphatheta2_32_1=10^(alphatheta2_32_1/10)*d_32^(-1.5);
                alphatheta1_32_2=10^(alphatheta1_32_2/10)*d_32^(-1.5);
                alphatheta2_32_2=10^(alphatheta2_32_2/10)*d_32^(-1.5);
                
                alphatheta1_6(:,k2+1)=alphatheta1_32;
                alphatheta2_6(:,k2+1)=alphatheta2_32;
                h1_32=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                h2_32=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                h_32=[alphatheta1_32*h1_32;alphatheta2_32*h2_32];
                h_32_1=[alphatheta1_32_1*h1_32;alphatheta2_32_1*h2_32];
                h_32_2=[alphatheta1_32_2*h1_32;alphatheta2_32_2*h2_32];
                
                H_32(:,k2+1)=h_32;
                H_32_1(:,k2+1)=h_32_1;
                H_32_2(:,k2+1)=h_32_2;
                angles((i-1)*K+k2+1, 11) = phi_32;
                angles((i-1)*K+k2+1, 12) = theta3_32;
                
                k2=k2+1;
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%以第3个小区为基准
        k3 = 0;
        A_3 = zeros(K, 2);
        H_33 = zeros(M,K);      %第3个小区基站到其自己用户的信道Hjj
        H_13 = zeros(M,K);      %第1个小区基站到第3个小区用户的信道Hj'j
        H_23 = zeros(M,K);      %第2个小区基站到第3个小区用户的信道Hj'j
        H_33_1 = zeros(M,K);
        H_33_2 = zeros(M,K);
        H_13_1 = zeros(M,K);
        H_13_2 = zeros(M,K);
        H_23_1 = zeros(M,K);
        H_23_2 = zeros(M,K);
        
        alphatheta1_7=zeros(1,K);
        alphatheta2_7=zeros(1,K);
        alphatheta1_8=zeros(1,K);
        alphatheta2_8=zeros(1,K);
        alphatheta1_9=zeros(1,K);
        alphatheta2_9=zeros(1,K);
        while k3 < K                                                           %第3个小区用户位置
            x3 = [0 r1];
            y3 = [-r1/2 r1];
            x3 = (max(x3)-min(x3))*rand(1,1)+min(x3);
            y3 = (max(y3)-min(y3))*rand(1,1)+min(y3);
            
            if (r2<=sqrt(x3^2+y3^2)) && (sqrt(x3^2+y3^2)<=r1) && ((y3<0 && atan(abs(y3/x3))<=pi/6) || y3>0)
                A_3(k3+1,1)=x3;
                A_3(k3+1,2)=y3;
                phi_33 = atan(y3/x3)-pi/6;                                     %第3个小区基站到其用户的信道Hj'j'
                theta3_33 = (pi/2)-atan((sqrt(x3^2+y3^2))/h);
                d_33=sqrt(h^2+x3^2+y3^2);
                alphatheta1_33= -min(min(12*(phi_33/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_33 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_33 = -min(min(12*(phi_33/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_33 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta1_33_1= -min(min(12*(phi_33/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_33 - theta3)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_33_1 = -min(min(12*(phi_33/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_33 - theta4)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta1_33_2= -min(min(12*(phi_33/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_33 - theta5)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_33_2 = -min(min(12*(phi_33/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_33 - theta6)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                
                alphatheta1_33=10^(alphatheta1_33/10)*d_33^(-1.5);
                alphatheta2_33=10^(alphatheta2_33/10)*d_33^(-1.5);
                alphatheta1_33_1=10^(alphatheta1_33_1/10)*d_33^(-1.5);
                alphatheta2_33_1=10^(alphatheta2_33_1/10)*d_33^(-1.5);
                alphatheta1_33_2=10^(alphatheta1_33_2/10)*d_33^(-1.5);
                alphatheta2_33_2=10^(alphatheta2_33_2/10)*d_33^(-1.5);
                
                alphatheta1_7(:,k3+1)=alphatheta1_33;
                alphatheta2_7(:,k3+1)=alphatheta2_33;
                h1_33=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                h2_33=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                h_33=[alphatheta1_33*h1_33;alphatheta2_33*h2_33];
                h_33_1=[alphatheta1_33_1*h1_33;alphatheta2_33_1*h2_33];
                h_33_2=[alphatheta1_33_2*h1_33;alphatheta2_33_2*h2_33];
                
                H_33(:,k3+1)=h_33;
                H_33_1(:,k3+1)=h_33_1;
                H_33_2(:,k3+1)=h_33_2;
                angles((i-1)*K+k3+1, 13) = phi_33;
                angles((i-1)*K+k3+1, 14) = theta3_33;
                alphas(k3+1,7) = alphatheta1_33;
                alphas(k3+1,8) = alphatheta2_33;
                
                phi_13 = -atan(abs(x3-r1)/(abs(y3-r1*sqrt(3))));               %第1个小区基站到第3个小区用户的信道Hjj'
                theta3_13 = (pi/2)-atan((sqrt((x3-r1)^2+(y3-r1*sqrt(3))^2))/h);
                d_13=sqrt(h^2+(x3-r1)^2+(y3-r1*sqrt(3))^2);
                alphatheta1_13 = -min(min(12*(phi_13/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_13 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_13 = -min(min(12*(phi_13/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_13 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta1_13_1 = -min(min(12*(phi_13/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_13 - theta3)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_13_1 = -min(min(12*(phi_13/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_13 - theta4)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta1_13_2 = -min(min(12*(phi_13/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_13 - theta5)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_13_2 = -min(min(12*(phi_13/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_13 - theta6)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                
                alphatheta1_13=10^(alphatheta1_13/10)*d_13^(-1.5);
                alphatheta2_13=10^(alphatheta2_13/10)*d_13^(-1.5);
                alphatheta1_13_1=10^(alphatheta1_13_1/10)*d_13^(-1.5);
                alphatheta2_13_1=10^(alphatheta2_13_1/10)*d_13^(-1.5);
                alphatheta1_13_2=10^(alphatheta1_13_2/10)*d_13^(-1.5);
                alphatheta2_13_2=10^(alphatheta2_13_2/10)*d_13^(-1.5);
                
                alphatheta1_8(:,k3+1)=alphatheta1_13;
                alphatheta2_8(:,k3+1)=alphatheta2_13;
                h1_13=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                h2_13=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                h_13=[alphatheta1_13*h1_13;alphatheta2_13*h2_13];
                h_13_1=[alphatheta1_13_1*h1_13;alphatheta2_13_1*h2_13];
                h_13_2=[alphatheta1_13_2*h1_13;alphatheta2_13_2*h2_13];
                
                H_13(:,k3+1)=h_13;
                H_13_1(:,k3+1)=h_13_1;
                H_13_2(:,k3+1)=h_13_2;
                angles((i-1)*K+k3+1, 15) = phi_13;
                angles((i-1)*K+k3+1, 16) = theta3_13;
                
                phi_23 = pi/6+atan((y3/(x3-2*r1)));                            %第2个小区基站到第3个小区用户的信道Hj''j'
                theta3_23 = (pi/2)-atan((sqrt((x3-2*r1)^2+y3^2))/h);
                d_23=sqrt(h^2+(x3-2*r1)^2+y3^2);
                alphatheta1_23= -min(min(12*(phi_23/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_23 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_23 = -min(min(12*(phi_23/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_23 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta1_23_1= -min(min(12*(phi_23/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_23 - theta3)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_23_1= -min(min(12*(phi_23/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_23 - theta4)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta1_23_2= -min(min(12*(phi_23/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_23 - theta5)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                alphatheta2_23_2= -min(min(12*(phi_23/phi3dB)^2, SLLaz)...
                    + min(12*((theta3_23 - theta6)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                
                alphatheta1_23=10^(alphatheta1_23/10)*d_23^(-1.5);
                alphatheta2_23=10^(alphatheta2_23/10)*d_23^(-1.5);
                alphatheta1_23_1=10^(alphatheta1_23_1/10)*d_23^(-1.5);
                alphatheta2_23_1=10^(alphatheta2_23_1/10)*d_23^(-1.5);
                alphatheta1_23_2=10^(alphatheta1_23_2/10)*d_23^(-1.5);
                alphatheta2_23_2=10^(alphatheta2_23_2/10)*d_23^(-1.5);
                
                alphatheta1_9(:,k3+1)=alphatheta1_23;
                alphatheta2_9(:,k3+1)=alphatheta2_23;
                h1_23=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                h2_23=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                h_23=[alphatheta1_23*h1_23;alphatheta2_23*h2_23];
                h_23_1=[alphatheta1_23_1*h1_23;alphatheta2_23_1*h2_23];
                h_23_2=[alphatheta1_23_2*h1_23;alphatheta2_23_2*h2_23];
                
                H_23(:,k3+1)=h_23;
                H_23_1(:,k3+1)=h_23_1;
                H_23_2(:,k3+1)=h_23_2;
                angles((i-1)*K+k3+1, 17) = phi_23;
                angles((i-1)*K+k3+1, 18) = theta3_23;
                k3=k3+1;
            end
        end
        
        H1 = H_11 + H_12 + H_13;      %Hjj
        H2 = H_22 + H_21 + H_23;      %Hj'j'
        H3 = H_33 + H_31 + H_32;      %Hj'j'
        H1_1 = H_11_1 + H_12_1 + H_13_1;      %Hjj
        H2_1 = H_22_1 + H_21_1 + H_23_1;      %Hj'j'
        H3_1 = H_33_1 + H_31_1 + H_32_1;      %Hj'j'
        H1_2 = H_11_2 + H_12_2 + H_13_2;      %Hjj
        H2_2 = H_22_2 + H_21_2 + H_23_2;      %Hj'j'
        H3_2 = H_33_2 + H_31_2 + H_32_2;      %Hj'j'
        
        F1 = H1*(H1'*H1)^(-1);     %Fj
        F2 = H2*(H2'*H2)^(-1);     %Fj'
        F3 = H3*(H3'*H3)^(-1);     %Fj'
        F1_1 = H1_1*(H1_1'*H1_1)^(-1);     %Fj
        F2_1 = H2_1*(H2_1'*H2_1)^(-1);     %Fj'
        F3_1 = H3_1*(H3_1'*H3_1)^(-1);     %Fj'
        F1_2 = H1_2*(H1_2'*H1_2)^(-1);     %Fj
        F2_2 = H2_2*(H2_2'*H2_2)^(-1);     %Fj'
        F3_2 = H3_2*(H3_2'*H3_2)^(-1);     %Fj'
        
        rho1 = P/trace((H1'*H1)^(-1));  %rhoj
        rho2 = P/trace((H2'*H2)^(-1));  %rhoj'
        rho3 = P/trace((H3'*H3)^(-1));  %rhoj'
        rho1_1 = P/trace((H1_1'*H1_1)^(-1));  %rhoj
        rho2_1 = P/trace((H2_1'*H2_1)^(-1));  %rhoj'
        rho3_1 = P/trace((H3_1'*H3_1)^(-1));  %rhoj'
        rho1_2 = P/trace((H1_2'*H1_2)^(-1));  %rhoj
        rho2_2 = P/trace((H2_2'*H2_2)^(-1));  %rhoj'
        rho3_2 = P/trace((H3_2'*H3_2)^(-1));  %rhoj'
        
        Intfe1 = zeros(K, 1);
        Intfe2 = zeros(K, 1);
        Intfe1_1 = zeros(K, 1);
        Intfe2_1 = zeros(K, 1);
        Intfe1_2 = zeros(K, 1);
        Intfe2_2 = zeros(K, 1);
        gamma = zeros(K,1);
        gamma1=zeros(K,1);
        gamma_1 = zeros(K,1);
        gamma_2 = zeros(K,1);
        
        numer=zeros(K,1);
        denom=zeros(K,1);
        sum=0;
        sum_1=0;
        sum_2=0;
        sum0=0;
        sum1=0;
        sum2=0;
        sum3=0;
        %         sum4=0;
        %         sum5=0;
        for k=1:K
            sum1=sum1+1/(M1*((alphatheta1_1(:,k))^2+(alphatheta1_5(:,k))^2+(alphatheta1_8(:,k))^2)+...
                M2*((alphatheta2_1(:,k))^2+(alphatheta2_5(:,k))^2+(alphatheta2_8(:,k))^2));
            sum2=sum2+1/(M1*((alphatheta1_4(:,k))^2+(alphatheta1_2(:,k))^2+(alphatheta1_9(:,k))^2)+...
                M2*((alphatheta2_4(:,k))^2+(alphatheta2_2(:,k))^2+(alphatheta2_9(:,k))^2));
            sum3=sum3+1/(M1*((alphatheta1_7(:,k))^2+(alphatheta1_3(:,k))^2+(alphatheta1_6(:,k))^2)+...
                M2*((alphatheta2_7(:,k))^2+(alphatheta2_3(:,k))^2+(alphatheta2_6(:,k))^2));
        end
        rho11=P/sum1;
        rho22=P/sum2;
        rho33=P/sum3;
        
        for k = 1 : K
            %式（41）：极限值
            numer(k,1)=((M1*((alphatheta1_1(:,k))^2)+M2*((alphatheta2_1(:,k))^2))/(M1*((alphatheta1_1(:,k))^2+(alphatheta1_5(:,k))^2+(alphatheta1_8(:,k))^2)+...   %式（36）
                M2*((alphatheta2_1(:,k))^2+(alphatheta2_5(:,k))^2+(alphatheta2_8(:,k))^2)))^2 * rho11;
            denom(k,1)=((M1*(alphatheta1_2(:,k))^2+M2*(alphatheta2_2(:,k))^2)/(M1*((alphatheta1_4(:,k))^2+(alphatheta1_2(:,k))^2+(alphatheta1_9(:,k))^2)+...
                M2*((alphatheta2_4(:,k))^2+(alphatheta2_2(:,k))^2+(alphatheta2_9(:,k))^2)))^2  * rho22 +...
                ((M1*((alphatheta1_3(:,k))^2)+M2*((alphatheta2_3(:,k))^2))/(M1*((alphatheta1_7(:,k))^2+(alphatheta1_3(:,k))^2+(alphatheta1_6(:,k))^2)+...
                M2*((alphatheta2_7(:,k))^2+(alphatheta2_3(:,k))^2+(alphatheta2_6(:,k))^2)))^2*rho33;
            gamma1(k,1)=numer(k,1)/(denom(k,1)+1);
            sum0=sum0+log2(1+gamma1(k,1));
            
            %式（35）：实际值
            Intfe1(k, 1) = (norm(H_21(:,k)' * F2))^2 * rho2...
                + (norm(H_31(:,k)' * F3))^2 * rho3;
            Intfe1_1(k, 1) = (norm(H_21_1(:,k)' * F2_1))^2 * rho2_1...
                + (norm(H_31_1(:,k)' * F3_1))^2 * rho3_1;
            Intfe1_2(k, 1) = (norm(H_21_2(:,k)' * F2_2))^2 * rho2_2...
                + (norm(H_31_2(:,k)' * F3_2))^2 * rho3_2;
            
            for k4=1:K
                if k4~=k
                    Intfe2(k4,1) = Intfe2(k4,1) + rho1*(abs(H_11(:,k)'*F1(:,k4)))^2;
                    Intfe2_1(k4,1) = Intfe2_1(k4,1) + rho1_1*(abs(H_11_1(:,k)'*F1_1(:,k4)))^2;
                    Intfe2_2(k4,1) = Intfe2_2(k4,1) + rho1_2*(abs(H_11_2(:,k)'*F1_2(:,k4)))^2;
                end
            end
            gamma(k, 1) = rho1*(abs(H_11(:,k)'*F1(:,k)))^2/(Intfe1(k,1)+Intfe2(k4,1)+1);
            gamma_1(k, 1) = rho1_1*(abs(H_11_1(:,k)'*F1_1(:,k)))^2/(Intfe1_1(k,1)+Intfe2_1(k4,1)+1);
            gamma_2(k, 1) = rho1_2*(abs(H_11_2(:,k)'*F1_2(:,k)))^2/(Intfe1_2(k,1)+Intfe2_2(k4,1)+1);
            sum=sum+log2(1+gamma(k, 1));
            sum_1=sum_1+log2(1+gamma_1(k, 1));
            sum_2=sum_2+log2(1+gamma_2(k, 1));
            
            %      %式42：下界值
            %             sum4=sum4+(M1*(alphatheta1_2(:,k))^2+M2*(alphatheta2_2(:,k))^2)^2;
            %             sum5=sum5+(1/(M1*(alphatheta1_4(:,k))^2+M2*(alphatheta2_4(:,k))^2))^2;
        end
        %         R(:,M/40)=R(:,M/40)+K*log2(1+(1/((L-1)*(sum4/K)*(sum5/K))));
        sumrate(:,K/5)=sumrate(:,K/5)+sum;
        sumrate1(:,K/5)=sumrate1(:,K/5)+sum0;
        sumrate_1(:,K/5)=sumrate_1(:,K/5)+sum_1;
        sumrate_2(:,K/5)=sumrate_2(:,K/5)+sum_2;
    end
end
sumrate(:,:)=sumrate(:,:)/Ni;
sumrate1(:,:)=sumrate1(:,:)/Ni;
sumrate_1(:,:)=sumrate_1(:,:)/Ni;
sumrate_2(:,:)=sumrate_2(:,:)/Ni;
% R(:,:)=R(:,:)/Ni;

plot(Qss,sumrate(:,:),'ks-'),xlabel('Number of the user K'),ylabel('Sum rate[bps/Hz]');
hold on;
plot(Qss,sumrate1(:,:),'gd-'),xlabel('Number of the user K'),ylabel('Sum rate[bps/Hz]');
plot(Qss,sumrate_1(:,:),'bx-'),xlabel('Number of the user K'),ylabel('Sum rate[bps/Hz]');
plot(Qss,sumrate_2(:,:),'cp-'),xlabel('Number of the user K'),ylabel('Sum rate[bps/Hz]');
% plot(Qss,R(:,:),'ro-'),xlabel('Number of the BS antenna M'),ylabel('Sum rate[bps/Hz]');
grid on;
legend('Actual value','Limit value','Center value','Two sector','Location','NorthWest');

% plot(A_1(:,1),A_1(:,2),'O');                                             %用户位置图
% hold on;
% plot(A_2(:,1)+r1,A_2(:,2)-r1*sqrt(3),'O');
% plot(A_3(:,1)-r1,A_3(:,2)-r1*sqrt(3),'O');