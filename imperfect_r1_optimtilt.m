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
Qs = (1:10);
Qss=20*Qs;
P=10e+9;
R1=zeros(1,10);
R2=zeros(1,10);
for p=1:length(Qs)
    r1=20*p;
    sumrate = zeros(N,N);
    sumrate1 = zeros(N,N);
    for j = 1:N
        theta1 = ((j-1)/N)*(pi/2);
        for q = 1:N
            if q<j
                continue;
            end
            theta2 = ((q-1)/N)*(pi/2);
            H1 = zeros (M,K);
            H2 = zeros (M,K);
            H3 = zeros (M,K);
            F1 = zeros (M,K);
            F2 = zeros (M,K);
            F3 = zeros (M,K);
            rho11=0;
            rho22=0;
            rho33=0;
            for i = 1:Ni
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%以第1个小区为基准
                k1 = 0;
                A_1 = zeros (K,2);
                H_11 = zeros(M,K);                %第1个小区基站到其用户的信道Hjj
                H_21 = zeros(M,K);                %第2个小区基站到第一个小区用户的信道Hj'j
                H_31 = zeros(M,K);                %第3个小区基站到第一个小区用户的信道Hj'j
                alphatheta1_1=zeros(1,K);
                alphatheta2_1=zeros(1,K);
                alphatheta1_2=zeros(1,K);
                alphatheta2_2=zeros(1,K);
                alphatheta1_3=zeros(1,K);
                alphatheta2_3=zeros(1,K);
                alphas = zeros(K,8);
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
                        alphatheta1_11=10^(alphatheta1_11/10)*d_11^(-1.5);
                        alphatheta2_11=10^(alphatheta2_11/10)*d_11^(-1.5);
                        alphatheta1_1(:,k1+1)=alphatheta1_11;
                        alphatheta2_1(:,k1+1)=alphatheta2_11;
                        h1_11=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                        h2_11=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                        h_11=[alphatheta1_11*h1_11;alphatheta2_11*h2_11];
                        H_11(:,k1+1)=h_11;
                        
                        phi_21 = -((pi/3)-atan(abs(x1-r1)/(y1+(sqrt(3))*r1)));         %第2个小区基站到第一个小区用户的信道Hj'j
                        theta3_21 = (pi/2)-atan((sqrt((x1-r1)^2+(y1+(sqrt(3))*r1)^2))/h);
                        d_21=sqrt(h^2+(x1-r1)^2+(y1+(sqrt(3))*r1)^2);
                        alphatheta1_21 = -min(min(12*(phi_21/phi3dB)^2, SLLaz)...
                            + min(12*((theta3_21 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                        alphatheta2_21 = -min(min(12*(phi_21/phi3dB)^2, SLLaz)...
                            + min(12*((theta3_21 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                        alphatheta1_21=10^(alphatheta1_21/10)*d_21^(-1.5);
                        alphatheta2_21=10^(alphatheta2_21/10)*d_21^(-1.5);
                        alphatheta1_2(:,k1+1)=alphatheta1_21;
                        alphatheta2_2(:,k1+1)=alphatheta2_21;
                        h1_21=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                        h2_21=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                        h_21=[alphatheta1_21*h1_21;alphatheta2_21*h2_21];
                        H_21(:,k1+1)=h_21;
                        alphas(k1+1,1) = alphatheta1_21;
                        alphas(k1+1,2) = alphatheta2_21;
                        
                        phi_31 =atan((y1+(sqrt(3))*r1)/(x1+r1))-pi/6;                  %第3个小区基站到第一个小区用户的信道Hj'j
                        theta3_31 = (pi/2)-atan((sqrt((x1+r1)^2+(y1+(sqrt(3))*r1)^2))/h);
                        d_31=sqrt(h^2+(x1+r1)^2+(y1+(sqrt(3))*r1)^2);
                        alphatheta1_31 = -min(min(12*(phi_31/phi3dB)^2, SLLaz)...
                            + min(12*((theta3_31 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                        alphatheta2_31 = -min(min(12*(phi_31/phi3dB)^2, SLLaz)...
                            + min(12*((theta3_31 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                        alphatheta1_31=10^(alphatheta1_31/10)*d_31^(-1.5);
                        alphatheta2_31=10^(alphatheta2_31/10)*d_31^(-1.5);
                        alphatheta1_3(:,k1+1)=alphatheta1_31;
                        alphatheta2_3(:,k1+1)=alphatheta2_31;
                        h1_31=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                        h2_31=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                        h_31=[alphatheta1_31*h1_31;alphatheta2_31*h2_31];
                        H_31(:,k1+1)=h_31;
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
                        alphatheta1_22=10^(alphatheta1_22/10)*d_22^(-1.5);
                        alphatheta2_22=10^(alphatheta2_22/10)*d_22^(-1.5);
                        alphatheta1_4(:,k2+1)=alphatheta1_22;
                        alphatheta2_4(:,k2+1)=alphatheta2_22;
                        h1_22=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                        h2_22=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                        h_22=[alphatheta1_22*h1_22;alphatheta2_22*h2_22];
                        H_22(:,k2+1)=h_22;
                        alphas(k2+1,5) = alphatheta1_22;
                        alphas(k2+1,6) = alphatheta2_22;
                        
                        
                        phi_12 = atan((x2+r1)/(abs(y2-r1*sqrt(3))));                   %第1个小区基站到第2个小区用户的信道Hjj'
                        theta3_12 = (pi/2)-atan((sqrt((x2+r1)^2+(y2-r1*sqrt(3))^2))/h);
                        d_12=sqrt(h^2+(x2+r1)^2+(y2-r1*sqrt(3))^2);
                        alphatheta1_12 = -min(min(12*(phi_12/phi3dB)^2, SLLaz)...
                            + min(12*((theta3_12 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                        alphatheta2_12 = -min(min(12*(phi_12/phi3dB)^2, SLLaz)...
                            + min(12*((theta3_12 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                        alphatheta1_12=10^(alphatheta1_12/10)*d_12^(-1.5);
                        alphatheta2_12=10^(alphatheta2_12/10)*d_12^(-1.5);
                        alphatheta1_5(:,k2+1)=alphatheta1_12;
                        alphatheta2_5(:,k2+1)=alphatheta2_12;
                        h1_12=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                        h2_12=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                        h_12=[alphatheta1_12*h1_12;alphatheta2_12*h2_12];
                        H_12(:,k2+1)=h_12;
                        
                        phi_32 = atan(y2/(x2+2*r1))-pi/6;                              %第3个小区基站到第2个小区用户的信道Hj''j'
                        theta3_32 = (pi/2)-atan((sqrt((x2+2*r1)^2+y2^2))/h);
                        d_32=sqrt(h^2+(x2+2*r1)^2+y2^2);
                        alphatheta1_32= -min(min(12*(phi_32/phi3dB)^2, SLLaz)...
                            + min(12*((theta3_32 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                        alphatheta2_32 = -min(min(12*(phi_32/phi3dB)^2, SLLaz)...
                            + min(12*((theta3_32 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                        alphatheta1_32=10^(alphatheta1_32/10)*d_32^(-1.5);
                        alphatheta2_32=10^(alphatheta2_32/10)*d_32^(-1.5);
                        alphatheta1_6(:,k2+1)=alphatheta1_32;
                        alphatheta2_6(:,k2+1)=alphatheta2_32;
                        h1_32=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                        h2_32=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                        h_32=[alphatheta1_32*h1_32;alphatheta2_32*h2_32];
                        H_32(:,k2+1)=h_32;
                        
                        k2=k2+1;
                    end
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%以第3个小区为基准
                k3 = 0;
                A_3 = zeros(K, 2);
                H_33 = zeros(M,K);      %第3个小区基站到其自己用户的信道Hjj
                H_13 = zeros(M,K);      %第1个小区基站到第3个小区用户的信道Hj'j
                H_23 = zeros(M,K);      %第2个小区基站到第3个小区用户的信道Hj'j
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
                        alphatheta1_33=10^(alphatheta1_33/10)*d_33^(-1.5);
                        alphatheta2_33=10^(alphatheta2_33/10)*d_33^(-1.5);
                        alphatheta1_7(:,k3+1)=alphatheta1_33;
                        alphatheta2_7(:,k3+1)=alphatheta2_33;
                        h1_33=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                        h2_33=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                        h_33=[alphatheta1_33*h1_33;alphatheta2_33*h2_33];
                        H_33(:,k3+1)=h_33;
                        alphas(k3+1,7) = alphatheta1_33;
                        alphas(k3+1,8) = alphatheta2_33;
                        
                        phi_13 = -atan(abs(x3-r1)/(abs(y3-r1*sqrt(3))));               %第1个小区基站到第3个小区用户的信道Hjj'
                        theta3_13 = (pi/2)-atan((sqrt((x3-r1)^2+(y3-r1*sqrt(3))^2))/h);
                        d_13=sqrt(h^2+(x3-r1)^2+(y3-r1*sqrt(3))^2);
                        alphatheta1_13 = -min(min(12*(phi_13/phi3dB)^2, SLLaz)...
                            + min(12*((theta3_13 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                        alphatheta2_13 = -min(min(12*(phi_13/phi3dB)^2, SLLaz)...
                            + min(12*((theta3_13 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                        alphatheta1_13=10^(alphatheta1_13/10)*d_13^(-1.5);
                        alphatheta2_13=10^(alphatheta2_13/10)*d_13^(-1.5);
                        alphatheta1_8(:,k3+1)=alphatheta1_13;
                        alphatheta2_8(:,k3+1)=alphatheta2_13;
                        h1_13=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                        h2_13=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                        h_13=[alphatheta1_13*h1_13;alphatheta2_13*h2_13];
                        H_13(:,k3+1)=h_13;
                        
                        phi_23 = pi/6+atan((y3/(x3-2*r1)));                            %第2个小区基站到第3个小区用户的信道Hj''j'
                        theta3_23 = (pi/2)-atan((sqrt((x3-2*r1)^2+y3^2))/h);
                        d_23=sqrt(h^2+(x3-2*r1)^2+y3^2);
                        alphatheta1_23= -min(min(12*(phi_23/phi3dB)^2, SLLaz)...
                            + min(12*((theta3_23 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                        alphatheta2_23 = -min(min(12*(phi_23/phi3dB)^2, SLLaz)...
                            + min(12*((theta3_23 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                        alphatheta1_23=10^(alphatheta1_23/10)*d_23^(-1.5);
                        alphatheta2_23=10^(alphatheta2_23/10)*d_23^(-1.5);
                        alphatheta1_9(:,k3+1)=alphatheta1_23;
                        alphatheta2_9(:,k3+1)=alphatheta2_23;
                        h1_23=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                        h2_23=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                        h_23=[alphatheta1_23*h1_23;alphatheta2_23*h2_23];
                        H_23(:,k3+1)=h_23;
                        
                        k3=k3+1;
                    end
                end
                
                H1 = H_11 + H_12 + H_13;      %Hjj
                H2 = H_22 + H_21 + H_23;      %Hj'j'
                H3 = H_33 + H_31 + H_32;      %Hj'j'
                
                F1 = H1*(H1'*H1)^(-1);     %Fj
                F2 = H2*(H2'*H2)^(-1);     %Fj'
                F3 = H3*(H3'*H3)^(-1);     %Fj'
                
                rho1 = P/trace((H1'*H1)^(-1));  %rhoj
                rho2 = P/trace((H2'*H2)^(-1));  %rhoj'
                rho3 = P/trace((H3'*H3)^(-1));  %rhoj'
                
                Intfe1 = zeros(K, 1);
                Intfe2 = zeros(K, 1);
                gamma = zeros(K,1);
                gamma1=zeros(K,1);
                numer=zeros(K,1);
                denom=zeros(K,1);
                sum=0;
                sum0=0;
                sum1=0;
                sum2=0;
                sum3=0;
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
                    numer(k,1)=((M1*((alphatheta1_1(:,k))^2)+M2*((alphatheta2_1(:,k))^2))/(M1*((alphatheta1_1(:,k))^2+(alphatheta1_5(:,k))^2+(alphatheta1_8(:,k))^2)+...   %式（36）
                        M2*((alphatheta2_1(:,k))^2+(alphatheta2_5(:,k))^2+(alphatheta2_8(:,k))^2)))^2 * rho11;
                    denom(k,1)=((M1*(alphatheta1_2(:,k))^2+M2*(alphatheta2_2(:,k))^2)/(M1*((alphatheta1_4(:,k))^2+(alphatheta1_2(:,k))^2+(alphatheta1_9(:,k))^2)+...
                        M2*((alphatheta2_4(:,k))^2+(alphatheta2_2(:,k))^2+(alphatheta2_9(:,k))^2)))^2  * rho22 +...
                        ((M1*((alphatheta1_3(:,k))^2)+M2*((alphatheta2_3(:,k))^2))/(M1*((alphatheta1_7(:,k))^2+(alphatheta1_3(:,k))^2+(alphatheta1_6(:,k))^2)+...
                        M2*((alphatheta2_7(:,k))^2+(alphatheta2_3(:,k))^2+(alphatheta2_6(:,k))^2)))^2*rho33;
                    gamma1(k,1)=numer(k,1)/(denom(k,1)+1);
                    sum0=sum0+log2(1+gamma1(k,1));
                    
                    Intfe1(k, 1) = (norm(H_21(:,k)' * F2))^2 * rho2...             %式（35）
                        + (norm(H_31(:,k)' * F3))^2 * rho3;
                    for k4=1:K
                        if k4~=k
                            Intfe2(k4,1) = Intfe2(k4,1) + rho1*(abs(H_11(:,k)'*F1(:,k4)))^2;
                        end
                    end
                    gamma(k, 1) = rho1*(abs(H_11(:,k)'*F1(:,k)))^2/(Intfe1(k,1)+Intfe2(k4,1)+1);
                    sum=sum+log2(1+gamma(k, 1));
                end
                sumrate(j,q)=sumrate(j,q)+sum;
                sumrate1(j,q)=sumrate1(j,q)+sum0;
            end
        end
    end
    [sumrate1maxa, sumrate1maxb] = find(sumrate1==max(max(sumrate1)));
    R1(1,p)=((sumrate1maxa-1)/N)*(pi/2);
    R2(1,p)=((sumrate1maxb-1)/N)*(pi/2);
end
plot(Qss, R1,'rs-'),xlabel('Cell radius r1[m]'),ylabel('Tilt angle[rad]');
hold on;
plot(Qss, R2,'go-'),xlabel('Cell radius r1[m]'),ylabel('Tilt angle[rad]');
grid on;
legend('\theta_1','\theta_2','Location','NorthEast');

% theta1=(((1:N)-1)/N)*(pi/2);
% theta2=(((1:N)-1)/N)*(pi/2);
% mesh(theta1,theta2,sumrate1/Ni);