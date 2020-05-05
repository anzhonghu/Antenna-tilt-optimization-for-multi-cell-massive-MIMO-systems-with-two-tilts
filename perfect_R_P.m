%容量随功率的变化
clear all;
close all;
h = 30;
r1 = 120;
r2 = 10;
M1 = 50;
M2 = 50;
M = 100;
K = 10;
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
sumrate2 = zeros(1,5);
sumrate3 = zeros(1,5);
th1 = pi/2 - atan(r1/h);
th2 = pi/2 - atan(r2/h);
dthe = th2 - th1;
Qs = (8:12);
Ni = 1e2;
for i=1:Ni
    k1 = 0;
    A_0 = zeros (K, 2);
    H_0 = zeros (M,K);
    H_2 = zeros (M,K);
    H_3 = zeros (M,K);
    H_01 = zeros (M,K);
    H_21 = zeros (M,K);
    H_31 = zeros (M,K);
    H_02 = zeros (M,K);
    H_22 = zeros (M,K);
    H_32 = zeros (M,K);
    sum1 = 0;
    while k1<K
        x0 = [-(sqrt(3)/2)*r1 (sqrt(3)/2)*r1];
        y0 = [0 -r1];
        x0 = (max(x0)-min(x0))*rand(1,1)+min(x0);
        y0 = (max(y0)-min(y0))*rand(1,1)+min(y0);
        
        if (r2<=sqrt(x0^2+y0^2)) && (sqrt(x0^2+y0^2)<=r1) && atan((abs(x0))/(abs(y0)))<=pi/3
            A_0(k1+1,1)=x0;
            A_0(k1+1,2)=y0;
            phi_0 = atan(x0/(abs(y0)));                                    %Hjj
            theta3_0 = (pi/2)-atan(sqrt(x0^2+y0^2)/h);
            d_0=sqrt(h^2+x0^2+y0^2);
            alphatheta1_0= -min(min(12*(phi_0/phi3dB)^2, SLLaz)...
                + min(12*((theta3_0 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta2_0 = -min(min(12*(phi_0/phi3dB)^2, SLLaz)...
                + min(12*((theta3_0 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta1_01= -min(min(12*(phi_0/phi3dB)^2, SLLaz)...
                + min(12*((theta3_0 - theta3)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta2_01 = -min(min(12*(phi_0/phi3dB)^2, SLLaz)...
                + min(12*((theta3_0 - theta4)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta1_02= -min(min(12*(phi_0/phi3dB)^2, SLLaz)...
                + min(12*((theta3_0 - theta5)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta2_02 = -min(min(12*(phi_0/phi3dB)^2, SLLaz)...
                + min(12*((theta3_0 - theta6)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            
            alphatheta1_0=10^(alphatheta1_0/10) * d_0^(-1.5);
            alphatheta2_0=10^(alphatheta2_0/10) * d_0^(-1.5);
            alphatheta1_01=10^(alphatheta1_01/10) * d_0^(-1.5);
            alphatheta2_01=10^(alphatheta2_01/10) * d_0^(-1.5);
            alphatheta1_02=10^(alphatheta1_02/10) * d_0^(-1.5);
            alphatheta2_02=10^(alphatheta2_02/10) * d_0^(-1.5);
            
            h1_0=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
            h2_0=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
            h_0=[alphatheta1_0*h1_0;alphatheta2_0*h2_0];
            h_01=[alphatheta1_01*h1_0;alphatheta2_01*h2_0];
            h_02=[alphatheta1_02*h1_0;alphatheta2_02*h2_0];
            
            H_0(:,k1+1)=h_0;
            H_01(:,k1+1)=h_01;
            H_02(:,k1+1)=h_02;
            
            phi_2 = -((pi/3)-atan(abs(x0-r1)/(y0+(sqrt(3))*r1)));              % Hj'j
            theta3_2 = (pi/2)-atan(sqrt((x0-r1)^2+(y0+(sqrt(3))*r1)^2)/h);
            d_2=sqrt(h^2+(x0-r1)^2+(y0+(sqrt(3))*r1)^2);
            alphatheta1_2 = -min(min(12*(phi_2/phi3dB)^2, SLLaz)...
                + min(12*((theta3_2 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta2_2 = -min(min(12*(phi_2/phi3dB)^2, SLLaz)...
                + min(12*((theta3_2 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta1_21 = -min(min(12*(phi_2/phi3dB)^2, SLLaz)...
                + min(12*((theta3_2 - theta3)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta2_21 = -min(min(12*(phi_2/phi3dB)^2, SLLaz)...
                + min(12*((theta3_2 - theta4)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta1_22 = -min(min(12*(phi_2/phi3dB)^2, SLLaz)...
                + min(12*((theta3_2 - theta5)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta2_22 = -min(min(12*(phi_2/phi3dB)^2, SLLaz)...
                + min(12*((theta3_2 - theta6)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            
            alphatheta1_2=10^(alphatheta1_2/10) * d_2^(-1.5);
            alphatheta2_2=10^(alphatheta2_2/10) * d_2^(-1.5);
            alphatheta1_21=10^(alphatheta1_21/10) * d_2^(-1.5);
            alphatheta2_21=10^(alphatheta2_21/10) * d_2^(-1.5);
            alphatheta1_22=10^(alphatheta1_22/10) * d_2^(-1.5);
            alphatheta2_22=10^(alphatheta2_22/10) * d_2^(-1.5);
            
            h1_2=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
            h2_2=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
            
            h_2=[alphatheta1_2.*h1_2;alphatheta2_2.*h2_2];
            h_21=[alphatheta1_21.*h1_2;alphatheta2_21.*h2_2];
            h_22=[alphatheta1_22.*h1_2;alphatheta2_22.*h2_2];
            
            H_2(:,k1+1)=h_2;
            H_21(:,k1+1)=h_21;
            H_22(:,k1+1)=h_22;
            
            phi_3 =-((pi/3)-atan((x0+r1)/(abs(y0+(sqrt(3))*r1))));               %Hj'j
            theta3_3 = (pi/2)-atan(sqrt((x0+r1)^2+(y0+(sqrt(3))*r1)^2)/h);
            d_3=sqrt(h^2+(x0+r1)^2+(y0+(sqrt(3))*r1)^2);
            alphatheta1_3 = -min(min(12*(phi_3/phi3dB)^2, SLLaz)...
                + min(12*((theta3_3 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta2_3 = -min(min(12*(phi_3/phi3dB)^2, SLLaz)...
                + min(12*((theta3_3 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta1_31 = -min(min(12*(phi_3/phi3dB)^2, SLLaz)...
                + min(12*((theta3_3 - theta3)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta2_31 = -min(min(12*(phi_3/phi3dB)^2, SLLaz)...
                + min(12*((theta3_3 - theta4)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta1_32 = -min(min(12*(phi_3/phi3dB)^2, SLLaz)...
                + min(12*((theta3_3 - theta5)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta2_32 = -min(min(12*(phi_3/phi3dB)^2, SLLaz)...
                + min(12*((theta3_3 - theta6)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            
            alphatheta1_3=10^(alphatheta1_3/10) * d_3^(-1.5);
            alphatheta2_3=10^(alphatheta2_3/10) * d_3^(-1.5);
            alphatheta1_31=10^(alphatheta1_31/10) * d_3^(-1.5);
            alphatheta2_31=10^(alphatheta2_31/10) * d_3^(-1.5);
            alphatheta1_32=10^(alphatheta1_32/10) * d_3^(-1.5);
            alphatheta2_32=10^(alphatheta2_32/10) * d_3^(-1.5);
            
            h1_3=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
            h2_3=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
            
            h_3=[alphatheta1_3.*h1_3;alphatheta2_3.*h2_3];
            h_31=[alphatheta1_31.*h1_3;alphatheta2_31.*h2_3];
            h_32=[alphatheta1_32.*h1_3;alphatheta2_32.*h2_3];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            sum1 = sum1 + 1 /(M1*(alphatheta1_0^2) + M2*(alphatheta2_0^2));
            H_3(:,k1+1)=h_3;
            H_31(:,k1+1)=h_31;
            H_32(:,k1+1)=h_32;
            
            k1 = k1+1;
        end
    end
    
    k2 = 0;                             %第一个Hj'j'
    A_1 = zeros (K, 2);
    H_1=zeros(M,K);
    H_11=zeros(M,K);
    H_12=zeros(M,K);
    while k2<K
        x1 = [-(sqrt(3)/2)*r1 (sqrt(3)/2)*r1];
        y1 = [0 -r1];
        x1 = (max(x1)-min(x1))*rand(1,1)+min(x1);
        y1 = (max(y1)-min(y1))*rand(1,1)+min(y1);
        
        if (r2<=sqrt(x1^2+y1^2)) && (sqrt(x1^2+y1^2)<=r1) && atan((abs(x1))/(abs(y1)))<=pi/3
            A_1(k2+1,1)=x1;
            A_1(k2+1,2)=y1;
            phi_1 = atan(x1/(abs(y1)));
            theta3_1 = (pi/2)-atan(sqrt(x1^2+y1^2)/h);
            d_1=sqrt(h^2+x1^2+y1^2);
            alphatheta1_1 = -min(min(12*(phi_1/phi3dB)^2, SLLaz)...
                + min(12*((theta3_1 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta2_1 = -min(min(12*(phi_1/phi3dB)^2, SLLaz)...
                + min(12*((theta3_1 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta1_11 = -min(min(12*(phi_1/phi3dB)^2, SLLaz)...
                + min(12*((theta3_1 - theta3)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta2_11 = -min(min(12*(phi_1/phi3dB)^2, SLLaz)...
                + min(12*((theta3_1 - theta4)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta1_12 = -min(min(12*(phi_1/phi3dB)^2, SLLaz)...
                + min(12*((theta3_1 - theta5)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta2_12 = -min(min(12*(phi_1/phi3dB)^2, SLLaz)...
                + min(12*((theta3_1 - theta6)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            
            alphatheta1_1=10^(alphatheta1_1/10)*d_1^(-1.5);
            alphatheta2_1=10^(alphatheta2_1/10)*d_1^(-1.5);
            alphatheta1_11=10^(alphatheta1_11/10)*d_1^(-1.5);
            alphatheta2_11=10^(alphatheta2_11/10)*d_1^(-1.5);
            alphatheta1_12=10^(alphatheta1_12/10)*d_1^(-1.5);
            alphatheta2_12=10^(alphatheta2_12/10)*d_1^(-1.5);
            
            h1_1=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
            h2_1=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
            
            h_1=[alphatheta1_1.*h1_1;alphatheta2_1.*h2_1];
            h_11=[alphatheta1_11.*h1_1;alphatheta2_11.*h2_1];
            h_12=[alphatheta1_12.*h1_1;alphatheta2_12.*h2_1];
            
            H_1(:,k2+1)=h_1;
            H_11(:,k2+1)=h_11;
            H_12(:,k2+1)=h_12;
            
            k2=k2+1;
        end
    end
    
    k4 = 0;                                %第二个Hj'j'
    A_4 = zeros (K, 2);
    H_4 = zeros(M,K);
    H_41 = zeros(M,K);
    H_42 = zeros(M,K);
    while k4<K
        x4 = [-(sqrt(3)/2)*r1 (sqrt(3)/2)*r1];
        y4 = [0 -r1];
        x4 = (max(x4)-min(x4))*rand(1,1)+min(x4);
        y4 = (max(y4)-min(y4))*rand(1,1)+min(y4);
        
        if (r2<=sqrt(x4^2+y4^2)) && (sqrt(x4^2+y4^2)<=r1) && atan((abs(x4))/(abs(y4)))<=pi/3
            A_4(k4+1,1)=x4;
            A_4(k4+1,2)=y4;
            phi_4 = atan(x4/(abs(y4)));
            theta3_4 = (pi/2)-atan(sqrt(x4^2+y4^2)/h);
            d_4=sqrt(h^2+x4^2+y4^2);
            alphatheta1_4 = -min(min(12*(phi_4/phi3dB)^2, SLLaz)...
                + min(12*((theta3_4 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta2_4 = -min(min(12*(phi_4/phi3dB)^2, SLLaz)...
                + min(12*((theta3_4 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta1_41 = -min(min(12*(phi_4/phi3dB)^2, SLLaz)...
                + min(12*((theta3_4 - theta3)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta2_41 = -min(min(12*(phi_4/phi3dB)^2, SLLaz)...
                + min(12*((theta3_4 - theta4)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta1_42 = -min(min(12*(phi_4/phi3dB)^2, SLLaz)...
                + min(12*((theta3_4 - theta5)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            alphatheta2_42 = -min(min(12*(phi_4/phi3dB)^2, SLLaz)...
                + min(12*((theta3_4 - theta6)/theta3dB)^2, SLLel), SLLtot) + alphamax;
            
            alphatheta1_4=10^(alphatheta1_4/10)*d_4^(-1.5);
            alphatheta2_4=10^(alphatheta2_4/10)*d_4^(-1.5);
            alphatheta1_41=10^(alphatheta1_41/10)*d_4^(-1.5);
            alphatheta2_41=10^(alphatheta2_41/10)*d_4^(-1.5);
            alphatheta1_42=10^(alphatheta1_42/10)*d_4^(-1.5);
            alphatheta2_42=10^(alphatheta2_42/10)*d_4^(-1.5);
            
            h1_4=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
            h2_4=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
            
            h_4=[alphatheta1_4.*h1_4;alphatheta2_4.*h2_4];
            h_41=[alphatheta1_41.*h1_4;alphatheta2_41.*h2_4];
            h_42=[alphatheta1_42.*h1_4;alphatheta2_42.*h2_4];
            
            H_4(:,k4+1)=h_4;
            H_41(:,k4+1)=h_41;
            H_42(:,k4+1)=h_42;
            
            k4=k4+1;
        end
    end
    
    for q=1:length(Qs)
        P=10^Qs(q);
        Q0=H_0'*H_0;
        rho1=P/trace(Q0^(-1));
        Q1=H_1'*H_1;
        rho2=P/trace(Q1^(-1));
        Q4=H_4'*H_4;
        rho3=P/trace(Q4^(-1));
        
        Q01=H_01'*H_01;
        rho11=P/trace(Q01^(-1));
        Q11=H_11'*H_11;
        rho21=P/trace(Q11^(-1));
        Q41=H_41'*H_41;
        rho31=P/trace(Q41^(-1));
        
        Q02=H_02'*H_02;
        rho12=P/trace(Q02^(-1));
        Q12=H_12'*H_12;
        rho22=P/trace(Q12^(-1));
        Q42=H_42'*H_42;
        rho32=P/trace(Q42^(-1));
        
        Intfe = zeros(K, 1);
        Intfe1 = zeros(K, 1);
        Intfe2 = zeros(K, 1);
        
        gamma = zeros(K, 1);
        gamma1 = zeros(K, 1);
        gamma2 = zeros(K, 1);
        
        sum=0;
        sum2=0;
        sum3=0;
        for k = 1 : K
            Intfe(k, 1) = (norm(H_2(:,k)' * H_1 / (H_1' * H_1)))^2 * rho2...
                +norm(H_3(:,k)' * H_4 / (H_4' * H_4))^2 * rho3;
            Intfe1(k, 1) = (norm(H_21(:,k)' * H_11 / (H_11' * H_11)))^2 * rho21...
                +norm(H_31(:,k)' * H_41 / (H_41' * H_41))^2 * rho31;
            Intfe2(k, 1) = (norm(H_22(:,k)' * H_12 / (H_12' * H_12)))^2 * rho22...
                +norm(H_32(:,k)' * H_42 / (H_42' * H_42))^2 * rho32;
            
            gamma(k, 1) = rho1/(Intfe(k, 1)+1);
            gamma1(k, 1) = rho11/(Intfe1(k, 1)+1);
            gamma2(k, 1) = rho12/(Intfe2(k, 1)+1);
            
            %gamma(k, 1) = rho1;
            sum=sum+log2(1+gamma(k, 1));
            sum2=sum2+log2(1+gamma1(k, 1));
            sum3=sum3+log2(1+gamma2(k, 1));
        end
        sumrate(:,q)=sumrate(:,q)+sum;
        sumrate1(:,q)=sumrate1(:,q)+K*log2(1+P/sum1);
        sumrate2(:,q)=sumrate2(:,q)+sum2;                  %倾角指向中心
        sumrate3(:,q)=sumrate3(:,q)+sum3;
    end
end
sumrate(:,:)=sumrate(:,:)/Ni;
sumrate1(:,:)=sumrate1(:,:)/Ni;
sumrate2(:,:)=sumrate2(:,:)/Ni;
sumrate3(:,:)=sumrate3(:,:)/Ni;

plot(Qs,sumrate(:,:),'-ks'),xlabel('SNR[dB]'),ylabel('Sum rate[bps/Hz]');
hold on;
plot(Qs,sumrate1(:,:),'-gd'),xlabel('SNR[dB]'),ylabel('Sum rate[bps/Hz]');
plot(Qs,sumrate2(:,:),'bx-'),xlabel('SNR[dB]'),ylabel('Sum rate[bps/Hz]');
plot(Qs,sumrate3(:,:),'cp-'),xlabel('SNR[dB]'),ylabel('Sum rate[bps/Hz]');


sum1 = 0;
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
g = sum1;

R1=zeros(1,5);
R=0;
for q=1:length(Qs)
    P=10^Qs(q);
    R=K*log2(1+(P/g));
    R1(:,q)=R;
end
plot(Qs,R1(:,:),'-ro'),xlabel('SNR[dB]'),ylabel('Sum rate[bps/Hz]');
grid on;
legend('Actual value','Limit value','Center value','Two sector','Lower bound','Location','NorthWest');