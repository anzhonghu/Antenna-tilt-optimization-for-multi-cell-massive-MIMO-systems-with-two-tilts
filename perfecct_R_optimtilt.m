%角度与容量关系
clear all;
close all;
h = 30;
r1 = 120;
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
sumrate = zeros(N,N);
sumrate1 = zeros(N,N);
th1 = pi/2 - atan(r1/h);
th2 = pi/2 - atan(r2/h);
dthe = th2 - th1;
Ni = 1e2;
R1=zeros(N,N);
P=10e+9;

for j = 1:N
    theta1 = ((j-1)/N)*(pi/2);
    for q = 1:N
        theta2 = ((q-1)/N)*(pi/2);
        for i=1:Ni
            k1 = 0;
            A_0 = zeros (K, 2);
            H_0=zeros(M,K);
            H_2=zeros(M,K);
            H_3=zeros(M,K);
            sum1 = 0;
            while k1<K
                x0 = [-(sqrt(3)/2)*r1 (sqrt(3)/2)*r1]; 
                y0 = [0 -r1];  
                x0 = (max(x0)-min(x0))*rand(1,1)+min(x0);
                y0 = (max(y0)-min(y0))*rand(1,1)+min(y0);

                if (r2<=sqrt(x0^2+y0^2)) && (sqrt(x0^2+y0^2)<=r1) && atan((abs(x0))/(abs(y0)))<=pi/3
                    A_0(k1+1,1)=x0;
                    A_0(k1+1,2)=y0;
                    phi_0 = atan(x0/(abs(y0)));                     %Hjj
                    theta3_0 = (pi/2)-atan(sqrt(x0^2+y0^2)/h);
                    d_0=sqrt(h^2+x0^2+y0^2);
                    alphatheta1_0= -min(min(12*(phi_0/phi3dB)^2, SLLaz)...
                        + min(12*((theta3_0 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                    alphatheta2_0 = -min(min(12*(phi_0/phi3dB)^2, SLLaz)...
                        + min(12*((theta3_0 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                    alphatheta1_0=10^(alphatheta1_0/10) * d_0^(-1.5);
                    alphatheta2_0=10^(alphatheta2_0/10) * d_0^(-1.5);
                    h1_0=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                    h2_0=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                    h_0=[alphatheta1_0*h1_0;alphatheta2_0*h2_0];
                    H_0(:,k1+1)=h_0;

                    phi_2 = -((pi/3)-atan(abs(x0-r1)/(y0+(sqrt(3))*r1)));              % Hj'j
                    theta3_2 = (pi/2)-atan(sqrt((x0-r1)^2+(y0+(sqrt(3))*r1)^2)/h);
                    d_2=sqrt(h^2+(x0-r1)^2+(y0+(sqrt(3))*r1)^2);
                    alphatheta1_2 = -min(min(12*(phi_2/phi3dB)^2, SLLaz)...
                        + min(12*((theta3_2 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                    alphatheta2_2 = -min(min(12*(phi_2/phi3dB)^2, SLLaz)...
                        + min(12*((theta3_2 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                    alphatheta1_2=10^(alphatheta1_2/10) * d_2^(-1.5);
                    alphatheta2_2=10^(alphatheta2_2/10) * d_2^(-1.5);
                    h1_2=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                    h2_2=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                    h_2=[alphatheta1_2.*h1_2;alphatheta2_2.*h2_2];
                    H_2(:,k1+1)=h_2;

                    phi_3 =-((pi/3)-atan((x0+r1)/(abs(y0+(sqrt(3))*r1))));               %Hj'j
                    theta3_3 = (pi/2)-atan(sqrt((x0+r1)^2+(y0+(sqrt(3))*r1)^2)/h);
                    d_3=sqrt(h^2+(x0+r1)^2+(y0+(sqrt(3))*r1)^2);
                    alphatheta1_3 = -min(min(12*(phi_3/phi3dB)^2, SLLaz)...
                        + min(12*((theta3_3 - theta1)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                    alphatheta2_3 = -min(min(12*(phi_3/phi3dB)^2, SLLaz)...
                        + min(12*((theta3_3 - theta2)/theta3dB)^2, SLLel), SLLtot) + alphamax;
                    alphatheta1_3=10^(alphatheta1_3/10) * d_3^(-1.5);
                    alphatheta2_3=10^(alphatheta2_3/10) * d_3^(-1.5);
                    h1_3=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                    h2_3=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                    h_3=[alphatheta1_3.*h1_3;alphatheta2_3.*h2_3];
                    %%%%%%%%%%%%%%%%%%%%%%%%%%
                    sum1 = sum1 + 1 /(M1*(alphatheta1_0^2) + M2*(alphatheta2_0^2)); 
                    H_3(:,k1+1)=h_3;
                    k1 = k1+1;
                end
            end

            k2 = 0;                             %第一个Hj'j'
            A_1 = zeros (K, 2);
            H_1=zeros(M,K);
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

                    alphatheta1_1=10^(alphatheta1_1/10)*d_1^(-1.5);
                    alphatheta2_1=10^(alphatheta2_1/10)*d_1^(-1.5);

                    h1_1=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                    h2_1=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                    h_1=[alphatheta1_1.*h1_1;alphatheta2_1.*h2_1];
                    H_1(:,k2+1)=h_1;
                    k2=k2+1;
                end
            end

            k4 = 0;                                %第二个Hj'j'
            A_4 = zeros (K, 2);
            H_4 = zeros(M,K);
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

                    alphatheta1_4=10^(alphatheta1_4/10)*d_4^(-1.5);
                    alphatheta2_4=10^(alphatheta2_4/10)*d_4^(-1.5);

                    h1_4=sqrt(1/2)*(randn(M1,1)+sqrt(-1)*randn(M1,1));
                    h2_4=sqrt(1/2)*(randn(M2,1)+sqrt(-1)*randn(M2,1));
                    h_4=[alphatheta1_4.*h1_4;alphatheta2_4.*h2_4];
                    H_4(:,k4+1)=h_4;
                    k4=k4+1;
                end
            end

             
            Q0=H_0'*H_0;
            rho1=P/trace(Q0^(-1));
            Q1=H_1'*H_1;
            rho2=P/trace(Q1^(-1));
            Q4=H_4'*H_4;
            rho3=P/trace(Q4^(-1));

            Intfe = zeros(K, 1);
            gamma = zeros(K, 1);
            sum=0;
            for k = 1 : K
                Intfe(k, 1) = (norm(H_2(:,k)' * H_1 / (H_1' * H_1)))^2 * rho2...
                    +norm(H_3(:,k)' * H_4 / (H_4' * H_4))^2 * rho3;
                gamma(k, 1) = rho1/(Intfe(k, 1)+1);
                %gamma(k, 1) = rho1;
                sum=sum+log2(1+gamma(k, 1));
            end
            sumrate(j,q)=sumrate(j,q)+sum;
            sumrate1(j,q)=sumrate1(j,q)+K*log2(1+P/sum1);
        end
    end
end
sumrate=sumrate/Ni;
sumrate1=sumrate1/Ni;

%figure(1)
%theta1=(((1:N)-1)/N)*(pi/2);
%theta2=(((1:N)-1)/N)*(pi/2);
%mesh(theta1,theta2,sumrate);
%plot(Qs,sumrate(j,q)/Ni,'k-'),xlabel('Pj'),ylabel('sumrate');
%hold on;
%plot(Qs,sumrate1(j,q)/Ni,'g-'),xlabel('Pj'),ylabel('sumrate');

    

for j = 1:N
    theta1 = ((j-1)/N)*(pi/2);
    for q = 1:N
        theta2 = ((q-1)/N)*(pi/2);
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
       % g(j,k) = sum1;
       % disp([j,k])
        R=K*log2(1+(P/sum1));
        R1(j,q)=R;
    end
end

theta1=(((1:N)-1)/N)*(pi/2);
theta2=(((1:N)-1)/N)*(pi/2);
figure(1)
mesh(theta1,theta2,R1);
figure(2)
mesh(theta1,theta2,sumrate);
figure(3)
mesh(theta1,theta2,sumrate1);






%plot(Qs,R1(:,:),'r'),xlabel('Pj'),ylabel('sumrate');