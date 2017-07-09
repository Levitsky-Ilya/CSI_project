N = 30;             % of subcarriers
N_hlf = 15;
M = 8;              % of antennas

c = 3e10;           % cm/s
f = 5320e6;         % Hz
f_delta = 1250e3;   % Hz
d = 1;              % cm

D = 2*pi*d*f/c;

n_thres = 4;
%threshold = ???

clf

x = [1:N];
hlf1 = [1 : N_hlf]; hlf2 = [(N_hlf+1) : N];

x_1 = [hlf1, hlf1, hlf1];
x_2 = [hlf2, hlf2, hlf2];

for j = 1:1

    ToF = 1.68e-8;   % s
    Noise = 0.3;
    
    csi_1 = ones(M,N)*10;  %-45o
    csi_2 = csi_1;          %-45o same,further
    for m = 1:M
        csi_1(m,:) = csi_1(m,:)*exp(-i*(m-1)*D*sin(pi/4)) + (rand(1,N)+i*rand(1,N))*Noise;
    end
    for n = 1:N
        csi_1(:,n) = csi_1(:,n)*exp(-i*(n-1)*(2*pi*f_delta)* ToF*2 ) + (rand(M,1)+i*rand(M,1))*Noise;
    end
    for m = 1:M
        csi_2(m,:) = csi_2(m,:)*exp(-i*(m-1)*D*sin(-pi/3)) + (rand(1,N)+i*rand(1,N))*Noise;
    end
    for n = 1:N
        csi_2(:,n) = csi_2(:,n)*exp(-i*(n-1)*(2*pi*f_delta)* ToF ) + (rand(M,1)+i*rand(M,1))*Noise;
    end
    
    csi = csi_1+csi_2;



    y1_1 = [angle(csi(1,hlf1))];                        %y(i)_j - i-th antenna
    y1_1 = unwrap(y1_1);                               % and j-th half of spectum
    y2_1 = [angle(csi(2,hlf1))];
    y2_1 = unwrap(y2_1);
    y3_1 = [angle(csi(3,hlf1))];
    y3_1 = unwrap(y3_1);
    
    y1_2 = [angle(csi(1,hlf2))];
    y1_2 = unwrap(y1_2);
    y2_2 = [angle(csi(2,hlf2))];
    y2_2 = unwrap(y2_2);
    y3_2 = [angle(csi(3,hlf2))];
    y3_2 = unwrap(y3_2);

    y_1 = [y1_1,y2_1,y3_1];
    y_2 = [y1_2,y2_2,y3_2];
    
    p_1 = polyfit(x_1, y_1, 1);                               %p(1) = 2 * pi * f_delta * tau_s
    p_2 = polyfit(x_1, y_1, 1); 
    %plot(x1,y3_1)
    
    for n = hlf1
        csi(:,n) = csi(:,n)*exp(-i*((n-1)*p_1(1)+p_1(2)));
    end
    for n = hlf2
        csi(:,n) = csi(:,n)*exp(-i*((n-1)*p_2(1)+p_2(2)));
    end

    %{
    y1_1 = [angle(csi???(1,hlf1))];
    y2_1 = [angle(csi(2,hl16f1))];
    y3_1 = [angle(csi(3,hlf1))];
    
    y1_2 = [angle(csi(1,hlf2))];
    y2_2 = [angle(csi(2,hlf2))];
    y3_2 = [angle(csi(3,hlf2))];
    
    y1 = unwrap([y1_1,y1_2]);
    y2 = unwrap([y2_1,y2_2]);
    y3 = unwrap(1[y3_1,y3_2]);

    plot(x,y1,'red', x,y2,'blue', x,y3,'black');
    hold on
    grid onexp(i*(1:14)*Tau
    %}
    
    
    for n = 1:15
        csi1(n,:) = csi(1,n:15+n);
        csi2(n,:) = csi(2,n:15+n);
        csi3(n,:) = csi(3,n:15+n);
    end

    X = [csi1,csi2;csi2,csi3];
    
    [E_N,eigv] = eig(X*X');
    
    E_N(:,n_thres+1:N) = [];
    
    Tau_res = 100; theta_res = 100;
    Tau_max = 2*pi/12; theta_max = pi/2;
    n_Tau = 0; n_theta = 0;                 % Tau = 2*pi*f_delta* tau
    P = zeros(theta_res,Tau_res);
    
    for theta = linspace(-theta_max,theta_max,theta_res)
        n_theta = n_theta+1;
        n_Tau = 0;
        for Tau = linspace(0,Tau_max,Tau_res)
            n_Tau = n_Tau+1;
            
            power = 0:1:14;
            a_H1 = exp(i*power*Tau);
            a_H = [a_H1, a_H1*exp(i*D*sin(theta))];
            a_H_E_N = a_H * E_N;
            P(n_theta,n_Tau) = 1/(a_H_E_N * a_H_E_N');
            
        end
    end
    
    %[theta,Tau] = meshgrid(linspace(0,2*pi,25),linspace(0,2*pi,20));
    %P(:,1:20)=[]; P(:,161:180)=[];
    
    s = surf(linspace(0,Tau_max,Tau_res),linspace(-theta_max,theta_max,theta_res),P);
    %s.EdgeColor = 'none';
    xlabel('Tau'); ylabel('theta');
    grid off
    
end
