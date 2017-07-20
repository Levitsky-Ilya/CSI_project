%%//////////////////////////////////////////////%%
%%       Theoretical testing programm           %%
%%             with modeled CSI                 %%
%%                                              %%
%% Properties: reprogrammable amount of         %%
%%  subcarriers and antennas. For now           %%
%%  parameters of smoothed CSI and steering     %%
%%  vector are edited manually.                 %%
%%                                              %%
%% Programmed by Levitsky Ilya 2017             %%
%%//////////////////////////////////////////////%%

%%%=========== SETTING UP PARAMETERS ===========%%%

%%%%%%%%%%%%%%% Reciever Parameters %%%%%%%%%%%%%%%
N = 116;            % of subcarriers
N_hlf = floor(N/2);
M = 2;              % of antennas
n_thres = 4;
packets = 10;
steering_length = 39;
steering_width = 2;

Nsl = N-steering_length+1;
Msw = M-steering_width+1;
X_length = Nsl*Msw;
X_width  = steering_length*steering_width;

%%%%%%%%%%%%%%% Physical Constants %%%%%%%%%%%%%%%
c = 3e10;           % cm/s
f = 5320e6;         % Hz
f_delta = 1250e3;   % Hz
d = 1;              % cm

D = 2*pi*d*f/c;
D_1 = 2*pi*f_delta/c;

Noise = 0.03;

%%%%%%%%%%%%%%%% CSI Processing %%%%%%%%%%%%%%
x = 1:N;
hlf1 = 1 : N_hlf; hlf2 = (N_hlf+1) : N;

x_1 = []; x_2 = [];
for m = 1:M
  x_1 = [x_1,hlf1];
  x_2 = [x_2,hlf2];
end

%%%%%%%%%%%%%%%% Problem boundariues %%%%%%%%%%%%%%
dist_res = 100; theta_res = 200;
dist_min = 0; dist_max = 2000; theta_max = pi/2;    % Maximum distance is 24000 cm
theta_max = pi/2;

dist_cover = dist_max-dist_min;
theta_cover = 2*theta_max;

Tau_max = D_1 * dist_max;
Tau_min = D_1 * dist_min;

eps_theta = theta_cover/(theta_res+1);
eps_dist = dist_cover/(dist_res+1);

%%%%%%%%%%%%%%%% Maxima Searching Param. %%%%%%%%%%%%%%
dist_units = 4; theta_units = 4;
pts_in_unit = 5;
points = dist_units * theta_units * pts_in_unit;
    
sector_length = (Tau_max-Tau_min)/dist_units;
sector_width  = 2*theta_max/theta_units;

dist_coords = [];
theta_coords = [];

for j = 1:packets

    %ToF = 1.68e-8;   % s
    
    %{
    % Test 1 %
    sin1 = -60/335.4;   dist1 = 335.4;  ampl1 = 5;
    sin2 = 480/582.5;   dist2 = 582.5;  ampl2 = 5/5;
    sin3 = -360/488.4;  dist3 = 488.4;  ampl3 = 5/5;
    sin4 = -60/454.0;   dist4 = 454.0;  ampl4 = 5/5;
    %}
    
    % Test 2 %
    sin1 = 60/335.1;    dist1 = 335.1;  ampl1 = 5;
    sin2 = 220/413.4;   dist2 = 413.4;  ampl2 = 5/5;
    sin3 = -340/488.0;  dist3 = 488.0;  ampl3 = 5/5;
    sin4 = 60/553.3;    dist4 = 553.3;  ampl4 = 5/5;
    
    %{
    % Test 3 %
    sin1 = 0.01;        dist1 = 200;    ampl1 = 5;
    sin2 = 270/336;     dist2 = 336;    ampl2 = 5/5;
    sin3 = -290/352.3;	dist3 = 352.3;  ampl3 = 5/5;
    sin4 = 0.005;       dist4 = 700;    ampl4 = 5/5;
    %}
    %{
    % Test 4 %
    sin1 = 90/174;      dist1 = 174;    ampl1 = 5;
    sin2 = 190/242.1;   dist2 = 242.1;  ampl2 = 5/5;
    sin3 = -370/339.2;	dist3 = 339.2;  ampl3 = 5/5;
    sin4 = 130/761.2;   dist4 = 761.2;  ampl4 = 5/5;
    %}
    
    csi_1 = ones(M,N)*ampl1; 
    csi_2 = ones(M,N)*ampl2;
    csi_3 = ones(M,N)*ampl3;
    csi_4 = ones(M,N)*ampl4;
    
    for m = 1:M
        csi_1(m,:) = csi_1(m,:)*exp(-1i*(m-1)*D* sin1) + (rand(1,N)+1i*rand(1,N))*Noise;
    end
    for n = 1:N
        csi_1(:,n) = csi_1(:,n)*exp(-1i*(n-1)*D_1* dist1) + (rand(M,1)+1i*rand(M,1))*Noise;
    end
    
    for m = 1:M
        csi_2(m,:) = csi_2(m,:)*exp(-1i*(m-1)*D* sin2) + (rand(1,N)+1i*rand(1,N))*Noise;
    end
    for n = 1:N
        csi_2(:,n) = csi_2(:,n)*exp(-1i*(n-1)*D_1* dist2-1i*pi) + (rand(M,1)+1i*rand(M,1))*Noise;
    end
    
    for m = 1:M
        csi_3(m,:) = csi_3(m,:)*exp(-1i*(m-1)*D* sin3) + (rand(1,N)+1i*rand(1,N))*Noise;
    end
    for n = 1:N
        csi_3(:,n) = csi_3(:,n)*exp(-1i*(n-1)*D_1* dist3-1i*pi) + (rand(M,1)+1i*rand(M,1))*Noise;
    end
    
    for m = 1:M    
        csi_4(m,:) = csi_4(m,:)*exp(-1i*(m-1)*D* sin4-1i*pi) + (rand(1,N)+1i*rand(1,N))*Noise;
    end
    for n = 1:N
        csi_4(:,n) = csi_4(:,n)*exp(-1i*(n-1)*D_1* dist4) + (rand(M,1)+1i*rand(M,1))*Noise;
    end
    
    
    csi = csi_1+csi_2+csi_3+csi_4;


         
    %{
    
    y_1 = zeros(1,N_hlf*M);
    y_2 = zeros(1,N_hlf*M);
    
    for m = 1:M
      top = N_hlf*m;
      bottom = top + 1 - N_hlf;
      
      y_1(bottom:top) = unwrap(angle(csi(m,hlf1)));
      y_2(bottom:top) = unwrap(angle(csi(m,hlf2)));
      
    end
    
    p_1 = polyfit(x_1, y_1, 1);                               %p(1) = 2 * pi * f_delta * tau_s
    p_2 = polyfit(x_1, y_2, 1); 
    %}
    
    %{
    figure(1)
    for m = 1:M
      top = N_hlf*m;
      bottom = top + 1 - N_hlf;
      
      plot(x,[y_1(bottom:top),y_2(bottom:top)]);
      hold on;
    end
    grid on; hold off;
    xlabel('Subcarrier index'); ylabel('phase');
    %}
    %{
    for n = hlf1
        csi(:,n) = csi(:,n)*exp(-1i*((n-1)*p_1(1)+p_1(2)));
    end
    for n = hlf2    
        csi(:,n) = csi(:,n)*exp(-1i*((n-1)*p_2(1)+p_2(2)));
    end
    %}
    %{
    figure(2)       ;???!!!
    for m = 1:M
      Y(m,:) = unwrap([angle(csi(m,hlf1)),angle(csi(m,hlf2))]);
      plot(x,Y(m,:));
      hold on;
    end
    grid on; hold off;
    %}
    
    
    X_sections = zeros(M,steering_length,Nsl);
    
    for m = 1:M
      for n = 1:steering_length
        X_sections(m,n,:) = csi(m,n:Nsl+n-1);
      end
    end

    X = [];
    for n = 1:steering_width
      
      X_1 = [];
      for m = n:Msw+n-1
        X_1 = [X_1,X_sections(m,:,:)];
      end
      X_1 = reshape(X_1,steering_length,Nsl*Msw);
      X = [X;X_1];
      
    end
    
    [E_N,eigv] = eig(X*X');
    E_N(:,n_thres+1:X_width) = [];
    
    %{
    P1 = zeros(theta_res,dist_res);
    
    n_theta = 0;
    for theta = linspace(-theta_max,theta_max,theta_res)
        n_theta = n_theta+1;
        
        n_dist = 0;
        for Tau = linspace(Tau_min,Tau_max,dist_res)
            n_dist = n_dist+1;
            
            power = 0:1:(steering_length-1);
            a_H1 = exp(1i*power*Tau);
            a_H = [];
            for m = 1:steering_width
              slm = steering_length*m;
              a_H(slm-steering_length+1:slm) = a_H1*exp((m-1)*1i*D*sin(theta));
            end
            a_H_E_N = a_H * E_N;
            P1(n_theta,n_dist) = 1/(a_H_E_N * a_H_E_N');
            
        end
    end
  
    figure(3);
    s = surf(linspace(dist_min,dist_max,dist_res),linspace(-theta_max,theta_max,theta_res),P1);
    xlabel('Distance'); ylabel('theta');
    hold on;
    %}
    
    point = 1; theta_chunk = zeros(1); dist_chunk = zeros(1);
    Tau_low = Tau_min;
    
    for Tau_unit = 1:dist_units
        
        theta_low = - theta_max;
        for theta_unit = 1:theta_units
            
            for N_point = 1:pts_in_unit
          
                Tau0   = rand*sector_length + Tau_low;
                theta0 = rand*sector_width + theta_low;
                P_par = @(theta_Tau)P(theta_Tau,E_N,steering_length,steering_width);
                theta_Tau = [theta0,Tau0];
                [theta_Tau0,Pval,exitflag,output] = fminsearch(P_par,theta_Tau);

                if (exitflag > 0)
                  
                  while (theta_Tau0(1) > theta_max+eps_theta || ...
                      theta_Tau0(1) < -theta_max-eps_theta)
                    
                    theta_shift = sign(theta_Tau0(1))*pi;
                    theta_Tau0(1) = -(theta_Tau0(1)-theta_shift);
                  end %while
                   
                  theta_Tau0(2) = theta_Tau0(2) / D_1;
                
                  %%!!!HARD.Make better!!!
                  if min(abs(theta_chunk-theta_Tau0(1))) > eps_theta ||...
                     min(abs(dist_chunk-theta_Tau0(2))) > eps_dist
               
                    point = point+1;
                    theta_chunk(point) = theta_Tau0(1);
                    dist_chunk(point) = theta_Tau0(2);
              
                  end
                  %%!!!HARD.Make better!!!
                  
                else 1
                end
                
            end %for N_point
            theta_low = theta_low + sector_width;
        end %for theta_unit
        Tau_low = Tau_low + sector_length;
    end %for Tau_unit
    
    theta_chunk(1)=[]; dist_chunk(1) = [];
    theta_coords = [theta_coords,theta_chunk];
    dist_coords = [dist_coords,dist_chunk];
    
end %for j
hold off

%%%=========== CLUSTERING BEGGINING ===========%%%

replics = 5; clusters = 4;

Y = [dist_coords;theta_coords];
Y(1,:) = Y(1,:)/dist_cover;         % prescaling to use..
Y(2,:) = Y(2,:)/theta_cover;        % ..sqeuclidean distance
Y = Y';
%opts = statset('Display','final');
[idx,C] = kmeans(Y,clusters,'Replicates',replics);%,'Options',opts);

% postscaling
C(:,1) = C(:,1)*dist_cover; C(:,2) = C(:,2)*theta_cover;


figure;
colour = ['r.','b.','g.','y.'];
plot(dist_coords(idx==1),theta_coords(idx==1),'r.','MarkerSize',12)
hold on
plot(dist_coords(idx==2),theta_coords(idx==2),'b.','MarkerSize',12)
plot(dist_coords(idx==3),theta_coords(idx==3),'g.','MarkerSize',12)
plot(dist_coords(idx==4),theta_coords(idx==4),'y.','MarkerSize',12)
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3)
title 'Cluster Assignments and Centroids'
hold off


cluster_size = zeros(1,clusters);
dist_var = zeros(1,clusters);
theta_var = zeros(1,clusters);

%w_s = ?
%w_theta = ?
%w_dist = ?
%w_short = ?
%{
for j = 1:clusters
    
    sz = size(dist_coords(idx==j));
    cluster_size(j) = sz(2);
    dist_var(j) = var(dist_coords(idx==j));
    theta_var(j) = var(theta_coords(idx==j));
    Likelihood(j) = exp(w_s*cluster_size(j) - w_theta*theta_var(j) - ...
                    w_dist*dist_var(j) - w_short*C(j,1));
    
end
%}