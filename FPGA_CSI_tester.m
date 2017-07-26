%%%//////////////////////////////////////////////%%%
%%%       Programm for localisation of           %%%
%%%       transmitter using modeled CSI          %%%
%%%                                              %%%
%%% Programmed by Levitsky Ilya 2017             %%%
%%%//////////////////////////////////////////////%%%

%%%=============== CONTROL PANEL ===============%%%
PACK_PROC   = 1;      % PACKET PROCCESSING enable: 1
PACKET_NO   = 0;      % Range for packets: 0 - all, other - PACKET_NO:PACKET_NO
CSI_EXTRACT = 1;      % Getting CSI enable: 1
MOD_PHASE   = 1;      % Modifying CSI phase enable: 1
PHASE_PLOT  = 0;      % Unwrapped and modified phase plots enable: 1
BUILD_X     = 1;      % Building X enable: 1
MAXIMA_P    = 0;      % Searching for Maxima enable: 1
P_PLOT      = 1;      % Function P plot enable: 1
CLUST_LIKE  = 0;      % CLUSTERING and LIKEL.EST. enable: 1

%%%=========== SETTING UP PARAMETERS ===========%%%

%%%%%%%%%%%%%%% Reciever Parameters %%%%%%%%%%%%%%%
N = 117;            % of subcarriers
M = 2;              % of antennas
n_thres = 4;
steering_len = 39;
steering_wid = 2;

Nsl = N-steering_len+1;
Msw = M-steering_wid+1;
X_len = Nsl*Msw;
X_wid  = steering_len*steering_wid;

%%%%%%%%%%%%%%% Channel Parameters %%%%%%%%%%%%%%%
sections = 4;  %continious sections of subcarriers where info is sent
% boundaries of these sections
cont_sect_bound = [1 26; 28 53; 65 90; 92 117];
section_len = 26;

%%%%%%%%%%%%%%% Physical Constants %%%%%%%%%%%%%%%
c = 3e10;           % cm/s
f = 2.437e9;        % Hz
f_delta = 312500;   % Hz
d = 1;           % cm

D = 2*pi*d*f/c;
D_1 = 2*pi*f_delta/c;   %Tau = D_1*distance

%%%%%%%%%%%%%%%% CSI Processing %%%%%%%%%%%%%%
epsilon = 1e-4;
subcarriers = 1:N;
section = zeros(sections,section_len);

x = [];
for j = 1:sections
    x = [x,cont_sect_bound(j,1):cont_sect_bound(j,2)];
end

xM = [];
for m = 1:M
    xM = [xM,x];
end

%%%%%%%%%%%%%%%% Problem boundariues %%%%%%%%%%%%%%
dist_res = 100; theta_res = 400;
dist_min = -500; dist_max = 2000;    % Maximum distance is 24000 cm
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
    
sector_len = (Tau_max-Tau_min)/dist_units;
sector_wid  = 2*theta_max/theta_units;

dist_coords = [];
theta_coords = [];

%%%=========== PACKET PROCCESSING BEGGINING ===========%%%
packets = 20;

    
if PACK_PROC
    Noise = 0.00001;
    
for packet = 1:packets
    
    %{
    % Test 1 %
    sin1 = -60/335.4;   dist1 = 335.4;  ampl1 = 5;
    sin2 = 480/582.5;   dist2 = 582.5;  ampl2 = 5/5;
    sin3 = -360/488.4;  dist3 = 488.4;  ampl3 = 5/5;
    sin4 = -60/454.0;   dist4 = 454.0;  ampl4 = 5/5;
    %}
    
    % Test 2 %
    sin1 = 60/335.1;    dist1 = 335.1;  ampl1 = 0.003;
    sin2 = 220/413.4;   dist2 = 413.4;  ampl2 = 0.003/5;
    sin3 = -340/488.0;  dist3 = 488.0;  ampl3 = 0.003/5;
    sin4 = 60/553.3;    dist4 = 553.3;  ampl4 = 0.003/5;
    
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
    
    
    csi_all(packet,:,:) = csi_2;
end %for j
    
end %PACK_PROC


if PACKET_NO
  packet_range = PACKET_NO;
else
  packet_range = 1:packets;
end %PACKETS_NO
for j = packet_range

%%%%%%%%%%%%%%%% Getting CSI %%%%%%%%%%%%%%    
if CSI_EXTRACT        
    csi = csi_all(j,:,:);
    csi = reshape(csi,M,N);
end %CSI_EXTRACT 
%%%%%%%%%%%%%%%% Modifying CSI phase %%%%%%%%%%%%%%
if MOD_PHASE 
    yM = [];
    for m = 1:M
        yM1 = angle(csi(m,x));
        yM1 = unwrap(yM1);
        yM = [yM,yM1];
    end
    
    
    y_shift = 0;
    x_low = 1;
    for k = 1:sections-1
        x_up = x_low + section_len - 1;
        
        p = polyfit(x(x_low:x_up),yM(x_low:x_up),1);
        
        turn = 2; % moving up mode.
        y2 = yM(x_up+1); y1 = yM(x_up);
        x2 = x(x_up+1);  x1 = x(x_up);
        K = (y2-y1)/(x2-x1);
        
        while turn ~= 0
            if turn == 2    % moving up mode.
                K_new = ((y2 + 2*pi + y_shift) - y1)/(x2-x1);
                if abs(K_new-p(1)) < abs(K-p(1))
                    K = K_new; y_shift = y_shift+2*pi;
                else
                    turn = 1;   % moving down mode.
                end
            elseif turn == 1    % moving down mode.
                K_new = ((y2 - 2*pi + y_shift) - y1)/(x2-x1);
                if abs(K_new-p(1)) < abs(K-p(1))
                    K = K_new; y_shift = y_shift-2*pi;
                else
                    turn = 0;
                end
            elseif turn ~= 0
                turn = 0;
            end %if turn
        end %while turn
        
        if abs(y_shift) > epsilon
            for m = 1:M
                xlm = x_low+section_len+(m-1)*length(x);
                xum = xlm + section_len - 1;
                yM(xlm:xum) = yM(xlm:xum)+y_shift;
            end
        end
        
        x_low = x_up+1;
 
    end %for k

    
    p = polyfit(xM, yM, 1);                 % p(1) = 2 * pi * f_delta * tau_s   
    for n = x
        csi(:,n) = csi(:,n)*exp(-1i*((n-1)*p(1)+p(2)));
    end
      
if PHASE_PLOT   
    figure(1)
    for m = 1:M
        plot(x,yM((1:length(x))+(m-1)*length(x)));
        hold on;
    end 
    plot(x, x*p(1)+p(2));
    grid on; hold off;
    xlabel('Subcarrier index'); ylabel('phase');
    
    figure(2)   % not working proparly!
    for m = 1:M
        y = unwrap(angle(csi(m,x))); 
        plot(x,y);
        hold on;
    end
    grid on; hold off;
    xlabel('Subcarrier index'); ylabel('phase');

end %PHASE_PLOT 
    
end %MOD_PHASE     
%%%%%%%%%%%%%%%% Building X - smoothed CSI %%%%%%%%%%%%%%
if BUILD_X
    X_sections = zeros(M,steering_len,Nsl);
    
    for m = 1:M
      for n = 1:steering_len
        X_sections(m,n,:) = csi(m,n:Nsl+n-1);
      end
    end

    X = [];
    for n = 1:steering_wid
      
      X_1 = [];
      for m = n:Msw+n-1
        X_1 = [X_1,X_sections(m,:,:)];
      end
      X_1 = reshape(X_1,steering_len,Nsl*Msw);
      X = [X;X_1];
      
    end
    [E_N,eigv] = eig(X*X');
    E_N = E_N(:,1:n_thres);
end %BUILD_X    
%%%%%%%%%%%%%%%% Searching for Maxima %%%%%%%%%%%%%%
if P_PLOT
    P_data = zeros(theta_res,dist_res);
    n_theta = 0; 
    theta_grid = linspace(-theta_max,theta_max,theta_res);
    Tau_grid = linspace(Tau_min,Tau_max,dist_res);
    dist_grid = linspace(dist_min,dist_max,dist_res);
    for theta = theta_grid
        n_theta = n_theta + 1;
        
        n_Tau = 0;
        for Tau = Tau_grid
            n_Tau = n_Tau + 1;
            
            P_data(n_theta,n_Tau) = -P([theta,Tau],E_N,steering_len,steering_wid,f,d);
        end
    end
    
    surf(dist_grid,theta_grid,P_data);
end %P_PLOT
   
if MAXIMA_P  
    point = 1; theta_chunk = zeros(1); dist_chunk = zeros(1);
    Tau_low = Tau_min;
    for Tau_unit = 1:dist_units
        
        theta_low = - theta_max;
        for theta_unit = 1:theta_units
            
            for N_point = 1:pts_in_unit
          
                Tau0   = rand*sector_len + Tau_low;
                theta0 = rand*sector_wid + theta_low;
                P_par = @(theta_Tau)P(theta_Tau,E_N,steering_len,steering_wid,f,d);
                [theta_Tau0,Pval] = fminsearch(P_par,[theta0,Tau0]);

                while theta_Tau0(1) > theta_max+eps_theta || ...
                      theta_Tau0(1) < -theta_max-eps_theta
                    
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
                
            end %for N_point
            theta_low = theta_low + sector_wid;
        end %for theta_unit
        Tau_low = Tau_low + sector_len;
    end %for Tau_unit
    
    theta_chunk(1)=[]; dist_chunk(1) = [];
    theta_coords = [theta_coords,theta_chunk];
    dist_coords = [dist_coords,dist_chunk];
    
end %MAXIMA_P
end %for j
hold off;


%%%=========== CLUSTERING BEGGINING ===========%%%
if CLUST_LIKE
replics = 5; clusters = 8;

Y = [dist_coords;theta_coords];
Y(1,:) = Y(1,:)/dist_cover;         % prescaling to use..
Y(2,:) = Y(2,:)/theta_cover;        % ..sqeuclidean distance
Y = Y';
%opts = statset('Display','final');
[idx,C] = kmeans(Y,clusters,'Replicates',replics);%,'Options',opts);

% postscaling
C(:,1) = C(:,1)*dist_cover; C(:,2) = C(:,2)*theta_cover;

%%%=========== LIKELIHOOD ESTIMATION ===========%%%
cluster_size = zeros(1,clusters);
dist_var = zeros(1,clusters);
theta_var = zeros(1,clusters);
Likelihood = zeros(1,clusters);

w_s = 0.003;
w_theta = 0.002;
w_dist = 0.002;
w_short = 0.001;

for j = 1:clusters
    
    sz = size(dist_coords(idx==j));
    cluster_size(j) = sz(2);
    dist_var(j) = var(dist_coords(idx==j));
    theta_var(j) = var(theta_coords(idx==j));
    Likelihood(j) = exp(w_s*cluster_size(j) - w_theta*theta_var(j) - ...
                    w_dist*dist_var(j) - w_short*abs(C(j,1)));
end
[Most_like,max_like_idx] = max(Likelihood);
OUT = [C'; w_s*cluster_size; -w_theta*theta_var; ...
    -w_dist*dist_var; -w_short*(C(:,1)'); Likelihood]

figure(4);
colour = 'kx';
for j = 1:clusters
    plot(dist_coords(idx==j),theta_coords(idx==j),'.','MarkerSize',12);
    
    if j == max_like_idx
        colour = 'rx';
    else
        colour = 'kx';
    end
    plot(C(j,1),C(j,2),colour,'MarkerSize',15,'Linewid',3);
    hold on;
end

title 'Cluster Assignments and Centroids';
hold off;
end %CLUST_LIKE
