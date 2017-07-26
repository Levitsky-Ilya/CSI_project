%%%//////////////////////////////////////////////%%%
%%%       Programm for localisation of           %%%
%%%   transmitter using real CSI, obtained CSI   %%%
%%%                 by FPGA                      %%%
%%%            version for 802.11a               %%%
%%%                                              %%%
%%% Programmed by Levitsky Ilya 2017             %%%
%%%//////////////////////////////////////////////%%%

%%%=============== CONTROL PANEL ===============%%%
PACK_PROC   = 1;      % PACKET PROCCESSING enable: 1
PACKET_NO   = 5:6;      % Range for packets: 0 - all, other - PACKET_NO
CSI_EXTRACT = 1;      % Getting CSI enable: 1
MOD_PHASE   = 1;      % Modifying CSI phase enable: 1
PHASE_PLOT  = 1;      % Unwrapped and modified phase plots enable: 1
BUILD_X     = 1;      % Building X enable: 1
MAXIMA_P    = 0;      % Searching for Maxima enable: 1
P_PLOT      = 0;      % Function P plot enable: 1
CLUST_LIKE  = 0;      % CLUSTERING and LIKEL.EST. enable: 1

%%%=========== SETTING UP PARAMETERS ===========%%%

%%%%%%%%%%%%%%% Reciever Parameters %%%%%%%%%%%%%%%
N = 53;            % of subcarriers
M = 2;              % of antennas
n_thres = 4;
steering_len = 18;
steering_wid = 2;

Nsl = N-steering_len+1;
Msw = M-steering_wid+1;
X_len = Nsl*Msw;
X_wid  = steering_len*steering_wid;

%%%%%%%%%%%%%%% Channel Parameters %%%%%%%%%%%%%%%
%sections = 2;  %continious sections of subcarriers where info is sent
% boundaries of these sections
%cont_sect_bound = [1 26; 28 53];
strange_subc = 27;
section_len = 26;
x_low = 1;
x_up = x_low + section_len - 1;   

%%%%%%%%%%%%%%% Physical Constants %%%%%%%%%%%%%%%
c = 3e10;           % cm/s
f = 2.437e9;        % Hz
f_delta = 312500;   % Hz
d = 11.5;           % cm

D = 2*pi*d*f/c;
D_1 = 2*pi*f_delta/c;   %Tau = D_1*distance

%%%%%%%%%%%%%%%% CSI Processing %%%%%%%%%%%%%%
epsilon = 1e-4;
subcarriers = 1:N;
%section = zeros(sections,section_len);

xM = [];
for m = 1:M
    xM = [xM,subcarriers];
end

%%%%%%%%%%%%%%%% Problem boundariues %%%%%%%%%%%%%%
dist_res = 100; theta_res = 200;
dist_min = 0; dist_max = 1000;    % Maximum distance is 24000 cm
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
if PACK_PROC

fileID = fopen('CSI_output_a.txt'); packets_max = 100; time_intv_max = 100000;

j = 0; packet = 0; 
csi_all = zeros(packets_max,M,N);
timestamp = zeros(1,2);

while ~feof(fileID)
    j = j+1;
    C_info(j,:) = textscan(fileID,'Timestamp:%d, Antenna %d');
    C_nsubc = textscan(fileID,repmat('%d',[1,N]),1,'CollectOutput',1);
    C_subc{j} = textscan(fileID,repmat('%f',[1,N]),1,'CollectOutput',1);
    
    
    antenna = C_info{j,2};
    timestamp(antenna) = C_info{j,1};
    
    if isempty(timestamp(antenna))
        continue
    end
    csi_all(packet+1,antenna,:) = C_subc{1,j}{1,1};
    
    if (csi_all(packet+1,3-antenna,1) ~= 0 && ...
            abs(timestamp(3-antenna)-timestamp(antenna)) < time_intv_max)
        packet = packet+1;
        if packet == packets_max
            break;
        end
    end
end %while
fclose(fileID);

packets = packet;
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
        yM1 = angle(csi(m,subcarriers));
        yM1(strange_subc) = (yM1(strange_subc+1)+yM1(strange_subc-1))/2;
        yM1 = unwrap(yM1);
        yM = [yM,yM1];
    end
    
    for m = 1:M %every single antenna!
    
        xlm = x_low + (m-1)*N;
        xum = xlm + section_len - 1;
        p = polyfit(subcarriers(x_low:x_up),yM(xlm:xum),1);

        turn = 2; % moving up mode.
        y2 = yM(xum+2); y1 = yM(xum);
        x2 = subcarriers(x_up+2);  x1 = subcarriers(x_up);
        K = (y2-y1)/(x2-x1); y_shift = 0;

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
            xlms = xlm+section_len+1;
            xums = xlms+section_len-1;
            yM(xlms:xums) = yM(xlms:xums)+y_shift;
        end

        yM(xum+1) = (yM(xum+2) + yM(xum))/2;

    end %for m
    
    p = polyfit(xM, yM, 1);                 % p(1) = 2 * pi * f_delta * tau_s   
    for n = subcarriers
        csi(:,n) = csi(:,n)*exp(-1i*((n-1)*p(1)+p(2)));
    end
    for m = 1:M
        csi(m,strange_subc) = (csi(m,strange_subc-1)+csi(m,strange_subc+1))/2;
    end
      
    
if PHASE_PLOT  
    colour = ['r','b'];
    figure(1)
    for m = 1:M
        plot(subcarriers,yM((1:N)+(m-1)*N),colour(m));
        hold on;
    end 
    %plot(subcarriers, subcarriers*p(1)+p(2));
    grid on;
    xlabel('Subcarrier index'); ylabel('phase');

    figure(2)   % not working proparly!
    for m = 1:M
        y = unwrap(angle(csi(m,subcarriers))); 
        plot(subcarriers,y,colour(m));
        hold on;
    end
    grid on;
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
w_dist = 0.00002;
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
    -w_dist*dist_var; -w_short*abs(C(:,1)'); Likelihood];
'dist'
OUT(1,:)
'theta'
OUT(2,:)
'weigthed cluster size' 
OUT(3,:)
'weigthed theta_varian' 
OUT(4,:)
'weigthed dist_varian ' 
OUT(5,:)
'weigthed distance    '
OUT(6,:)
'Total likelihood     '
OUT(7,:)

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
