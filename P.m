function P = P(theta_Tau,E_N)

c = 3e10;           % cm/s
f = 5320e6;         % Hz
d = 1;              % cm

D = 2*pi*d*f/c;
theta = theta_Tau(1);
Tau = theta_Tau(2);

power = 0:1:14;
a_H1 = exp(1i*power*Tau);
a_H = [a_H1, a_H1*exp(1i*D*sin(theta))];
a_H_E_N = a_H * E_N;
P = -1/(a_H_E_N * a_H_E_N');

end