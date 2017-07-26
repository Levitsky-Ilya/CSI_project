function P = P(theta_Tau,E_N,steering_length,steering_width,f,d)

c = 3e10;           % cm/s

D = 2*pi*d*f/c;
theta = theta_Tau(1);
Tau = theta_Tau(2);

power = 0:1:steering_length-1;
a_H1 = exp(1i*power*Tau);
a_H = [];
for m = 1:steering_width
	slm = steering_length*m;
	a_H(slm-steering_length+1:slm) = a_H1*exp((m-1)*1i*D*sin(theta));
end
a_H_E_N = a_H * E_N;
P = -1/(a_H_E_N * a_H_E_N');

end