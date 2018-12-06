clear; close all; clc;
addpath(genpath('~/rvctools'));
addpath('functions');

%--- Dataset parameters
gravity = 9.81007;          % Gravity magnitude, m/s^2

%--- Container of the results
N = 10;
XX = zeros(N, 4);

%--- Initialization
x = [ 0.999988079071  -0.000324862776  -0.000137553128  -0.004879518412 ]';

dx = zeros(3,1);
P = [[ 0.00025591533631086349487304687500  0.00000000000872402237256375556739  0.00000003580366936262180388439446];[0.00000000000873242970989007716298  0.00025591524899937212467193603515  0.00000012538744442736060591414570];[0.00000003580369778433123428840190  0.00000012538750127077946672216057  0.00077918911119922995567321777343 ]];





qq = zeros(4, N);
bb = zeros(3, N);

noise_gyro = 0.0000036509;
noise_accel = 0.0105740093;

w     = [ 0.000815593114  -0.000174769943  0.000599211256 ]';
acc_m = [ 0.005143792834  0.007435118780  9.815170288086 ]';
deltat = 0.071883000433;

% EKF - Transition matrices
A_x = 1/2*Omega(w);
F_x = eye(4) + A_x*deltat;
x = F_x*x;
x = x/norm(x);

% ESKF - Transition matrices
F_dx = eye(3) + ...
    skew(w)*deltat;
F_i = eye(3);
Q_i = noise_gyro*deltat^2*eye(3);

% Error state and covariance update
dx = F_dx*dx;
P = F_dx*P*F_dx' + F_i*Q_i*F_i';

% 2. Correction
acc = acc_m;
acc = acc / norm(acc);

% Measurement Jacobians
H_x = 2*[-x(3)    x(4)    -x(1)   x(2); ...
    x(2)    x(1)     x(4)   x(3); ...
    x(1)   -x(2)    -x(3)   x(4)];
X_dx = Qmat(x);
H_dx = H_x*X_dx;

Na = (noise_accel/norm(acc))^2*eye(3);
Na = noise_accel*eye(3);
Z = H_dx*P*H_dx' + Na;
K = P*H_dx'/Z;

R = fromqtoR(x);
acc_est = -R'*[0;0;-gravity];
acc_est = acc_est/norm(acc_est);
z = acc - acc_est;
dx = K*z;
S = eye(3) - K*H_dx;
P = S*P*S' + K*Z*K';

% Reset operation
x(1:4) = leftQuaternion(x(1:4))*[1; dx(1:3)/2];

G = eye(3) - skew(0.5*dx(1:3));
P = G*P*G';

dx = zeros(3,1);