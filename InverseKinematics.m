clc; clear all; close all, format compact, format shortG

# Manipulator Parameter
eta = [0, 2*pi/3, -2*pi/3];     % angular spacing between the unit vector ui and vi
alph1 = pi/2;                   % first link (proximal)
alph2 = pi/2;                   % second link (distal)
gama1 = pi/2;                   % angle between ui
gama2 = pi/2;                   % angle between vi
beta1 = asin( (2*sqrt(3)/3)*sin(gama1/2) );
beta2 = asin( (2*sqrt(3)/3)*sin(gama2/2) );

# Orientation of the end-effector. It requires the function:
# RotMatrix: give the 3x3 rotation matrix from an axis and an angle. 
# angleYXY: give the phi, theta, psi Eurler angle from the YXY convention. 
u  = [0,0,1]';
th = 15*pi/180;
R = RotMatrix(u,th)
q = angleYXY(R)

% Inertial frame: z is vertical, u1 belongs to xz and x is horizontal, y complete the frame. 
% Reference configuration for the unit vector of the effector: 
% v1, v2, v3 dont necessarily built an orthogonal frame. 
% v10: in the same plane as u1
% v20: in the same plane as u2
% v30: in the same plane as u3

# Unit axis
x = [1,0,0]';
y = [0,1,0]';
z = [0,0,1]';

# Calculation of the unit vector  u and vi
for i = 1:3
    u(:,i)  = RotMatrix(z,eta(i)) * RotMatrix(y,-beta1) * (-z); 
    v0(:,i) = RotMatrix(z,eta(i)) * RotMatrix(y, beta2) * ( z); 
end

# Calculation of the unit vector v after a rotation R
for i = 1:3
    v(:,i) = R * v0(:,i);
end
# Calculation of U, V, W
for i = 1:3
    U(i) = -sin(eta(i)) * sin(alph1) * v(1,i) + ...
            cos(eta(i)) * sin(alph1) * v(2,i); 
    V(i) =  cos(eta(i)) * sin(alph1) * cos(beta1) * v(1,i) + ...
            sin(eta(i)) * sin(alph1) * cos(beta1) * v(2,i) + ...
                          sin(alph1) * sin(beta1) * v(3,i); 
    W(i) =  cos(eta(i)) * cos(alph1) * sin(beta1) * v(1,i) + ...
            sin(eta(i)) * cos(alph1) * sin(beta1) * v(2,i) + ...
                          cos(alph1) * cos(beta1) * v(3,i) - cos(alph2); 
    A(i) = W(i) - U(i); 
    B(i) = V(i); 
    C(i) = W(i) + U(i); 
    th_plus(i) = 2*atan2(-B(i) + sqrt(B(i)^2-A(i)*C(i)), A(i)) ;
    th_minus(i)= 2*atan2(-B(i) - sqrt(B(i)^2-A(i)*C(i)), A(i)) ;
end
th_plus
th_minus

