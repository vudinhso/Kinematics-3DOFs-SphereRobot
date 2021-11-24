# Determine the Rotation Matrix from th1, th2, th3. 
# Use Matrix multiplication and Euler ZXY

# unit axis 
x = [1,0,0]';
y = [0,1,0]';
z = [0,0,1]';

# Parameters of the mechanism
eta = [0, 2*pi/3, -2*pi/3];     % angular spacing between the unit vector ui and vi
alph1 = pi/2;                   % first link (proximal)
alph2 = pi/2;                   % second link (distal)
gama1 = pi/2;                   % angle between ui
gama2 = pi/2;                   % angle between vi
beta1 = asin( (2*sqrt(3)/3)*sin(gama1/2) );
beta2 = asin( (2*sqrt(3)/3)*sin(gama2/2) );

# Calculation of the unit vector u, v0, and w
for i = 1:3
    u(:,i)  = RotMatrix(z,eta(i)) * RotMatrix(y,-beta1) * (-z); 
    v0(:,i) = RotMatrix(z,eta(i)) * RotMatrix(y, beta2) * ( z); 
    w(:,i)  = RotMatrix(z,eta(i)) * ...
              RotMatrix(y,pi/2-beta1) * ...
              RotMatrix(x,th(i)) * ...
              RotMatrix(z,alph1) * x;
end

# Orientation reference. phi1, phi2, phi3 is used for the Euler angle ZXY representation
phi1 = 0; phi2 = 30; phi3 = 15; q = [phi1, phi2, phi3]; 
R_ref = RotZXY(q); 

# Initial guess for the orientation. It should not be too far from the actual orientation of the effector. 
q = qref + .1*rand(3,1); 

# Initialization of the error. 
err = 1;

# Newton-Gauss loop
j=0; iter = 15;
while (j<iter) && (err > 1e-4)
    j=j+1
    R = RotZXY(q);
    # Constraint equation f
    for i = 1:3
        f(i,1) = w(:,i)' * R * v0(:,i) - cos(alph2);
    end
    err = norm(f);
    # Jacobian matrix J
    Q    =  RotMatrix(y, q(3))' * VecToso3(x) * RotMatrix(y, q(3));
    dRd1 =  VecToso3(z) * R;
    dRd2 =  R * Q;
    dRd3 =  R * VecToso3(y);
    J = [   w(:,1)' * dRd1 * v0(:,1), ...  
            w(:,1)' * dRd2 * v0(:,1), ...
            w(:,1)' * dRd3 * v0(:,1); 
            w(:,2)' * dRd1 * v0(:,2), ...
            w(:,2)' * dRd2 * v0(:,2), ...
            w(:,2)' * dRd3 * v0(:,2); 
            w(:,3)' * dRd1 * v0(:,3), ...
            w(:,3)' * dRd2 * v0(:,3), ...
            w(:,3)' * dRd3 * v0(:,3)];
    q = q - J\f;
end

# q is containing 3 parameters (Euler angle ZXY)
# It can be converted into rotation matrix 
R = RotZXY(q); 
