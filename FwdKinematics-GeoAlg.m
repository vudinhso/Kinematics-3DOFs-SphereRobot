# Determine the Rotation Matrix from th1, th2, th3. 
# Use Geometric Algebra to avoid Euler angle sinularity (gimbal lock)
# Equivalent to quaternions

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
phi1 = 0; phi2 = 30; phi3 = 15;
R_ref = RotZXY([phi1, phi2, phi3]); 

# Initial guess for the orientation. It should not be too far from the actual orientation of the effector. 
q = angleGA(Rref)+.01*rand(4,1);

# Initialization of the error. 
err = 1;

# Newton-Gauss loop
j=0; iter = 15;
while (j<iter) && (err > 1e-4)
    j=j+1; 
    R = RotGA(q); 
    a = q(1); b = q(2); c = q(3); d = q(4);
    # Constraint equation f
    for i = 1:3
        f(i,1) = w(:,i)' * R * v0(:,i) - cos(alph2);
    end
    # bivector constraint equation
    f(4,1) = a^2 + b^2 + c^2 + d^2 -1;
    err = norm(f);
    # Jacobian matrix J
    dRda = [    0, -2*d, 2*c; 
                2*d, 0, -2*b; 
                -2*c, 2*b, 0]; 
    dRdb = [    0, 2*c, 2*d; 
                2*c, -4*b, -2*a; 
                2*d, 2*a, -4*b]; 
    dRdc = [    -4*c, 2*b, 2*a; 
                2*b, 0, 2*d; 
                -2*a, 2*d, -4*c];
    dRdd = [    -4*d, -2*a, 2*b; 
                2*a, -4*d, 2*c; 
                2*b, 2*c, 0];
    J11 = w(:,1)' * dRda * v0(:,1);
    J21 = w(:,2)' * dRda * v0(:,2);
    J31 = w(:,3)' * dRda * v0(:,3);
    J41 = 2*a;

    J12 = w(:,1)' * dRdb * v0(:,1);
    J22 = w(:,2)' * dRdb * v0(:,2);
    J32 = w(:,3)' * dRdb * v0(:,3);
    J42 = 2*b;

    J13 = w(:,1)' * dRdc * v0(:,1);
    J23 = w(:,2)' * dRdc * v0(:,2);
    J33 = w(:,3)' * dRdc * v0(:,3);
    J43 = 2*c;

    J14 = w(:,1)' * dRdd * v0(:,1);
    J24 = w(:,2)' * dRdd * v0(:,2);
    J34 = w(:,3)' * dRdd * v0(:,3);
    J44 = 2*d;

    J = [   J11, J12, J13, J14; 
            J21, J22, J23, J24; 
            J31, J32, J33, J34; 
            J41, J42, J43, J44]; 
    q = q - J\f;
end

# q is containing 4 parameters (bivector-quaternion)
# It is not so intuitive, thus it can be converted into rotation matrix 
R = RotGA(q); 
