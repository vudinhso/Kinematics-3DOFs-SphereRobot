# Kinematics-3DOFs-SphereRobot
AgileEye Kinematics - Forward and Inverse Kinematics. 
- InverseKinematics: calculate (th1, th2, th3) from a rotation matrix
- angleYXY: input: rotation matrix. output phi, theta, psi YXY Euler angle. 
- angleGA: input: rotation matrix. output: bivector (quaternion)
- RotMatrix: input: axis, angle. output: rotation matrix (3x3). 
- RotLog: input: rotation matrix. output axis angle. 
- RotZXY.m: input [phi1, phi2, phi3] (3x1). output rotation matrix (3x3).  
- RotGA.m: input bivector (4x1). output rotation matrix (3x3)
- FwdKinematics-GeoAlg.m: Forward kinematics with Geometric Algebra (no Gimbal lock)
- FwdKinematics-ZXY.m: Forward kinematics with ZXY (potential gimbal lock)
