%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%Motion estimation using 4 points algorithm%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%                 vs                       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%      2 points algorithm with IMU         %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%        planar scene : homography         %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name : Meng Di 
% Date ?19/01/2018
%
%
%

close all
clear all

%test 1 : example with a particular data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    data generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lowerbound=-300;
upperbound=300;
nbpoints = 50
zposition = 500
P=[lowerbound + (upperbound-lowerbound).*rand(nbpoints,1),...
    lowerbound + (upperbound-lowerbound).*rand(nbpoints,1),...
    zposition*ones(nbpoints,1),...
    ones(nbpoints,1)];

% the data belong on a plane  of equation N'X+d=0 d=-zposition N = [0,0,1]'

%camera parameter (camera is calibrated)
f=1; u0 = 0; v0 = 0;
K=[f 0 u0;0 f v0;0 0 1];
K1=[f 0 u0 0;0 f v0 0;0 0 1 0];

%%%%%%%%%%%%%%%%%%%%%%%    camera 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% camera position at time 1
tx1=10; ty1=-2; tz1=20;
% rotation angles in degree
roll1=5; pitch1=5; yaw1=30;

% rotation of the camera 1

Rp1=[cos(deg2rad(pitch1)) 0 sin(deg2rad(pitch1));...
    0 1 0;...
    -sin(deg2rad(pitch1)) 0 cos(deg2rad(pitch1))];
Ry1=[cos(deg2rad(yaw1)) -sin(deg2rad(yaw1)) 0;...
    sin(deg2rad(yaw1)) cos(deg2rad(yaw1)) 0;...
    0 0 1];
Rr1=[1 0 0;...
    0 cos(deg2rad(roll1)) -sin(deg2rad(roll1));...
    0 sin(deg2rad(roll1)) cos(deg2rad(roll1))];
R1=Ry1*Rp1*Rr1;

T1 = [tx1,ty1,tz1]';

% camera image 1 :
for i =1 : nbpoints
    P1(i,:) = K1*[R1' , -R1'*T1; 0 0 0 1]*P(i,:)';
end

%%%%%%%%%%%%%%%%%%%%%%%    camera 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% camera position at time 2
tx2=10; ty2=-6; tz2=10;
% rotation angles in degree
roll2=3; pitch2=5; yaw2=5;

% rotation of the camera 2

Rp2=[cos(deg2rad(pitch2)) 0 sin(deg2rad(pitch2));...
    0 1 0;...
    -sin(deg2rad(pitch2)) 0 cos(deg2rad(pitch2))];
Ry2=[cos(deg2rad(yaw2)) -sin(deg2rad(yaw2)) 0;...
    sin(deg2rad(yaw2)) cos(deg2rad(yaw2)) 0;...
    0 0 1];
Rr2=[1 0 0;...
    0 cos(deg2rad(roll2)) -sin(deg2rad(roll2));...
    0 sin(deg2rad(roll2)) cos(deg2rad(roll2))];
R2 =Ry2*Rp2*Rr2;

T2 = [tx2,ty2,tz2]';

% camera image 2 :

for i =1 : nbpoints
    P2(i,:) = K1*[R2' , -R2'*T2; 0 0 0 1]*P(i,:)';
end

%%%%%%%%%%%%%%%%%%%%%%%    display data     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% hold on
% plot3(P(:,1),P(:,2),P(:,3),'*')
% c1 = plot3(tx1,ty1,tz1,'g*');
% Rx = R1*[100,0,0]';
% line([tx1,tx1+Rx(1)],[ty1,ty1+Rx(2)],[tz1,tz1+Rx(3)],'Color','b')
% Ry = R1*[0,100,0]';
% line([tx1,tx1+Ry(1)],[ty1,ty1+Ry(2)],[tz1,tz1+Ry(3)],'Color','g')
% Rz = R1*[0,0,100]';
% line([tx1,tx1+Rz(1)],[ty1,ty1+Rz(2)],[tz1,tz1+Rz(3)],'Color','r')
% c2 = plot3(tx2,ty2,tz2,'r*');
% Rx = R2*[100,0,0]';
% line([tx2,tx2+Rx(1)],[ty2,ty2+Rx(2)],[tz2,tz2+Rx(3)],'Color','b')
% Ry = R2*[0,100,0]';
% line([tx2,tx2+Ry(1)],[ty2,ty2+Ry(2)],[tz2,tz2+Ry(3)],'Color','g')
% Rz = R2*[0,0,100]';
% line([tx2,tx2+Rz(1)],[ty2,ty2+Rz(2)],[tz2,tz2+Rz(3)],'Color','r')
% axis equal
% legend([c1, c2], 'Camera 1', 'Camera 2')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    4 pts algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Theoretical Homography H :

N = [0,0,1]';
d = (N'*T1-zposition); 
% explain why d is expressed like this!
% -> The distance d is computed by substraction of 
% -> camera's height and points' height 


Rt = R2'*R1;
Tt = R2'*(T1-T2);
H = Rt-Tt*(R1'*N)'/d;
H=H/H(3,3)
HP1 = H*P1(1,:)';
HP1 = HP1/HP1(3)
GroundTruth = (P2(1,:)/P2(1,3))'


% Homography estimation
H4pt=homography2d(P1',P2');
H4pt=H4pt/H4pt(3,3)

% Homography decomposition
% solutions = invhomog(H4pt);
% 
% solutions(1).T(:,4)/solutions(1).T(3,4);
% solutions(1).T(1:3,1:3)
% 
% solutions(2).T(:,4)/solutions(2).T(3,4);
% solutions(2).T(1:3,1:3);



% HYPOTHESIS roll and pitch angles are known :

% Virtual image from camera 2 (Z axis correpsonds to the vertical)
PV1 = (Rp1*Rr1*P1(:,1:3)')';
PV2 = (Rp2*Rr2*P2(:,1:3)')';

% 2 points algorithm
N =[0 0 1]';

HV2 = Ry2'*Ry1-Ry2'*(T1-T2)/(N'*T1-zposition)*(Ry1'*N)'

HV2*PV1(2,:)'/norm(HV2*PV1(2,:)')
GroundTruth = (PV2(2,:)/norm(PV2(2,:)))'
TrueHomography = HV2/HV2(3,3)

H = homography2d2Pt(PV1',PV2');

EstimatedH=H/H(3,3)




%%%%%%%%%%%%%%%%Yaw Estimation (see exercice 1)

Yaw= atan2(H(2,1),H(1,1))
YawGroundTruth=deg2rad(yaw1-yaw2)

Azimuth=atan2(H(2,3),H(1,3));
SEalpha=-((sin(Yaw)/H(2,1))-1);
CEalpha=(-sin(Yaw)/H(2,1)*H(1,3))/cos(Azimuth);
Elevation=atan2(SEalpha,CEalpha);
T2E = [cos(Azimuth)*cos(Elevation);sin(Azimuth)*cos(Elevation);sin(Elevation)];


T2t=Ry2'*(T1-T2);
T2t/T2t(3)
T2E/T2E(3)



%***********************************************************************
%**************test 2 : example with different datas********************
%***********************************************************************
% propose a test with different camera position
% R1 = I, T1 = 0
% angles of rotation of camera 2 between 0° and 45°
% translation of 0 to 100

% camera position at time 3
tx3=18; ty3=-8; tz3=80;
% rotation angles in degree
roll3=5; pitch3=5; yaw3=8;

% rotation of the camera 2

Rp3=[cos(deg2rad(pitch3)) 0 sin(deg2rad(pitch3));...
    0 1 0;...
    -sin(deg2rad(pitch3)) 0 cos(deg2rad(pitch3))];
Ry3=[cos(deg2rad(yaw3)) -sin(deg2rad(yaw3)) 0;...
    sin(deg2rad(yaw3)) cos(deg2rad(yaw3)) 0;...
    0 0 1];
Rr3=[1 0 0;...
    0 cos(deg2rad(roll3)) -sin(deg2rad(roll3));...
    0 sin(deg2rad(roll3)) cos(deg2rad(roll3))];
R3 =Ry3*Rp3*Rr3;

T3 = [tx3,ty3,tz3]';

% camera image 3 :

for i =1 : nbpoints
    P3(i,:) = K1*[R3' , -R3'*T3; 0 0 0 1]*P(i,:)';
end

% Displaying the data
% figure
% hold on
% plot3(P(:,1),P(:,2),P(:,3),'*')
% c1 = plot3(tx1,ty1,tz1,'g*');
% Rx = R1*[100,0,0]';
% line([tx1,tx1+Rx(1)],[ty1,ty1+Rx(2)],[tz1,tz1+Rx(3)],'Color','b')
% Ry = R1*[0,100,0]';
% line([tx1,tx1+Ry(1)],[ty1,ty1+Ry(2)],[tz1,tz1+Ry(3)],'Color','g')
% Rz = R1*[0,0,100]';
% line([tx1,tx1+Rz(1)],[ty1,ty1+Rz(2)],[tz1,tz1+Rz(3)],'Color','r')
% c3 = plot3(tx3,ty3,tz3,'r*');
% Rx = R3*[100,0,0]';
% line([tx3,tx3+Rx(1)],[ty3,ty3+Rx(2)],[tz3,tz3+Rx(3)],'Color','b')
% Ry = R3*[0,100,0]';
% line([tx3,tx3+Ry(1)],[ty3,ty3+Ry(2)],[tz3,tz3+Ry(3)],'Color','g')
% Rz = R3*[0,0,100]';
% line([tx3,tx3+Rz(1)],[ty3,ty3+Rz(2)],[tz3,tz3+Rz(3)],'Color','r')
% axis equal
% legend([c1, c3], 'Camera 1', 'Camera 2')

%************************4pt algorithm for test 2***********************
% Theoretical Homography H :

N = [0,0,1]';
d2 = (N'*T1-zposition); 
% explain why d is expressed like this!
% -> The distance d is computed by substraction of 
% -> camera's height and points' height 


Rt = R3'*R1;
Tt = R3'*(T1-T3);
H = Rt-Tt*(R1'*N)'/d2;
H=H/H(3,3)
HP1 = H*P1(1,:)';
HP1 = HP1/HP1(3)
GroundTruth = (P3(1,:)/P3(1,3))'


% Homography estimation
H4pt=homography2d(P1',P3');
H4pt=H4pt/H4pt(3,3)
%***********************************************************************

%************************2pt algorithm for test 2***********************

% Virtual image from camera 2 (Z axis correpsonds to the vertical)
PV1 = (Rp1*Rr1*P1(:,1:3)')';
PV3 = (Rp3*Rr3*P3(:,1:3)')';

% 2 points algorithm
N =[0 0 1]';

HV3 = Ry3'*Ry1-Ry3'*(T1-T3)/(N'*T1-zposition)*(Ry1'*N)'

HV3*PV1(2,:)'/norm(HV3*PV1(2,:)')
GroundTruth = (PV3(2,:)/norm(PV3(2,:)))'
TrueHomography = HV3/HV3(3,3)

H = homography2d2Pt(PV1',PV3');

EstimatedH=H/H(3,3)

%Yaw Estimation

Yaw= atan2(H(2,1),H(1,1))
YawGroundTruth=deg2rad(yaw1-yaw3)

Azimuth=atan2(H(2,3),H(1,3));
SEalpha=-((sin(Yaw)/H(2,1))-1);
CEalpha=(-sin(Yaw)/H(2,1)*H(1,3))/cos(Azimuth);
Elevation=atan2(SEalpha,CEalpha);
T3E = [cos(Azimuth)*cos(Elevation);sin(Azimuth)*cos(Elevation);sin(Elevation)];


T3t=Ry3'*(T1-T3);
T3t/T3t(3)
T3E/T3E(3)


%***********************************************************************

%***********************************************************************
%********************test 3 : example with noise************************
%***********************************************************************
% propose a test with different camera position
% R1 = I, T1 = 0
% angles of rotation of camera 2 between 0° and 45°
% translation of 0  to 100
% AND white noise in image points of camera 2 between 0 to 1 pixel std
% We continue using the camera 2 from tets 2

% Adding white noise to image points of camera 2
% White noise between 0 to 1 std
P3_n = P3 + rand(size(P3));

%************************4pt algorithm for test 3***********************
% Theoretical Homography H :

N = [0,0,1]';
d2 = (N'*T1-zposition); 

Rt = R3'*R1;
Tt = R3'*(T1-T3);
H = Rt-Tt*(R1'*N)'/d2;
H=H/H(3,3)
HP1 = H*P1(1,:)';
HP1 = HP1/HP1(3)
GroundTruth = (P3_n(1,:)/P3_n(1,3))'


% Homography estimation
H4pt=homography2d(P1',P3_n');
H4pt=H4pt/H4pt(3,3)
%***********************************************************************

%************************2pt algorithm for test 2***********************

% Virtual image from camera 2 (Z axis correpsonds to the vertical)
PV1 = (Rp1*Rr1*P1(:,1:3)')';
PV3_n = (Rp3*Rr3*P3_n(:,1:3)')';

% 2 points algorithm
N =[0 0 1]';

HV3 = Ry3'*Ry1-Ry3'*(T1-T3)/(N'*T1-zposition)*(Ry1'*N)'

HV3*PV1(2,:)'/norm(HV3*PV1(2,:)')
GroundTruth = (PV3_n(2,:)/norm(PV3_n(2,:)))'
TrueHomography = HV3/HV3(3,3)

H = homography2d2Pt(PV1',PV3_n');

EstimatedH=H/H(3,3)
%***********************************************************************

%************************ransac algorithm for test 2***********************

[H, ~ ] = ransacfithomography(P1', P3_n', 0.000001);
H1 = H./H(3,3)
[H2pt, ~ ] = ransacfithomography(PV1', PV3_n', 0.000001);
H2pt = H2pt./H2pt(3,3)
%***********************************************************************

%***********************************************************************
%***********test 4 : example with noise on IMU information**************
%***********************************************************************
% propose a test with different camera position
% R1 = I, T1 = 0
% angles of rotation of camera 2 between 0° and 45°
% translation of 0 to 100
% AND white noise in image points of camera 2 between 0 to 1 pixel std
% AND white noise in IMU 2 between 0 to 2°

% camera position at time 3 (keep the same)
tx3=18; ty3=-8; tz3=80;
% rotation angles in degree (add noise)
roll4 = roll3 + 2*rand(); 
pitch4 = pitch3 + 2*rand(); 
yaw4 = yaw3 + 2*rand();

% rotation of the camera 2

Rp4=[cos(deg2rad(pitch3)) 0 sin(deg2rad(pitch3));...
    0 1 0;...
    -sin(deg2rad(pitch3)) 0 cos(deg2rad(pitch3))];
Ry4=[cos(deg2rad(yaw3)) -sin(deg2rad(yaw3)) 0;...
    sin(deg2rad(yaw3)) cos(deg2rad(yaw3)) 0;...
    0 0 1];
Rr4=[1 0 0;...
    0 cos(deg2rad(roll3)) -sin(deg2rad(roll3));...
    0 sin(deg2rad(roll3)) cos(deg2rad(roll3))];
R4 =Ry4*Rp4*Rr4;

T4 = [tx3,ty3,tz3]';

% camera image 3 :

for i =1 : nbpoints
    P4(i,:) = K1*[R4' , -R4'*T4; 0 0 0 1]*P(i,:)';
end

% add noise to image points
P4_n = P4 + rand(size(P4));

%************************4pt algorithm for test 4***********************
% Theoretical Homography H :

N = [0,0,1]';
d2 = (N'*T1-zposition); 

Rt = R4'*R1;
Tt = R4'*(T1-T4);
H = Rt-Tt*(R1'*N)'/d2;
H=H/H(3,3)
HP1 = H*P1(1,:)';
HP1 = HP1/HP1(3)
GroundTruth = (P4_n(1,:)/P4_n(1,3))'


% Homography estimation
H4pt=homography2d(P1',P4_n');
H4pt=H4pt/H4pt(3,3)
%***********************************************************************

%************************2pt algorithm for test 2***********************

% Virtual image from camera 2 (Z axis correpsonds to the vertical)
PV1 = (Rp1*Rr1*P1(:,1:3)')';
PV4_n = (Rp4*Rr4*P4_n(:,1:3)')';

% 2 points algorithm
N =[0 0 1]';

HV4 = Ry4'*Ry1-Ry4'*(T1-T4)/(N'*T1-zposition)*(Ry1'*N)'

HV4*PV1(2,:)'/norm(HV4*PV1(2,:)')
GroundTruth = (PV4_n(2,:)/norm(PV4_n(2,:)))'
TrueHomography = HV4/HV4(3,3)

H = homography2d2Pt(PV1',PV4_n');

EstimatedH=H/H(3,3)
%***********************************************************************

%************************ransac algorithm for test 4***********************

[H, ~ ] = ransacfithomography(P1', P4_n', 0.000001);
H1 = H./H(3,3)
[H2pt, ~ ] = ransacfithomography(PV1', PV4_n', 0.000001);
H2pt = H2pt./H2pt(3,3)
%***********************************************************************
