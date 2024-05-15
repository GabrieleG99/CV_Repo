close all;
clear all;
clc;

img = imread('../data/PalazzoTe.jpg');

img=rgb2gray(img);
img=histeq(img);
img=imrotate(img, -90);
load("a.mat");
load("b.mat");
load("c.mat");
load("d.mat");
load("C1.mat");
load("C2.mat");
load("l1.mat");
load("l2.mat");
%% load the  two conics coefficients into vectors 

[a1, b1, c1, d1, e1, f1] = deal(C1(1),C1(2)*2,C1(5),C1(3)*2,C1(6)*2,C1(9));
[a2, b2, c2, d2, e2, f2] = deal(C2(1),C2(2)*2,C2(5),C2(3)*2,C2(6)*2,C2(9));

%% intersect ellipses
syms x y;

eq1 = a1*x^2 + b1*x*y + c1*y^2 + d1*x + e1*y + f1;
eq2 = a2*x^2 + b2*x*y + c2*y^2 + d2*x + e2*y + f2;

eqns = [eq1 ==0, eq2 ==0];
S = solve(eqns, [x,y]);
%% solutions of the intersection

s1 = [double(S.x(1));double(S.y(1));1];
s2 = [double(S.x(2));double(S.y(2));1];
s3 = [double(S.x(3));double(S.y(3));1];
s4 = [double(S.x(4));double(S.y(4));1];

% s1 and s2 are complex conjugate, so they identify the image of the circular 
% points of C1 and C2

II = s1;
JJ = s2;
horizon=hcross(II,JJ); %compute the horizon
horizon = horizon./norm(horizon);
img_dual_con = II*JJ.' + JJ*II.'; % image of conic dual to circular points
img_dual_con = img_dual_con./norm(img_dual_con);
%% draw vanishing line associated to the plane orthogonal to the cylinder axis
figure;
imshow(img)
hline(horizon)
hold all
%% compute the rectifying homography
[U,D,U_T] = svd(img_dual_con)

% in order to take account of the numerical error a correction procedure
% has to be implemented. It consists in taking the square root of the
% matrix D and multiply it with the left singolar value because when we
% execute the SVD on the image of the conic dual to the circular points we
% obtain a decomposition of type:
%
%                      |s_1 0 0|     |√s_1 0  0||1 0 0||√s_1 0  0|
%       SVD(imDCCP) = U|0 s_2 0|U^T=U|0 √s_2  0||0 1 0||0 √s_2  0|U^T
%                      |0  0  0|     |0   0   1||0 0 0||0   0   1|
%

D(3,3) = 1;
A = sqrt(D)
H_R = inv(U * A);
H = inv(H_R);
%% compute image projection of cylinder axis a
horizon=horizon./horizon(3);
c1_center=inv(C1)*horizon; %compute center of c1
c2_center=inv(C2)*horizon; %compute center of c2

lc1c2=hcross(c1_center,c2_center); % compute cylinder axis a

c1_center=c1_center./c1_center(3);
c2_center=c2_center./c2_center(3);

plot([c1_center(1), c2_center( 1)], [c1_center(2), c2_center(2)], 'b');
scatter(c1_center(1),c1_center(2),'filled',MarkerFaceColor='b')
scatter(c2_center(1),c2_center(2),'filled',MarkerFaceColor='b')
text(c1_center(1), c1_center(2), 'c1', 'FontSize', 20, 'Color', 'r')
text(c2_center(1), c2_center(2), 'c2', 'FontSize', 20, 'Color', 'r')

% compute and plot the vanishing point of the cylinder axis
V = hcross(l1, l2);
V = V./V(3);
plot(V(1), V(2), 'b');
scatter(V(1),V(2),'filled',MarkerFaceColor='b');
text(V(1), V(2), 'V','FontSize', 20, 'Color', 'r');

h1 = H(:, 1);
h2 = H(:, 2);
%% computation of the calibration matrix K
syms w1 w2 w3 w4;

W = [w1 0 w2;
    0 1 w3;
    w2 w3 w4];

X = hcross(horizon, (W*V));

const_1 = h1.'*W*h2;
const_2 = h1.'*W*h1 - h2.'*W*h2;
const_3 = X(1,1);
const_4 = X(2,1);

S_w = solve([const_1==0, const_2==0, const_3==0, const_4==0], [w1,w2,w3,w4]);

w = double([S_w.w1 0 S_w.w2;
    0  1 S_w.w3;
    S_w.w2 S_w.w3 S_w.w4])

diac = inv(w);
K = chol(diac);
K=K./K(9);

%% computation of the homography wrt the axis plane

lab = hcross(a, b);
lcd = hcross(c, d);

v1 = hcross(lab, lcd);
v2 = hcross(l1, l2);

l_inf_axis = hcross(v1, v2) %vanishing line of the axis plane

syms x y;

X = [x; y; 1];

eqn1 = X.'*w*X;
eqn2 = l_inf_axis.'*X;

S2 = solve([eqn1==0, eqn2==0], [x, y]);

I_a = double([S2.x(1); S2.y(1); 1]);
J_a = double([S2.x(2); S2.y(2); 1]);

imDCCP = I_a*J_a.'+J_a*I_a.';
imDCCP = imDCCP./norm(imDCCP);

[u, dd, v] = svd(imDCCP);
dd(3,3) = 1;
AA = sqrt(dd);
H_R_a = inv(u*AA);
H_a = inv(H_R_a);

%% compute orientation of the camera

Q = inv(K)*H_a

R_t = [Q(:, 1), Q(:, 2), cross(Q(:, 1), Q(:, 2)), Q(:, 3);
    0 0 0 1]

%% compute cylinder axis direction in camera reference
c1_3d=R_t*[c1_center(1);c1_center(2);0;c1_center(3)]; % compute c1 center in camera reference
c2_3d=R_t*[c2_center(1);c2_center(2);0;c2_center(3)]; % compute c1 center in camera reference

c1_3d=c1_3d./c1_3d(4)
c2_3d=c2_3d./c2_3d(4)

figure(2);
scatter3(c1_3d(1),c1_3d(2),c1_3d(3),'filled')
text(c1_3d(1),c1_3d(2),c1_3d(3),'C1_center')

hold on

scatter3(c2_3d(1),c2_3d(2),c2_3d(3),'filled')
text(c2_3d(1),c2_3d(2),c2_3d(3),'C2_center')

scatter3(0,0,0,'r')
text(0,0,0,'camera')

plot3([c1_3d(1) c2_3d(1)],[c1_3d(2) c2_3d(2)],[c1_3d(3) c2_3d(3)],'color','b','LineWidth',2)
xlabel('X');ylabel("Y");zlabel("Z")

%% radius to distance ratio

a_3d = R_t*[a(1); a(2); 0; a(3)];
a_3d = a_3d./a_3d(4);

dist = norm(c1_3d - c2_3d);
rad = norm(a_3d-c1_3d);

ratio = rad / dist;







