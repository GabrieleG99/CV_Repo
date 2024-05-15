%% rectification stratified on axis plane

lab = hcross(a, b);
lcd = hcross(c, d);

v1 = hcross(lab, lcd);
v2 = hcross(l1, l2);

l_inf_axis = hcross(v1, v2)

figure(5);
hold all;
plot(v1(1), v1(2), 'r.', 'MarkerSize',50);
plot(v2(1), v2(2), 'r.', 'MarkerSize',50);
imshow(img)
line([v1(1), v2(1)], [v1(2), v2(2)], 'Color', 'Green', 'Linewidth', 3);

H_aff = [eye(2), zeros(2,1); l_inf_axis(:)']

fprintf('The vanish line is mapped to:\n');
disp(inv(H)'*l_inf_axis);

tform = projective2d(H');
J = imwarp(img, tform);
figure;
imshow(J);
imwrite(J, '../data/affRect.JPG')

%%
affImg = imread('../data/hist_gray_img.jpg');
figure; imshow(affImg);
numConstarints = 5;

hold all;
count = 1;
A = zeros(numConstarints, 6);

while (count <= numConstarints)
    figure(gcf);
    col = 'rgbcmykwrgbcmykw';
    segment1 = drawline('Color', col(count));
    segment2 = drawline('Color', col(count));

    l = segToLine(segment1.Position);
    m = segToLine(segment2.Position);

    A(count, :)= [l(1)*m(1), 0.5*(l(1)*m(2)+l(2)*m(1)), l(2)*m(2),...
        0.5*(l(1)*m(3)+l(3)*m(1)), 0.5*(l(2)*m(3)+l(3)*m(2)), l(3)*m(3)];
    count = count + 1;
end

[~,~,v] = svd(A);
sol = v(:, end);
imDCCP2 = [sol(1), sol(2)/2, sol(4)/2;
    sol(2)/2 sol(3) sol(5)/2;
    sol(4)/2 sol(5)/2 sol(6)];

[U,D,V] = svd(imDCCP2);
D(3,3) = 1;
A = U * sqrt(D);
C = [eye(2), zeros(2,1);zeros(1,3)];

min(norm(A*C*A' - imDCCP2), norm(A*C*A' + imDCCP2))

H_rect = inv(A);

tform = projective2d(H_rect');
IM = imwarp(affImg, tform);
figure; imshow(IM);