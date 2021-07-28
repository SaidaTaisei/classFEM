% 正方形物体の動的な弾塑性変形 (10&times;10)
% g, cm, msec

width = 30; height = 30; h = 1;
m = 10; n = 10;

E = 10.0; nu = 0.48;
c1 = 0.04*E;
c2 = 2000*E;
%c2 = 1200*E; % ytop = 25
rho = 1.00;

tp = 1e+3;
vp = 0.8*(height/3)/tp;
th = 1e+3;
%tf = 3e+3;
tf = 5e+3; % c2 = 2000*E のとき

alpha = 1e+6;

npoint = m*n;
[point,triangle] = rectangular_object(m,n,width,height);
         
[lambda,mu] = Lame_constants(E, nu);
[lambdav1,muv1] = Lame_constants(c1, nu);
[lambdav2,muv2] = Lame_constants(c2, nu);

[Jlambda, Jmu] = connection_matrices(point,triangle,h);
M = inertia_matrix(point,triangle,h,rho);

% 床に固定する節点
A_fixed = zeros(npoint*2,2*m);
for k=1:2*m; A_fixed(k,k)=1; end

% 押す節点
A_push = zeros(npoint*2,8); % P94 - P97
for k=1:8; A_push(94*2-1+(k-1),k) = 1; end

push = @(t) push_param (4, t, tp,vp,th);
dotpush = @(t) dotpush_param (4, t, tp,vp,th);

square_object_push = @(t,q) three_element_square_object_push_param(t,q,npoint,M, ...
    Jlambda,Jmu,lambda,mu,lambdav1,muv1,lambdav2,muv2, ...
    4,A_push,push,dotpush,m,A_fixed,alpha);
square_object_free = @(t,q) three_element_square_object_free_param(t,q,npoint,M, ...
    Jlambda,Jmu,lambda,mu,lambdav1,muv1,lambdav2,muv2, ...
    m,A_fixed,alpha);

% 上部を押している
interval = [0, tp+th];
qinit = zeros(npoint*8,1);
[time_push, q_push] = ode23tb(square_object_push, interval, qinit);

% 上部を解放した後
interval = [tp+th, tp+th+tf];
sz = size(q_push); tn = sz(1);
qinit = q_push(tn,:);
[time_free, q_free] = ode23tb(square_object_free, interval, qinit);

time = [time_push; time_free];
q = [q_push; q_free];

set(0,'defaultAxesFontSize',16);
set(0,'defaultTextFontSize',16);

figure('position', [0, 0, 600, 600]);
clf;
for t = 0:0.1*1000:tp+th+tf
    fprintf("時刻 %f\n", t);
    index = nearest_index(time, t);
    draw_triangles(point, triangle, q(index,1:npoint*2));
    hold off;
    xlim([-10,40]);
    ylim([-10,40]);
    pbaspect([1 1 1]);
    grid on;
    filename = strcat('three-element-10-10/deform-', num2str(floor(t),'%04d'), '.png');
    saveas(gcf, filename, 'png');
end

clf('reset');
ts = time(1);
te = time(end);
fr = 1;
clear M;
for t = 0:0.1*100:tp+th+tf
    index = nearest_index(time, t);
    draw_triangles(point, triangle, q(index,1:npoint*2));
    hold off;
    xlim([-10,40]);
    ylim([-10,40]);
    xticks([-10:10:40])
    yticks([-10:10:40])
    pbaspect([1 1 1]);
    title(['time ' num2str(t/1000,"%3.2f")]);
    grid on;
    drawnow;
    M(fr) = getframe(gcf);
    fr = fr + 1;
    disp(t/1000);
end
M(fr) = getframe(gcf);

v = VideoWriter('three_element_square_object_10_10', 'MPEG-4');
open(v);
writeVideo(v, M);
close(v);