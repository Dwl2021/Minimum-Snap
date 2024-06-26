clc;clear;close all;
numPoints = input('Enter the number of points you want to input: ');

while ~isnumeric(numPoints) || mod(numPoints, 1) ~= 0 || numPoints < 1
    disp('Please enter a positive integer.');
    numPoints = input('Enter the number of points you want to input: ');
end

% path = ginput(numPoints) * 100.0;
figure('Position', [500, 500, 1000, 1000]);
hold on;
axis([0 100 0 100]);
path = [];
for i=1:numPoints
    [x, y] = ginput(1);
    path = [[x,y];path];
    plot(x, y, 'o', 'MarkerSize', 20, 'MarkerEdgeColor', '#FFA500', 'MarkerFaceColor', '#FFA500');
end
hold on;

n_order       = 7;% order of poly
n_seg         = size(path,1)-1;% segment number
n_poly_perseg = (n_order+1); % coef number of perseg

ts = zeros(n_seg, 1);
% calculate time distribution in proportion to distance between 2 points
% dist     = zeros(n_seg, 1);
% dist_sum = 0;
% T        = 25;
% t_sum    = 0;
% 
% for i = 1:n_seg
%     dist(i) = sqrt((path(i+1, 1)-path(i, 1))^2 + (path(i+1, 2) - path(i, 2))^2);
%     dist_sum = dist_sum+dist(i);
% end
% for i = 1:n_seg-1
%     ts(i) = dist(i)/dist_sum*T;
%     t_sum = t_sum+ts(i);
% end
% ts(n_seg) = T - t_sum;

% or you can simply set all time distribution as 1
for i = 1:n_seg
    ts(i) = 1.0;
end

poly_coef_x = MinimumSnapQPSolver(path(:, 1), ts, n_seg, n_order);
poly_coef_y = MinimumSnapQPSolver(path(:, 2), ts, n_seg, n_order);


% display the trajectory
X_n = [];
Y_n = [];
k = 1;
tstep = 0.01;
for i=0:n_seg-1
    %#####################################################
    % STEP 3: get the coefficients of i-th segment of both x-axis
    % and y-axis
    Pxi = poly_coef_x((n_order+1)*i+1:(n_order+1)*i+n_order+1); 
    Pyi = poly_coef_y((n_order+1)*i+1:(n_order+1)*i+n_order+1);
    for t = 0:tstep:ts(i+1)
        %flip表示翻转Pxi中的元素
        %polyval（[P0,P1,P2],t）生成P0t^2+P1^t+P2,因此需要翻转原始数据
        X_n(k)  = polyval(flip(Pxi), t);
        Y_n(k)  = polyval(flip(Pyi), t);
        k = k + 1;
    end
end
 
plot(X_n, Y_n ,'Color','#DC143C','LineWidth',5);
hold on
scatter(path(1:size(path, 1), 1), path(1:size(path, 1), 2));

function poly_coef = MinimumSnapQPSolver(waypoints, ts, n_seg, n_order)
    start_cond = [waypoints(1), 0, 0, 0];
    end_cond   = [waypoints(end), 0, 0, 0];
    %#####################################################
    % STEP 1: compute Q of p'Qp
    Q = getQ(n_seg, n_order, ts);
    %#####################################################
    % STEP 2: compute Aeq and beq 
    [Aeq, beq] = getAbeq(n_seg, n_order, waypoints, ts, start_cond, end_cond);
    f = zeros(size(Q,1),1);
    poly_coef = quadprog(Q,f,[],[],Aeq, beq);
end