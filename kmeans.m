%
% Copyright (c) 2018, Vishal_S
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Title: K-means clustering (As part of BT3041)
% 
% Developer: Vishal S
% 
% Contact Info: vishalsubbu97@gmail.com
%

clear all
clf
clc

%% Load data 
clear all
clf
clc
hold on
filename = 'outlier.txt';
% read from the file
fileID = fopen(filename,'r');
formatSpec = '%f,%f,%f';
sizeA = [3 Inf];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);

%% Variable 
no_of_clusters = 4;
convergence = 1.0;
iteration = 1;

% Define variables of appropriate length
max_len = size(A,2);
mean_old = zeros(2,no_of_clusters);
mean_new = zeros(2,no_of_clusters);
cluster_id  = zeros(max_len,1);
sum  = zeros(2,no_of_clusters);
no_of_points = zeros(1,no_of_clusters);
initial_mean = randi([1 max_len],[no_of_clusters 1]);
%Define stopping criteria
sse_old = 0.0;
convergence_limit = 0.001;
max_iteration = 1000;

%% Kmeans algorithm
% initialising the clusters' means
for i = 1:no_of_clusters
    for j = 1:2
    mean_old(j,i) = A(j,initial_mean(i));
    end
end

while ((convergence > convergence_limit)&&(iteration < max_iteration))
    X = ['Iteration = ',num2str(iteration)];
    disp(X);
    % Assigning each point to a cluster
    dist = 0.0;
    for i = 1:max_len
        dist_x = (A(1,i)-mean_old(1,1))^2;
        dist_y = (A(2,i)-mean_old(2,1))^2;
        dist = dist_x + dist_y;
        index = 1;
        for j = 2:no_of_clusters
            dist_x = (A(1,i)-mean_old(1,j))^2;
            dist_y = (A(2,i)-mean_old(2,j))^2;
            if(dist > dist_x + dist_y)
                dist = dist_x + dist_y;
                index = j;
            end
        end
        cluster_id(i,1) = index;
    end
    
    for i = 1:no_of_clusters
        sum(1,i) = 0.0;
        sum(2,i) = 0.0;
        no_of_points(1,i) = 0.0;
    end
    % Calculating new mean
    for i = 1:max_len
        sum(1,i) = 0.0;
        sum(2,i) = 0.0;
        sum(1,cluster_id(i,1)) = sum(1,cluster_id(i,1)) +  A(1,i);
        sum(2,cluster_id(i,1)) = sum(2,cluster_id(i,1)) +  A(2,i);
        no_of_points(1,cluster_id(i,1)) =no_of_points(1,cluster_id(i,1)) + 1;
    end
    
    % Calculating the new means
    for i = 1:no_of_clusters
        for j = 1:2
            mean_new(j,i) = sum(j,i)/no_of_points(1,i);
        end
    end
    
    % Caclulating SSE
    sse_new = 0.0;
    for i = 1:max_len
        for j = 1:2
            sse_new = sse_new + ( A(j,i)-mean_new(j,cluster_id(i)))^2;
        end
    end
    convergence = abs((sse_new - sse_old)) / sse_new ;
    % update the variables
    iteration = iteration + 1;
    sse_old = sse_new;
    mean_old = mean_new;
end
%% Plot results
k=max(cluster_id);
Colors=hsv(k);
Legends = {};
for i=0:k
    A_i=A(:,cluster_id==i);
    if i~=0
        Style = 'x';
        MarkerSize = 8;
        Color = Colors(i,:);
        Legends{end+1} = ['Cluster #' num2str(i)];
    end
    if ~isempty(A_i)
        %scatter3(A_i(1,:),A_i(2,:),A_i(3,:),Style,'MarkerSize',MarkerSize,'MarkerFaceColor',Color);
        scatter(A_i(1,:),A_i(2,:),Style,'MarkerFaceColor',Color);
    end
    hold on;
end
Style = 'o';
MarkerSize = 10;
Color = [0 0 0];
Legends{end+1} = 'Centroids';
scatter(mean_new(1,:),mean_new(2,:),Style,'MarkerFaceColor',Color);      

hold off;
grid on;
legend(Legends);
legend('Location', 'NorthEastOutside');


%% Print Results 

X = ['SSE =',num2str(sse_new)];
disp(X);
Y = 'Centroids';
disp(Y);
for i = 1:no_of_clusters
    Y = ['      Centroid #',num2str(i),'=[',num2str(mean_new(1,i)),' ',num2str(mean_new(2,i)),']'];
    disp(Y)
end
