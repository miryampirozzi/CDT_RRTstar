close all
clear all
clear figure

mappp.points=[];
mappp.lines=[];

% mapp = import_data_catia_vrml('warehouse.wrl');
% mapp = import_data_catia_vrml('warehouse2.wrl');
mapp = import_data_catia_vrml('nuclearFacility.wrl');

%% Warehouse
% for i=1:size(mapp.points,2)
%     if mapp.points(1,i)==-0.02
%         mappp.points=[mappp.points; mapp.points(2,i) mapp.points(3,i)];
%     end
% end

%% Warehouse2

% for i=1:size(mapp.points,2)
%     if mapp.points(1,i)==-0.02
%         mappp.points=[mappp.points; mapp.points(2,i) mapp.points(3,i)]; % to put in vertical
%     end
% end

%% Nuclear facility
for i=1:size(mapp.points,2)
    if mapp.points(3,i)==0
        mappp.points=[mappp.points; mapp.points(1,i) mapp.points(2,i)]; % to put in vertical
    end
end


map.points=unique(mappp.points,'rows');
map.points=map.points*2000;

%% Warehouse
% for i=1:size(mapp.lines,2)
%     if mapp.points(1,mapp.lines(1,i))==-0.02 && mapp.points(1,mapp.lines(2,i))==-0.02
%         x1=mapp.points(2,mapp.lines(1,i)); 
%         y1=mapp.points(3,mapp.lines(1,i));
%         x2=mapp.points(2,mapp.lines(2,i)); 
%         y2=mapp.points(3,mapp.lines(2,i));
%         for j=1:size(map.points,1)
%             if x1*2000==map.points(j,1) && y1*2000==map.points(j,2)
%                 jj=j;
%             end
%         end
%         for z=1:size(map.points,1)
%             if x2*2000==map.points(z,1) && y2*2000==map.points(z,2)
%                 zz=z;
%             end
%         end
%         mappp.lines=[mappp.lines; jj zz];
% 
%     end
% end

%% Warehouse2 
% for i=1:size(mapp.lines,2)
%     if mapp.points(1,mapp.lines(1,i))==-0.02 && mapp.points(1,mapp.lines(2,i))==-0.02
%         x1=mapp.points(2,mapp.lines(1,i)); 
%         y1=mapp.points(3,mapp.lines(1,i));
%         x2=mapp.points(2,mapp.lines(2,i)); 
%         y2=mapp.points(3,mapp.lines(2,i));
%         for j=1:size(map.points,1)
%             if x1*2000==map.points(j,1) && y1*2000==map.points(j,2)
%                 jj=j;
%             end
%         end
%         for z=1:size(map.points,1)
%             if x2*2000==map.points(z,1) && y2*2000==map.points(z,2)
%                 zz=z;
%             end
%         end
%         mappp.lines=[mappp.lines; jj zz];
% 
%     end
% end
% % map.lines=unique(mappp.lines, 'rows');
% map.lines=mappp.lines;

%% Nuclear facility
for i=1:size(mapp.lines,2)
    if mapp.points(3,mapp.lines(1,i))==0 && mapp.points(3,mapp.lines(2,i))==0
        x1=mapp.points(1,mapp.lines(1,i)); 
        y1=mapp.points(2,mapp.lines(1,i));
        x2=mapp.points(1,mapp.lines(2,i)); 
        y2=mapp.points(2,mapp.lines(2,i));
        for j=1:size(map.points,1)
            if x1*2000==map.points(j,1) && y1*2000==map.points(j,2) 
                jj=j;
            elseif x1*2000==map.points(j,2) && y1*2000==map.points(j,1)
                jj=j;
            end
        end
        for z=1:size(map.points,1)
            if x2*2000==map.points(z,1) && y2*2000==map.points(z,2) 
                zz=z;
            elseif x2*2000==map.points(z,2) && y2*2000==map.points(z,1)
                zz=z;
            end
        end
        mappp.lines=[mappp.lines; jj zz];

    end
end
map.lines=unique(mappp.lines, 'rows');



%% Plot the map
% for i=1:size(map.lines,1)
%     x1=map.points(map.lines(i,1),1);
%     y1=map.points(map.lines(i,1),2);
%     x2=map.points(map.lines(i,2),1);
%     y2=map.points(map.lines(i,2),2);
% 
%     x=[x1 x2];
%     y=[y1 y2];
%     hold on;
%     plot(x, y, '-r','LineWidth',2);
% end


max_x=max(map.points(:,1));
min_x=min(map.points(:,1));

max_y=max(map.points(:,2));
min_y=min(map.points(:,2));

%% RRT*
% Define start and goal points

startNode = [0, -750];
goalNode = [2500, 2000];


maxIter=6000;
maxDistance=0.06*2000; %0.02
goalThreshold=0.06*2000; %0.02
rewiringRadius=0.1*2000; %0.05

%Choose the number of feasible paths to find between startNode and goalNode
num_paths=1;

paths = inf(num_paths,3);
rewiring_times=[];

circle_radius =  0.003*2000; 
safety_margin = circle_radius + circle_radius*0.05; %+il 5%


for loopp = 1:1 
    i=0;
    while i<num_paths
        tree = Tree();
        tree.root = Node();
        tree.root.position = startNode;
        tree.root.parent = [];
        tree.root.cost = 0;
        tree.nodes = [];
        tree.nodes=[tree.nodes; tree.root];
    
        toc1_=0;
        toc2_=0;
        %% Main loop
        for iter = 1:maxIter
            
            % Attempt to extend the tree towards the random point
            [tree, nearNodes,const, newNode, rewiredNodes, exceedNodes, toc1,toc2] = extendTree(max_x, min_x, max_y, min_y, tree, maxDistance, rewiringRadius, map, safety_margin); %create q_new nodo (con coordinate)
            toc1_=toc1_+toc1;
            toc2_=toc2_+toc2;
            if exceedNodes==true
                break;
            end
            
            if exceedNodes==false %if the max number of tree.nodes has not been reached
                % Check if reached goal
                check=goalNode-newNode.position;
                if norm(check) < goalThreshold
                    goalNode_ = newNode;
                    tic;
                    path = findPath(tree, goalNode_);
                    toc3=toc;
                    toc1_=toc1_+toc3;
                    break;
                end
            end
            
        end
        
        if exceedNodes==false
            i=i+1;
            % visualizeTree(tree);
            % hold on;
            if ~isempty(path)
                % plot(path(:,1), path(:,2), 'k', 'LineWidth', 2);
                % color = rand(1, 3);
                % h=plot(path(:,1), path(:,2),'Color', color, 'LineWidth', 2);
                % 
                % frame = getframe(gcf);
                % writeVideo(video, frame); 
    
                % pause(0.5);
                % delete(h);
            end
            % disp('After concatenation:');
            % for i=1:size(tree.nodes,1)
            %    disp(tree.nodes(i));
            % end
        end
        %% Computation of the total path distance
        if exceedNodes==false
            x = path(:, 1);
            y = path(:, 2);
            dx = diff(x);
            dy = diff(y);
            distances = sqrt(dx.^2 + dy.^2);
            total_distance = sum(distances);
        
            % fprintf('Elapsed time is %.2f seconds.\n', elapsedTime);
            % total_time=time1+time2;
            paths((((loopp-1)*num_paths)+i),1)=loopp;
            paths((((loopp-1)*num_paths)+i),2)=toc1_-toc2_;
            paths((((loopp-1)*num_paths)+i),3)=total_distance;
            rewiring_times=[rewiring_times; toc2_];
        end
    
    
    end
end

%% Remove rows in which 2nd e 3rd elements are nulle since the max number of tree.nodes was eccedeed
rowsToRemove = isinf(paths(:, 2)) & isinf(paths(:, 3));
paths(rowsToRemove, :) = [];

hold on;

% startNode=[-400, 8];
% goalNode=[600, 950];
% 
% plot(startNode(1), startNode(2), 'gx', 'MarkerSize', 6, 'LineWidth', 2);
% plot(goalNode(1), goalNode(2), 'bx', 'MarkerSize', 6, 'LineWidth', 2);


% excel_file = 'paths_normal_RRT.xlsx';
% writematrix(paths, excel_file);


% close(video);





