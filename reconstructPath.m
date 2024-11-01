function path = reconstructPath(triangulation, visited_order)
    edges=[];
    for i=1:size(visited_order,1)
        edges=[edges; triangulation.ConnectivityList(visited_order(i),1) triangulation.ConnectivityList(visited_order(i),2); triangulation.ConnectivityList(visited_order(i),2) triangulation.ConnectivityList(visited_order(i),3); triangulation.ConnectivityList(visited_order(i),3) triangulation.ConnectivityList(visited_order(i),1)];
    end
    hold on;
    f=plot(triangulation.Points(edges',1),triangulation.Points(edges',2),'-g','LineWidth',2);

    % frame = getframe(gcf);
    % writeVideo(video, frame); 
    % 
    % pause(1);
    % delete(f);
end