function [map, triangulation]=mapSelection(mapSelectionInput)
    mappp.points=[];
    mappp.lines=[];
    
    if mapSelectionInput==1
        mapp = import_data_catia_vrml('warehouse.wrl');
        for i=1:size(mapp.points,2)
            if mapp.points(1,i)==0.02
                mappp.points=[mappp.points; mapp.points(2,i) mapp.points(3,i)]; % to put in vertical
            end
        end
        
        map.points=unique(mappp.points,'rows'); % to select the unique points
        map.points=map.points*2000;
        
        for i=1:size(mapp.lines,2)
            if mapp.points(1,mapp.lines(1,i))==0.02 && mapp.points(1,mapp.lines(2,i))==0.02
                x1=mapp.points(2,mapp.lines(1,i)); 
                y1=mapp.points(3,mapp.lines(1,i));
                x2=mapp.points(2,mapp.lines(2,i)); 
                y2=mapp.points(3,mapp.lines(2,i));
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
        %% Triangulation
        triangulation = delaunayTriangulation(map.points);

    elseif mapSelectionInput==2
        mapp = import_data_catia_vrml('warehouse2.wrl');
        for i=1:size(mapp.points,2)
            if mapp.points(1,i)==-0.02
                mappp.points=[mappp.points; mapp.points(2,i) mapp.points(3,i)]; % to put in vertical
            end
        end

        map.points=unique(mappp.points,'rows'); % to select the unique points
        map.points=map.points*2000;

        for i=1:size(mapp.lines,2)
            if mapp.points(1,mapp.lines(1,i))==-0.02 && mapp.points(1,mapp.lines(2,i))==-0.02
                x1=mapp.points(2,mapp.lines(1,i));
                y1=mapp.points(3,mapp.lines(1,i));
                x2=mapp.points(2,mapp.lines(2,i));
                y2=mapp.points(3,mapp.lines(2,i));
                for j=1:size(map.points,1)
                    if x1*2000==map.points(j,1) && y1*2000==map.points(j,2)
                        jj=j;
                    end
                end
                for z=1:size(map.points,1)
                    if x2*2000==map.points(z,1) && y2*2000==map.points(z,2)
                        zz=z;
                    end
                end
                mappp.lines=[mappp.lines; jj zz];
    
            end
        end
        triangulation = delaunayTriangulation(map.points, [9 10; 10 70; 70 69; 69 10; 69 9; 7 8; 8 68; 68 67; 67 7; 8 67; 4 30; 12 29;  11 63; 48 66; 47 64; 64 63; 66 65; 3 65]);
    elseif mapSelectionInput==3
        mapp = import_data_catia_vrml('nuclearFacility.wrl');
        for i=1:size(mapp.points,2)
            if mapp.points(3,i)==0
                mappp.points=[mappp.points; mapp.points(1,i) mapp.points(2,i)]; % to put in vertical
            end
        end

        map.points=unique(mappp.points,'rows'); % to select the unique points
        map.points=map.points*2000;

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
        %% Triangulation
        triangulation = delaunayTriangulation(map.points, [77 78; 79 143; 110 162;109 110;189 161;112 164; 163 111; 114 166; 165 113; 109 161; 193 194; 194 188; 188 193; 193 199; 161 187; 188 161; 188 187]);
    end
  
    map.lines=mappp.lines;
end