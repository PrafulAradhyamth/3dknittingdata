function [tri_idx_new, level_set_new, level_field_new] = findContourPathKP(tri_idx,level_set, level_field, TR)
    function [visititedTriangles, visititedTrianglesPairsIndex, pathValid]=findContourPathNextStep(tri_idx, TR, pointsStartingTriangle, visititedTriangles, visititedTrianglesPairsIndex, noMinusPoints, hasLeftStartingPosition)        
        pathValid = false;
        lastVisitedTriangle = visititedTriangles(end);
        posLastVisitedTriangles = find(tri_idx == lastVisitedTriangle);
        posLastVisitedTrianglesPairsIndex = posLastVisitedTriangles - (posLastVisitedTriangles > trianglesPairsCount) * trianglesPairsCount;
        
        posCandidatesTrianglesPairsIndex=setdiff(posLastVisitedTrianglesPairsIndex,visititedTrianglesPairsIndex);
        for posCandidateTrianglesPairsIndex=posCandidatesTrianglesPairsIndex'
            visititedTrianglesTemp = visititedTriangles;
            value1 = tri_idx(posCandidateTrianglesPairsIndex,1);
            value2 = tri_idx(posCandidateTrianglesPairsIndex,2);
            
            if lastVisitedTriangle == value1;                         
                    nextTriangleToVisit = value2;
                else
                    nextTriangleToVisit = value1;
            end
            
			%Handling Minus Ends
            if ~noMinusPoints & nextTriangleToVisit < 0
                pointsCurrentMinus = TR.ConnectivityList(tri_idx(posCandidateTrianglesPairsIndex),:);
                
                if isempty(intersect(pointsStartingTriangle, pointsCurrentMinus))                
                    visititedTriangles(end+1) = nextTriangleToVisit;
                    visititedTrianglesPairsIndex(end+1) = posCandidateTrianglesPairsIndex;
                    pathValid=true;
                    %Ende gefunden
                end
			%Handling Ends which has no Minus - Ends if finds itself but only if has once left starting position
            elseif visititedTriangles(1) == nextTriangleToVisit
				if hasLeftStartingPosition 
					visititedTriangles(end+1) = nextTriangleToVisit;
                    visititedTrianglesPairsIndex(end+1) = posCandidateTrianglesPairsIndex;
                    pathValid=true;
                    %Ende gefunden					
				end
            else
				% Checking if starting position was left
				if noMinusPoints & ~hasLeftStartingPosition & isempty(intersect(TR.ConnectivityList(visititedTriangles(1),:), TR.ConnectivityList(nextTriangleToVisit,:)))
					hasLeftStartingPosition = true;
				end
                if not(any(visititedTriangles(:) == nextTriangleToVisit))
                    visititedTrianglesTemp = visititedTriangles;
                    visititedTrianglesTemp(end+1) = nextTriangleToVisit;
                    
                    visititedTrianglesPairsIndexTemp = visititedTrianglesPairsIndex;
                    visititedTrianglesPairsIndexTemp(end+1) = posCandidateTrianglesPairsIndex;
                    
                    [visititedTrianglesCheck, visititedTrianglesPairsIndexCheck, pathValid]=findContourPathNextStep(tri_idx, TR, pointsStartingTriangle, visititedTrianglesTemp, visititedTrianglesPairsIndexTemp, noMinusPoints, hasLeftStartingPosition);
                    
                    if pathValid
                        visititedTrianglesPairsIndex = visititedTrianglesPairsIndexCheck;
                        visititedTriangles = visititedTrianglesCheck;                        
                        break
                    end
                end                 
            end
        end
    end
    
    noMinusPoints = false;
    visititedTrianglesMinus = [];
    
    trianglesPairsCount = size(tri_idx,1);    
    startingTriangleIndex = find(tri_idx < 0)-trianglesPairsCount;
    
    % If there is no MinusOne Point then we take the first point as
    if isempty(startingTriangleIndex)
        startingTriangleIndex = [startingTriangleIndex;1]
        noMinusPoints = true;
    end
        
    tri_idx_new = [];
    level_set_new = [];
    level_field_new = cell(1,length(level_field));
    
    pathValid = false;

    for startingTriangleIndexCurrent = startingTriangleIndex'        
        if isempty(intersect(visititedTrianglesMinus, tri_idx(startingTriangleIndexCurrent,2)))
            pointsStartingTriangle = TR.ConnectivityList(tri_idx(startingTriangleIndexCurrent),:);
            visititedTriangles = [];	
            % If there Minus Points we add the Minus Triangle
            if noMinusPoints
                visititedTrianglesPairsIndex = [];                
            else
                visititedTrianglesPairsIndex = startingTriangleIndexCurrent;    
                visititedTriangles(end+1) = tri_idx(startingTriangleIndexCurrent,2);
            end
            
            visititedTriangles(end+1) = tri_idx(startingTriangleIndexCurrent);
            
            [visititedTriangles, visititedTrianglesPairsIndex, pathValid]=findContourPathNextStep(tri_idx, TR, pointsStartingTriangle, visititedTriangles, visititedTrianglesPairsIndex, noMinusPoints, false);
            if pathValid & isempty(intersect(tri_idx_new, tri_idx(visititedTrianglesPairsIndex)))
                tri_idx_new = cat(1, tri_idx_new, tri_idx(visititedTrianglesPairsIndex,:));
                level_set_new =  cat(2, level_set_new, level_set(:,visititedTrianglesPairsIndex));
                for j=1:length(level_field)
                    level_field_new{j} = cat(1, level_field_new{j}, level_field{j}(visititedTrianglesPairsIndex));
                end
            end
            visititedTrianglesMinus = cat(2, visititedTrianglesMinus, visititedTriangles(find(visititedTriangles < 0)));
        end
        if pathValid;
            %Use only first found path
            break;
        end
    end
end   