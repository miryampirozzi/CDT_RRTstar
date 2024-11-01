function intersect = segmentsIntersect(p1, q1, p2, q2)
    orient1 = computeOrientation(p1, q1, p2);

    orient2 = computeOrientation(p1, q1, q2);
  
    orient3 = computeOrientation(p2, q2, p1);

    orient4 = computeOrientation(p2, q2, q1);

    % Check the intersection conditions
    if (orient1 * orient2 < 0) && (orient3 * orient4 < 0)
        intersect = true;
    elseif (orient1 == 0) && isOnSegment(p1, q1, p2)
        % p2 in on p1-q1
        intersect = true;
    elseif (orient2 == 0) && isOnSegment(p1, q1, q2)
        % q2 in on p1-q1
        intersect = true;
    elseif (orient3 == 0) && isOnSegment(p2, q2, p1)
        % p1 in on p2-q2
        intersect = true;
    elseif (orient4 == 0) && isOnSegment(p2, q2, q1)
        % q1 in on p2-q2
        intersect = true;
    else
        intersect = false;
    end
end

function onSegment = isOnSegment(p, q, r)

    if (min(p(1), q(1)) <= r(1) && r(1) <= max(p(1), q(1))) && ...
       (min(p(2), q(2)) <= r(2) && r(2) <= max(p(2), q(2)))
        onSegment = true;
    else
        onSegment = false;
    end
end