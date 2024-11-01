function orientation = computeOrientation(p, q, r)
    v1 = q - p;
    v2 = r - q;
    orientation = det2D(v1, v2);
end