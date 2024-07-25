function distance = determineDistance(x, y)

distance = sqrt((x-x').^2 + (y-y').^2);
distance = distance(tril(true(size(distance)),-1));