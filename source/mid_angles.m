function angles = mid_angles(mididx, points)
vecs = points(mididx(2:end),:) - points(mididx(1:end-1),:);
angles = acos(( dot(vecs(1:end-1,:),vecs(2:end,:),2) ) ./ (vecnorm(vecs(1:end-1,:),2,2) .* vecnorm(vecs(2:end,:),2,2)));