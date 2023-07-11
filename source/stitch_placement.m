function jac = stitch_placement(jac)
% stitch placement optimization according to liu et al
% trline = [];
% for i=1:size(jac,1)-1
%     w1 = find(jac(i,:));
%     w2 = find(jac(i+1,:));
%     if isempty(intersect(w1,w2))
%         if max(w2)<min(w1)
%             jac(i,min(w2):min(w1)-1) = true;
%         else
%             jac(i,max(w1)+1:max(w2)) = true;
%         end
%         trline = cat(1,trline,i);
%     end
% end



% jacjac = [jac;jac];
% jacjac(1:2:end,:) = jac;
% jacjac(2:2:end,:) = jac;
% jac = jacjac;
% jac(trline*2,:) = [];
% 
% % odd to even
% for i=1:2:size(jac,1)-1
%     j0 = find(jac(i,:),1,'first');
%     j1 = find(jac(i+1,:),1,'first'); 
%     if j0~=j1
%          vec = min(j0,j1):max(j0,j1)-1;
%          jac(i+1,vec) = jac(i,vec);                  
%     end
%     j0 = find(jac(i,:),1,'last'); % !!!!
%     j1 = find(jac(i+1,:),1,'last'); 
%     if abs(j0-j1)>1
%         vec = [j0-1,j0+1];
%         [~,idx] = min(abs(vec-j1));
%         vec2 = min(j1+1,vec(idx)):max(j1,vec(idx));
%         jac(i+1,vec2) = jac(i,vec2);
%     end
% end
% 
% % even to odd
% for i=2:2:size(jac,1)-1
%     j0 = find(jac(i,:),1,'last');
%     j1 = find(jac(i+1,:),1,'last'); 
%     if j0~=j1
%          vec = min(j0,j1)+1:max(j0,j1);
%          jac(i+1,vec) = jac(i,vec);                  
%     end
%     j0 = find(jac(i,:),1,'first');
%     j1 = find(jac(i+1,:),1,'first'); 
%     if abs(j0-j1)>1
%         vec = [j0-1,j0+1];
%         [~,idx] = min(abs(vec-j1));
%         vec2 = min(j1,vec(idx)):max(j1-1,vec(idx));
%         jac(i+1,vec2) = jac(i,vec2);
%     end
% end

% no common wales
for i=1:size(jac,1)-1
    w1 = find(jac(i,:));
    w2 = find(jac(i+1,:));
    if isempty(intersect(w1,w2))
        if max(w2)<min(w1)
            jac(i,max(w2):min(w1)-1) = 1;
        else
            jac(i,max(w1)+1:min(w2)) = 1;
        end
    end
end

% double courses
jacjac = [jac;jac];
jacjac(1:2:end,:) = jac;
jacjac(2:2:end,:) = jac;
jac = jacjac;
% correct two decreases
% for i=1:size(jac,1)-1
%     j0 = find(jac(i,:),1,'last');
%     j1 = find(jac(i+1,:),1,'last'); 
%     if j0-j1==2
%         jac(i+1,j1+1) = true;
%     end
%     j0 = find(jac(i,:),1,'first');
%     j1 = find(jac(i+1,:),1,'first'); 
%     if j1-j0==2
%         jac(i+1,j1-1) = true;
%     end
% end 

figure;
[r,c] = find(jac==1);
scatter(c,r,'b.'); axis equal; hold on
[r,c] = find(jac==2);
scatter(c,r,'r.'); 
[r,c] = find(jac==3);
scatter(c,r,'m.');
grid on

% adjust even courses
for i=2:2:size(jac,1)-2
    j0 = find(jac(i,:),1,'last');
    j1 = find(jac(i+1,:),1,'last');
%     if ~ismember(j0-j1,[0,1])
    if j0~=j1%abs(j0-j1)>1
        if j1<j0
            jac(i,j1+1:j0) = 0;            
        else
            jac(i,j0+1:j1) = 1;
        end
        if jac(i+1,j1) == 3
            jac(i,j1) = 3;
        end
%          vec = min(j0,j1)+1:max(j0,j1);
%          jac(i,vec) = 2-jac(i,vec);
%          jac(i+2,vec) = ~jac(i+2,vec);         
    end   
end
% adjust odd courses
for i=1:2:size(jac,1)-3
    j0 = find(jac(i,:),1,'first');
    j1 = find(jac(i+1,:),1,'first');
%     if ~ismember(j0-j1,[0,1])
    if j0~=j1%abs(j0-j1)>1
        if j0<j1
            jac(i,j0:j1-1) = 0;        
        else
            jac(i,j1:j0-1) = 1;
        end
        if jac(i+1,j1) == 2
            jac(i,j1) = 2;
        end

%     if j0~=j1
%         if j0<j1
%             jac(i+1,j0:j1-1) = 2;        
%         else
%             jac(i+1,j1:j0-1) = 0;
%         end
%         if jac(i,j0) == 1
%             jac(i+1,j0) = 1;
%         end  
    
%          vec = min(j0,j1):max(j0,j1)-1;
%          jac(i,vec) = 2-jac(i,vec);
%          jac(i+2,vec) = ~jac(i+2,vec);         
    end   
end



figure;
[r,c] = find(jac==1);
scatter(c,r,'b.'); axis equal; hold on
[r,c] = find(jac==2);
scatter(c,r,'r.'); 
[r,c] = find(jac==3);
scatter(c,r,'m.');
grid on

