function HF = fill_holesDS(SV,SF)
  % FILL_HOLES Try to fill holes using a very simple technique, followed by a
  % "robust" back up. Neither is guaranteed to avoid things like
  % self-intersections, but the resulting mesh should be manifold at holes if
  % the input hole boundaries are manifold to begin with. Should probably try
  % to remove any duplicate vertices before calling. 
  %
  % Inputs:
  %   SV  #SV by 3 list of mesh vertex positions
  %   SF  #SF by 3 list of mesh triangle indices into SV
  % Outputs:
  %   HF  #HF by 3 list of indices into SV of triangles filling holes in
  %     (SV,SF).
  % 
  TR = triangulation(SF,SV);
  O = TR.freeBoundary();
%   O = outline(SF);
  CO = connected_components(O);
  CO = CO(O(:,1));
  [~,~,CO] = unique(CO);
  HF = [];
  for c = 1:max(CO)
    E = O(CO==c,:);
    [EV,IM] = remove_unreferenced(SV,E);
    J = [];
    J(IM) = 1:size(SV,1);
    E = IM(E);
    [~,A] = affine_fit(EV);
%     [EVV,EF] = triangle(EV*A,E,[],'Quiet');
    pts = EV*A;
    DT = delaunayTriangulation(pts(:,1),pts(:,2),E);
    
    [eu,cEu]=grIsEulerian(E);
    Hreord = E(cEu,:)';
    Hreord = Hreord(:);
    reordered = unique(Hreord,'stable');
    if eu==1
        reordered = [reordered; reordered(1)];
    end
    IC = incenter(DT);
    in = inpolygon(IC(:,1),IC(:,2),pts(reordered,1),pts(reordered,2));
    out = ~in;
    EVV = DT.Points;
    EF = DT.ConnectivityList;
    EF(out,:) = [];
    
    if size(EVV,1) ~= size(EV,1)
      EF = [];
      EE = E;
      while ~isempty(EE)
        [~,loop] = grIsEulerian(EE);
        %loop = full(outline_loopDS(EE));
        maxA = inf;
        maxEF = [];
        for s = 1:numel(loop)
          loop = loop([2:end 1]);
          loop_1 = loop([2:end 1]);
          loop_2 = loop_1([2:end 1]);
          sEF = [repmat(loop(1),numel(loop)-2,1) [loop_1(1:end-2) loop_2(1:end-2)]];
          A = sum(doublearea(EV,sEF));
          if A < maxA
            maxA = A;
            maxEF = sEF;
          end
        end
        EF = [EF;maxEF];
        EE = EE(~all(ismember(EE,loop),2),:);
      end
    end
    allE = [EF(:,2:3);EF(:,[3 1]);EF(:,1:2)];
    flip = ismember(E(1,:),allE,'rows');
    if flip
      EF = fliplr(EF);
    end
    HF = [HF;J(EF)];
  end

end
