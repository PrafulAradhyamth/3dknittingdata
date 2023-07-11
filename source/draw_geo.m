function lh = draw_geo(fig,bndflag,G,V,srci,desti,alg)
if bndflag
    disp('Boundary mode. Press S to change to seam mode.')
    P = shortestpath(G,srci,desti);
    P = P(:);
    lh = line(V(P,1),V(P,2),V(P,3),'color','r','linewidth',3,'UserData',[P,ones(length(P),1)]);
else
    disp('Seam mode. Press B to change to boundary mode.')
    src = {geodesic_create_surface_point('vertex',srci,V(srci,:))};
    dest = geodesic_create_surface_point('vertex',desti,V(desti,:));
    set(fig,'Pointer', 'watch');
    drawnow;
    geodesic_propagate(alg, src, {dest}, 2);
    set(fig,'Pointer', 'arrow');
    path = geodesic_trace_back(alg, dest);
    path = flipud(path);
    [x,y,z,P,type] = extract_coordinates_from_path(path);
    lh = line(x,y,z,'color','r','linewidth',3,'UserData',[P,type]);
end