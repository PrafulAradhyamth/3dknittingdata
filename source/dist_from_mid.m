function dmid = dist_from_mid(TR_start,mid)
[tm,bd,myParam] = setup_distance_computation(TR_start,mid);
SEmex  = SolveEikonal2Dmex(tm,bd,myParam,[]);
% SEmex  = SolveEikonal2D(tm,bd,myParam);
dmid  = SEmex.Compute_Soln;
