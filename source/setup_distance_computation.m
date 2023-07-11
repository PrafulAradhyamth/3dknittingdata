function [tm,bd,myParam] = setup_distance_computation(TR,boundary,bd_values)
tm.Vtx = TR.Points;
tm.DoFmap = TR.ConnectivityList;
tm.NegMask = false(length(TR.Points),1);

bd.Nodes = boundary; 
if nargin == 2
    bd.Data = zeros(length(bd.Nodes),1);
else
    bd.Data = bd_values;
end

myParam.Max_Tri_per_Star = 1200;
myParam.NumGaussSeidel = 150;
myParam.INF_VAL = 100000000;
myParam.TOL = 1e-12;

