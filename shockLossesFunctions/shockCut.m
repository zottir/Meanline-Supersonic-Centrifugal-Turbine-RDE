function [x_shock,y_shock,phi] = shockCut(x_shock,y_shock,phi,x_i,y_i)
% Function to cut the shock close to the intersection point in a stable
% way, indipendently on the shape of the shock

% Parametrized Curve to stabilize the cut (Cumulative Length)
s = [0, cumsum(sqrt(diff(x_shock).^2 + diff(y_shock).^2))];

% Remove duplicates in s 
[s_unique, idx_unique] = unique(s, 'stable');
x_shock = x_shock(idx_unique);
y_shock = y_shock(idx_unique);
phi = phi(idx_unique);
s = s_unique;

% Find the closest point
[~, idx_cut_up] = min((x_shock - x_i).^2 + (y_shock - y_i).^2);

% Densify shock curve until intersection point
s_dense_up = linspace(s(1), s(idx_cut_up), 500);

x_shock = interp1(s(1:idx_cut_up), x_shock(1:idx_cut_up), s_dense_up, 'linear');
y_shock = interp1(s(1:idx_cut_up), y_shock(1:idx_cut_up), s_dense_up, 'linear');
phi = interp1(s(1:idx_cut_up), phi(1:idx_cut_up), s_dense_up, 'linear');

end