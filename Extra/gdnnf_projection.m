function w = gdnnf_projection(y, v)
%gdnnf_projection  Auxiliary function for IR Tools
%
% w = gdnnf_projection(y, v)
%
% This function implements the projection needed in the projected 
% gradient methods.  Specifically, it solves the problem:
%
%     w = argmin||w - y||, such that w >= 0 and sum(w(:)) = v
%
% That is, we enforce nonnegativity on w, and preservation of volume of 
% the data.
%
% The solution to this problem is:
%
%     w = min(y - lambda, 0)
%
% where the Lagrange multiplier, lambda, satisfies:
%
%     sum_i max(y(i) - lambda, 0) = v,
%
% or lambda is a root of 
%
%     f(lambda) = sum_i max(y(i) - lambda, 0) - v = 0
%
% We find a good initial guess of the root, using the given y(i), and then
% use a couple Newton iterations to refine the approximation.A_times_vec.

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD Licence. A separate license file should be provided as part 
% of the package.

% Define the function f(lambda)
[M, N] = size(y);
y = y(:);
f = @(lambda) sum(max(y-lambda,0)) - v;

% Find a good initial guess of the root using the given y values.
% Specifically, f(lambda) is a monotonically decreasing function, so
% if y(1)<=y(2)<=...<=y(n), then we find the largest j such that 
%      sum_{j=k+1}^n (y(k) - y(j)) - v > 0
% and take the initial lambda to be y(j).
ys = sort(y,'ascend');
n = length(ys);
s = zeros(n,1);
s(1) = sum(ys(2:n)) - (n-1)*ys(1) - v;
for j = 1:n-1
  s(j+1) = s(j) - (n-j)*(ys(j+1)-ys(j));
end
idx = find(s>0, 1, 'last');
if isempty(idx)
  lambda = 0;
else
  lambda = ys(idx);
end
f0 = f(lambda);

% Now do a few Newton iterations.  Here the maximum number of iterations
% is set to 10, but generally only 2 or 3 iterations are needed.
for k = 1:10
  nmax = length(find(y-lambda>0));
  lambda = lambda + f0/nmax;
  f1 = f(lambda);
  if f0 - f1 < 1e-10
    break
  end
  f0 = f1;
end
w = max(y-lambda,0);
w = reshape(w, M, N);