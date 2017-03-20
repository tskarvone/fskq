function N = fss_numel(u, d)
% FSS_NUMEL - the number of elements in a fully symmetric set
%   Given a generator vector u, N = FSS_NUMEL(u) computes
%   the number of elements in the fully symmetric set generated
%   by u.
%
%   Given two arguments u and d, N = FSS_NUMEL(u, d) computes
%   the number of elements  of the symmetric set in dimension d 
%   generated by u appended with d-length(u) zeros.
%
% INPUT
%   - u   the generator vector
%   - d   dimension (optional), if not set, assumed equal to length(u)      
%
% OUTPUT
%   - U   the fully symmetric set generated by u

% Toni Karvonen, 2017

  % Some argument checks
  if nargin < 2
    d = length(u);
  end
  
  if length(u) > d
    error('Length of the generator vector u cannot exceed the given dimension d.');
  end
  
  % Compute the number of non-zero elements
  uu = abs(u(:));
  uu = sort(uu(uu > 1e-14));
  dnz = length(uu);
  
  if dnz == 0
    N = 1;
    return
  end
  
  % eln contains multiplicities of element in
  % the full generator (that is, also multiplicity
  % of zeroes)
  [uuq, iuu] = unique(uu);
  eln = [diff(iuu); dnz - iuu(end) + 1; d - dnz];
  
  % Compute the number of elements
  N = 2^dnz * factorial(d) / prod(factorial(eln));
  
end