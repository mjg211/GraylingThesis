function [eq3_13, eq3_12, eq3_14, eq3_15] = groupSeqCrossoverMatrices(D)

% Check input value of D is an integer
if (rem(D, 1) ~= 0)
  error('D must be an integer greater than or equal to 2.')
end

% Initialise required symbolic variables
e = sym('e');
b = sym('b');
l = sym('l');
n = sym('n');

% Determine different forms for equations (3.12) and (3.13) depending on the
% number of drugs remaining in the trial. Return each in a single matrix
eq3_13 = b + zeros((D - 1)*D, D);
eq3_12 = b + zeros((D - 1)*D, D);
for r = D:-1:2
  Sigma_r    = b^2 + e^2*eye(r);
  invSigma_r = inv(Sigma_r);
  Sigma_r    = [Sigma_r; zeros(D - r, r)];
  invSigma_r = [invSigma_r; zeros(D - r, r)];
  Sigma_r    = [Sigma_r, zeros(D, D - r)];
  invSigma_r = [invSigma_r, zeros(D, D - r)];
  eq3_13((1 + (D - r)*D):((D + 1 - r)*D), :) = invSigma_r;
  eq3_12((1 + (D - r)*D):((D + 1 - r)*D), :) = Sigma_r;
end

% Determine different forms for equation (3.14) depending on the number of
% drugs remaining in the trial. Return in a single matrix
eq3_14 = b + zeros((2*D - 1)*(D - 1), 2*D - 1);
for r = D:-1:2
  sumXTinvSigmaX = b + zeros(2*D - 1, 2*D - 1);
  for p = 1:(2*D - 1)
    for q = 1:(2*D - 1)
      if ((p == 1) && (q == 1))
        sumXTinvSigmaX(p, q) = r;
      elseif (((p == 1) && (1 < q)) || ((1 < p) && (q == 1)))
        sumXTinvSigmaX(p, q) = 1;
      elseif ((p == q) && (1 < p))
        sumXTinvSigmaX(p, q) = 1 + (r - 1)*b^2/e^2;
      elseif ((p ~= q) && (p > 1) && (p < D + 1) && (q > 1) && (q < D + 1))
        sumXTinvSigmaX(p, q) = -b^2/e^2;
      elseif ((p ~= q) && (p > D) && (q > D))
        sumXTinvSigmaX(p, q) = -b^2/e^2;
      else
        sumXTinvSigmaX(p, q) = 1/r;
      end
      if (((p > r) && (p < D + 1)) || ((q > r) && (q < D + 1)) || (p > D + r - 1) || (q > D + r - 1))
        sumXTinvSigmaX(p, q) = 0;
      end
    end
  end
  sumXTinvSigmaX = n*e^2*sumXTinvSigmaX/(e^2*(e^2 + r*b^2));
  eq3_14((1 + (2*D - 1)*(D - r)):((2*D - 1)*(D + 1 - r)), :) = sumXTinvSigmaX;
end

% Determine the form for equation (3.15) depending on the value of D
eq3_15 = b + zeros(2*D - 1, 2*D - 1);
for p = 1:(2*D - 1)
  for q = 1:(2*D - 1)
    if ((p == 1) && (q == 1))
      eq3_15(p, q) = b^2 + ((2*D - 1)/D)*e^2;
    elseif (((p == 1) && (q > 1)) || ((p > 1) && (q == 1)))
      eq3_15(p, q) = -e^2;
    elseif ((p == q) && (p > 1))
      eq3_15(p, q) = 2*e^2;
    elseif ((p ~= q) && (p > 1) && (p < D + 1) && (q > 1) && (q < D + 1))
      eq3_15(p, q) = e^2;
    elseif ((p ~= q) && (p > D) && (q > D))
      eq3_15(p, q) = e^2;
    else
      eq3_15(p, q) = 0;
    end 
  end
end
eq3_15 = eq3_15/(l*n);

end