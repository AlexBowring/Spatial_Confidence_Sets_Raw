function [C,h,Ph,F,Fa,Fc,k] = spm_reml_sc(YY,X,Q,N,hE,hC,A,K)
% ReML estimation of covariance components from y*y' - proper components
% FORMAT [C,h,Ph,F,Fa,Fc,k] = spm_reml_sc(YY,X,Q,N,[hE,hC,A,K]);
%
% YY  - (m x m) sample covariance matrix Y*Y'  {Y = (m x N) data matrix}
% X   - (m x p) design matrix
% Q   - {1 x q} covariance components
% N   - number of samples
%
% hE  - hyperprior expectation in log-space [default = -32]
% hC  - hyperprior covariance  in log-space [default = 256]
% A   - proportional hyperpriors [default = 1, yes]
% K   - number of iterations [default = 32]
%
% C   - (m x m) estimated errors = h(1)*Q{1} + h(2)*Q{2} + ...
% h   - (q x 1) ReML hyperparameters h
% Ph  - (q x q) conditional precision of log(h)
%
% F   - [-ve] free energy F = log evidence = p(Y|X,Q) = ReML objective
%
% Fa  - accuracy
% Fc  - complexity (F = Fa - Fc)
%
% k   - number of iterations required
%
% Performs a Fisher-Scoring ascent on F to find MAP variance parameter
% estimates.  NB: uses weakly informative log-normal hyperpriors.
% See also spm_reml for an unconstrained version that allows for negative
% hyperparameters
%
%__________________________________________________________________________
%
% SPM ReML routines:
%
%      spm_reml:    no positivity constraints on covariance parameters
%      spm_reml_sc: positivity constraints on covariance parameters
%      spm_sp_reml: for sparse patterns (c.f., ARD)
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_reml_sc.m 3264 2009-07-10 14:01:31Z karl $

% assume proportional hyperpriors not specified
%--------------------------------------------------------------------------
try, A; catch, A  = 1;  end
 
% assume a single sample if not specified
%--------------------------------------------------------------------------
try, N; catch, N  = 1;  end
 
% default number of iterations
%--------------------------------------------------------------------------
try, K; catch, K  = 64; end
 
% initialise h
%--------------------------------------------------------------------------
n     = length(Q{1});
m     = length(Q);
h     = zeros(m,1);
dh    = zeros(m,1);
dFdh  = zeros(m,1);
dFdhh = zeros(m,m);
 
% ortho-normalise X
%--------------------------------------------------------------------------
if isempty(X)
    X = sparse(n,0);
    R = speye(n,n);
else
    X = orth(full(X));
    R = speye(n,n) - X*X';
end
 
 
% initialise and specify hyperpriors
%==========================================================================

% scale Q and YY
%--------------------------------------------------------------------------
if A
    sY = trace(R*YY)/N/n;
    YY = YY/sY;
    for i = 1:m
        sh(i,1) = trace(R*Q{i})/n;
        Q{i}    = Q{i}/sh(i);
    end
else
    sY = 1;
    sh = 1;
end

% hyperpriors
%--------------------------------------------------------------------------
try, hE = hE(:);                               catch, hE = -32;   end
try, hP = inv(hC + speye(length(hC))/exp(16)); catch, hP = 1/256; end
 
% check sise
%--------------------------------------------------------------------------
if length(hE) < m, hE = hE(1)*ones(m,1);   end
if length(hP) < m, hP = hP(1)*speye(m,m);  end

% intialise h: so that sum(exp(h)) = 1
%--------------------------------------------------------------------------
if any(diag(hP) >  exp(16))
    h = hE;
end
 
% ReML (EM/VB)
%--------------------------------------------------------------------------
dF    = Inf;
as    = 1:m;
ds    = 1:m;
for k = 1:K
 
    % compute current estimate of covariance
    %----------------------------------------------------------------------
    C     = sparse(n,n);
    for i = as
        C = C + Q{i}*exp(h(i));
    end
    iC    = inv(C + speye(n,n)/exp(32));
 
    % E-step: conditional covariance cov(B|y) {Cq}
    %======================================================================
    iCX    = iC*X;
    if ~isempty(X)
        Cq = inv(X'*iCX);
    else
        Cq = sparse(0);
    end
 
    % M-step: ReML estimate of hyperparameters
    %======================================================================
 
    % Gradient dF/dh (first derivatives)
    %----------------------------------------------------------------------
    P     = iC - iCX*Cq*iCX';
    U     = speye(n) - P*YY/N;
    for i = as
 
        % dF/dh = -trace(dF/diC*iC*Q{i}*iC)
        %------------------------------------------------------------------
        PQ{i}   = P*Q{i};
        dFdh(i) = -trace(PQ{i}*U)*N/2;
 
    end
 
    % Expected curvature E{dF/dhh} (second derivatives)
    %----------------------------------------------------------------------
    for i = as
        for j = as
 
            % dF/dhh = -trace{P*Q{i}*P*Q{j}}
            %--------------------------------------------------------------
            dFdhh(i,j) = -trace(PQ{i}*PQ{j})*N/2;
            dFdhh(j,i) =  dFdhh(i,j);
 
        end
    end
 
    % modulate
    %----------------------------------------------------------------------
    dFdh  = dFdh.*exp(h);
    dFdhh = dFdhh.*(exp(h)*exp(h)');
 
    % add hyperpriors
    %----------------------------------------------------------------------
    e     = h     - hE;
    dFdh  = dFdh  - hP*e;
    dFdhh = dFdhh - hP;
 
    % Fisher scoring: update dh = -inv(ddF/dhh)*dF/dh
    %----------------------------------------------------------------------
    dh    = spm_dx(dFdhh(as,as),dFdh(as))*exp(-k/(K/2));
    h(as) = h(as) + dh;
    
    % Convergence (1% change in log-evidence)
    %======================================================================
    % bar(h);drawnow
    
    % convergence
    %----------------------------------------------------------------------
    dF    = dFdh(as)'*dh;
    fprintf('%-30s: %i %30s%e\n','  ReML Iteration',k,'...',full(dF));
    if dF < 1e-2
        break
    else
        % eliminate redundant components (automatic selection)
        %------------------------------------------------------------------
        as  = find(h > -16);
        as  = as(:)';
    end
end

% log evidence = ln p(y|X,Q) = ReML objective = F = trace(R'*iC*R*YY)/2 ...
%--------------------------------------------------------------------------
Ph    = -dFdhh;
if nargout > 3
 
    % tr(hP*inv(Ph)) - nh + tr(pP*inv(Pp)) - np (pP = 0)
    %----------------------------------------------------------------------
    Ft = trace(hP*inv(Ph)) - length(Ph) - length(Cq);
 
    % complexity - KL(Ph,hP)
    %----------------------------------------------------------------------
    Fc = Ft/2 + e'*hP*e/2 + spm_logdet(Ph*inv(hP))/2 - N*spm_logdet(Cq)/2;
 
    % Accuracy - ln p(Y|h)
    %----------------------------------------------------------------------
    Fa = Ft/2 - trace(C*P*YY*P)/2 - N*n*log(2*pi)/2  - N*spm_logdet(C)/2;
 
    % Free-energy
    %----------------------------------------------------------------------
    F  = Fa - Fc - N*n*log(sY)/2;
 
end


% return exp(h) hyperpriors and rescale
%--------------------------------------------------------------------------
h  = sY*exp(h)./sh;
C  = sY*C;
