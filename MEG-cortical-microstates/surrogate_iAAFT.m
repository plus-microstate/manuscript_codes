% --------------- iAAFT Surrogates ----------------------------------------
function x_surrogate = surrogate_iAAFT(x,maxiter,mvflag)

    % check inputs
    [m,n] = size(x) ; 
    
    if nargin < 2
        maxiter = 10 ; 
    end
    
    if nargin < 3
        mvflag = false ; 
    end
 
    % ranked distribution
    [sx,idx0] = sort(x) ; 
    % fourier amplitudes
    Fx = fft(x) ; 
    Af = abs(Fx) ; 
    Ph = angle(Fx) ; % only required for multivariate
   
    % randomly shuffle data
    sx0 = zeros(size(x)) ; 
    for i = 1:n
        idx = randperm(m) ; 
        sx0(:,i) = x(idx,i) ; 
    end

    % iterate algorithm
    % sx0: surrogate data at start of iteration
    % sx1: surrogate data at end of iteration
    % Sxn: fourier transform of surrogate data (n=0,1)
    for iter = 1:maxiter
        
        % fourier transform of shuffled data
        Sx0 = fft(sx0) ; 
        % phase stays the same, but amplitudes replaced with amplitudes of 
        % signal
        phi0 = angle(Sx0) ;
        
        if mvflag
            % Get angle of rotation for multivariate correlations
            dphi = phi0-Ph ;
            cosdphi = cos(dphi) ; 
            sindphi = sin(dphi) ; 
            tana = sum(sindphi,2)./sum(cosdphi,2) ; 
            al = atan(tana) ; 
            
            % Ensure correct quadrant
            al = al+pi*(sum(cos(al-dphi),2)<0) ;
            
            % Update angle
            phi0 = Ph+al ; 
            
        end
            
            
        Sx1 = Af.*exp(1i.*phi0) ; 
        % inverse transform 
        sx1 = real(ifft(Sx1)) ; 
        % replace the distribution of sx1 with the values from eeg
        [~,idx1] = sort(sx1) ; [~,idx1] = sort(idx1) ; 
        for i = 1:n
            sx1(:,i) = sx(idx1(:,i),i) ; 
        end 
        % Check for convergence
        if all(all(sx1==sx0))
            break
        else
            idx0 = idx1 ; 
        end
        % update iteration
        sx0 = sx1 ; 
    end
    x_surrogate = sx1 ; 
    
end