function [A,B,C,Diag,Time] = Parafac3w(varargin);
% function [A,B,C,Diag,Time] = PARAFAC3W (X,F,Options,Ai,Bi,Ci);
% 
% Description:
% Fits a PARAFAC model to a 3-way array.
% 
% Inputs: 
% X                  : 3-way array
% F                  : rank of the model
% Options (optional) : see PAROptions
% Ai,Bi,Ci (optional): cell vector with the initial estimations for the factors (each element contains a
%                      matrix with the loadings in the corresponding mode).
% 
% 
% Outputs:
% A,B,C: loading matrices
% Diag : structure with some diagnostic values
%           Diag.fit(1)     : value of the loss function after initialisation
%           Diag.fit(2)     : value of the loss function after fitting in the compressed space
%           Diag.fit(3)     : value of the loss function after fitting (refining) in the original space
% 
%           Diag.it(1)      : # iterations in fitting
%           Diag.it(2)      : # iterations in refining the model (if compression is used)
%
%           Diag.convergence: met convergence criteria (only for LM algorithm). See GNPARAFAC_LM and GNPARAFAC_PMF3
% Time : time spent in the different stages of the fitting algorithm
%           Time(1) : compression (if applied)
%           Time(2) : initialisation (if required)
%           Time(3) : fitting of the model
%           Time(4) : refining of the model (if compression is applied)
%           Time(5) : total time employed
% 
% 
% 
% Subroutines:
% Internal: svdf, tucker3, vec (returns vectorised argument)
% External: check, genparafac, gnparafac, initpar, kr, nshape, parafac, paroptions, scale_factors
% 
% 
% Author: 
% Giorgio Tomasi 
% Royal Agricultural and Veterinary University 
% MLI, LMT, Chemometrics group 
% Rolighedsvej 30 
% DK-1958 Frederiksberg C 
% Denmark 
% 
% Last modified: 10-Jan-2005 15:21
% 
% Contact: Giorgio Tomasi, gt@kvl.dk, Rasmus Bro, rb@kvl.dk 
%

% Check for minimal input
if ~nargin
    help parafac3w
    return
end

[X,F,Options,Start] = Check_GenParafac_Input(varargin{:});
if isempty(F)
   [A,B,C,Diag,Time] = deal([]);
   return
end

tot_t               = cputime; %Initialise timer for time consumption

% Some initial values
dimX   = size(X);                                             % Array's size vector
SSX    = tss(X,false);                                        % Array's total sum of squares
Time   = zeros(4,1);                                          % Initialise Outputs
Diag   = struct('fit',zeros(3,1),'it',zeros(4,1),'convergence',zeros(7,1));

% Dimensions of the subspaces used for compression
W = min(Options.compression.tuckerextra + F,dimX);

% Display some information before fitting
if isequal(Options.diagnostics,'on')
    ShowDiagnostics(Options,F,SSX,dimX,~isempty(Start))
end

% Display information while iterating
ShowIter = Options.display;
if isequal(Options.display,'none')
    ShowIter = NaN;
end

% Compress the array by fitting a Tucker3 model (uses Tucker.m function from the n-way toolbox)
if strcmp(Options.compress,'on')
    Time0            = cputime;
    TuckOpt          = [1e-5 0 0 0 NaN Options.compression.tuckeriter];                             % Option for fitting the Tucker3 model
    [CompBas,Gt,fit] = tucker(X,W,TuckOpt,[],[],{rand(dimX(1),F),rand(dimX(2),F),rand(dimX(3),F)}); % The Tucker algorithm is initialised with loading matrices of random numbers
    Time(1)          = cputime - Time0;
else
    Gt = X;                                         % Nothing is done to X
end

% Initialize PARAFAC
if isempty(Start)
    
    Time0                 = cputime;
    [CompIni{1:ndims(X)}] = InitPar(Gt,F,Options);  % Compute initial estimates for the PARAFAC loading matrices
    if isequal(Options.compress,'on')               
        for m = 1:ndims(X)
            Start{m} = CompBas{m} * CompIni{m};     % Expand initial values to the original array size
        end
    else
        Start = CompIni;                            % No expansion is needed if compression is not applied
    end
    Time(2,1) = cputime - Time0;                    % Time spent for initialisation
    
else
    
    if isequal(Options.compress,'on')
        for m = 1:ndims(X)
            CompIni{m} = CompBas{m}' * Start{m};    % Compress initial estimates if given as a input parameters
        end
    else
        CompIni = Start;                            % If compression is not applied, start values are used 'as is'
    end
    
end

% Fit rank F PARAFAC model (to the compressed array if compression is used)
Time0  = cputime;
switch Options.algorithm
    case 'als'      % Use Parafac.m from the n-way toolbox
        OptionsALS            = [Options.convcrit.relfit 0 0 0 ShowIter Options.convcrit.maxiter];
        [FacComp,Diag.it(1)]  = parafac(Gt,F,OptionsALS,[],CompIni,[]);
        if Diag.it(1) == Options.convcrit.maxiter
            Diag.convergence(6) = true;
        else
            Diag.convergence(1) = true;
        end
        
    case 'pmf3'     % Use the Positive Matrix Factorization for 3-way arrays algorithm (cf. reference)
        [FacComp,Diag] = PMF3(Gt,[],Options,CompIni{:});
        
    case 'dgn'      % Use a damped Gauss-Newton algorithm (cf. reference)
        [FacComp,Diag] = dGN(Gt,[],Options,CompIni{:});
        
    case 'swatld'   % Uses the Self-Weighted Alternating Trilinear Decomposition algorithm (cf. reference)
        [FacComp,Diag] = SWATLD(Gt,[],Options,CompIni{:});

    otherwise
        error('Reference to non existing algorithm')
        
end
Time(3,1)   = cputime - Time0;                              % Time for fitting the model
Diag.fit(1) = tss(X - nmodel(Start),false);                 % Fit using initial estimates

% Expand solution if compression is used
if isequal(Options.compress,'on')
    for m = 1:ndims(X)
        FacInit{m} = CompBas{m} * FacComp{m};
    end
else
    FacInit = FacComp;   
end

Diag.fit(2) = tss(X - nmodel(FacInit),false);               % Compute loss function value after fitting the model (of the expanded
                                                            % solution if compression is employed)

% Refine solution using PARAFAC-ALS if compression was applied
if isequal(Options.compress,'on')
    
    OptionsALS                      = [Options.convcrit.relfit 0 0 0 ShowIter Options.compression.refmaxiter];
    Time0                           = cputime;
    [FacFin,Diag.it(4),Diag.fit(3)] = parafac(X,F,OptionsALS,[],FacInit);
    Time(4,1)                       = cputime - Time0;
    if Diag.it(4) ~= Options.convcrit.maxiter
        Diag.convergence(7) = true;
    end
    
else
    FacFin = FacInit;
end

%Scale factors according to standard convention
[FacFin{end:-1:1}] = scale_factors(1,FacFin{end:-1:1}); 

%Sort the factors according to their norm (in decresing order)
[nil,Seq] = sort(-sum(FacFin{1}.^2));
[A,B,C]   = facperm(Seq,FacFin{:});

if strcmp(Options.diagnostics,'on')
    
    Precision = sprintf(' %%1.%ie',-log10(Options.convcrit.relfit));
    if isequal(lower(Options.compress),'on')
    
        if any(strcmpi({'dgn','pmf3'},Options.algorithm))
            fprintf('\n\n # Iterations\n Fitting   : %i (%i Hessian computations)\n Refining  : %i',Diag.it([1,2,4]))
        else
            fprintf('\n\n # Iterations\n Fitting   : %i\n Refining  : %i',Diag.it([1,4]))
        end
        fprintf(['\n\n Loss function value\n Init.     :',Precision],Diag.fit(1))
        fprintf(['\n Fitting   :',Precision],Diag.fit(2))
        fprintf(['\n Refining  :',Precision],Diag.fit(3))
        Precision = sprintf(' %%3.%if',-log10(Options.convcrit.relfit) - 2);
        fprintf(['\n\n Explained variance \n Fitting  :',Precision,'%%\n Refining :',Precision,'%%\n'],100*(1 - Diag.fit(2:3)/SSX))
        
    else
        
        if any(strcmpi({'dgn','pmf3'},Options.algorithm))
            fprintf('\n # Iterations\n Fitting   : %i (%i Hessian computations)\n Refining  : -',Diag.it([1,2]))
        else
            fprintf('\n # Iterations\n Fitting   : %i\n Refining  : -',Diag.it(1))
        end
        fprintf(['\n\n Loss function value\n Init.     :',Precision],Diag.fit(1))
        fprintf(['\n Fitting   :',Precision],Diag.fit(2))
        fprintf('\n Refining  : -')
        Precision = sprintf(' %%3.%if',-log10(Options.convcrit.relfit) - 2);
        fprintf(['\n\n Explained variance \n Fitting  :',Precision,'%%\n Refining : -'],100*(1 - Diag.fit(2)/SSX))
        
    end
    fprintf('\n\n Time consumption\n Compress  : %g s\n Init      : %g s\n Fitting   : %g s\n Refining  : %g s\n Total     : %g s',Time,cputime - tot_t)
    
end

