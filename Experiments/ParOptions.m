function Options = ParOptions(varargin);
% function Options = ParOptions (varargin);
% 
% Description:
% Defines the Options for GenPARAFAC and GNPARAFAC
% 
% Inputs: 
% varargin: it can be one structure (in which case the function controls only if the values in the
% structure are available) or a list of option names and obtions values
% Options = ParOptions (optionname1,optionvalue1,optionname2,....);
% In the latter case the non specified options will have default value.
% 
% 
% Outputs:
% Options: Options structure, fields:
%          algorithm      : dgn    - damped Gauss-Newton (Levenberg - Marquadt method)                      (default)
%                           pmf3   - PMF3 (modified damped Gauss-Newton)
%                           als    - Alternating Least Squares
%                           swatld - Self-Weighted Alternating TriLinear Decomposition
%          comp_refmaxiter: Max number of iteration in the refining step                                    (default: 10000)
%          comp_tuckextra : Number of extra components in Tucker3 compression model.                        (default: 2)
%          comp_tuckiter  : Max number of iterations in Compression                                         (default: 3)
%          compress       : on/off(default) - activate/deactivate compression
%          cc_fit         : Convergence criterion -> value of loss function compared to Frobenius norm of X (default: 10 * eps)
%          cc_grad        : Convergence criterion -> infinite norm of the gradient                          (default: 1e-9)
%          cc_maxiter     : Convergence criterion -> max number of iterations                               (default: 2500)
%          cc_par         : Convergence criterion -> relative norm of the parameters update vector          (default: 1e-8)
%          cc_relfit      : Convergence criterion -> relative loss function decrease                        (default: 1e-6)
%          diagnostics    : on(default)/off - show/do not show diagnostics
%          display        : n  -> show fit every n iterations                                               (default: 5)
%                           -1 -> nothing is shown while iterating
%          init_method    : random - random values                                                          (default)
%                           orth   - orthogonalised random loadings
%                           svd    - svd
%                           dtld   - DTLD/GRAM
%                           swatld - SWATLD
%                           best   - Try all of the above and continue with the set giving the best fit
%          init_iter      : Max number of iterations for initialisation                                     (default: 5)
%          init_tol       : Tolerance for initialisation (when iterative methods are employed)              (default: 1e-5)
%          lambdainit     : Initial value of the damping parameter with respect to the largest value 
%                           of diag(J'*J)                                                                   (default: 1)
%          lambdaudpar    : Damping parameter update terms                                                  (default: [2,1/3])
%          reg_gamma      : Regularisation parameter (Gamma in the reference)                               (default: 1)
%          reg_gammaupdate: Update value for Gamma (i.e. Gamma is divided by reg_gammaupdate when updated)  (default: 5)
%          reg_iter       : Consecutive iterations with relative fit decrease lower than a threshold 
%                           defined based on cc_relfit (cf. reference for the PMF3 algorithm)               (default: 3)
%          reg_n_updates  : Number of updates for the regularisation term                                   (default: 5)
%
% Author: 
% Giorgio Tomasi 
% Royal Agricultural and Veterinary University 
% Rolighedsvej 30 
% DK-1958 Frederiksberg C 
% Denmark 
% 
% Last modified: 10-Jan-2005 15:24
% 
% Contact: Giorgio Tomasi, gt@kvl.dk; Rasmus Bro, rb@kvl.dk 
%

if nargin == 1 % Check options.
    Op = fieldnames(varargin{1});
    for i = 1:length(Op)
        
        V = varargin{1}.(Op{i});
        switch lower(Op{i})            
            case 'algorithm'
                if ~any(strcmp({'als','dgn','pmf3','swatld'},lower(V)))
                    error('Invalid assignment for ''Algorithm''. It must be: ''dGN'' (default), ''ALS'', ''pmf3'', ''SWATLD''')
                end
                
            case 'compress'   
                if ~any(strcmp({'on','off'},lower(V)))
                    error('Invalid assignment for ''Compress''. It can either be ''on'' or ''off''')
                end
                
            case 'compression'
                ComOp = fieldnames(V);
                for j = 1:length(ComOp)
                    
                    U = V.(ComOp{j});
                    switch lower(ComOp{j})
                        case 'tuckeriter'
                            if ~isa(U,'double') | length(U) ~= 1
                                error('Invalid property assignment for ''Comp_TuckerIter''. Value must be a scalar.')
                            elseif rem(U,1)
                                error('Invalid property assignment for ''Comp_TuckerIter''. Value must be an integer.')
                            end         
                            
                        case 'tuckerextra'    
                            if ~isa(U,'double') | length(U) ~= 1
                                error('Invalid property assignment for ''Comp_TuckExtra''. Value must be a scalar.')
                            elseif rem(U,1)
                                error('Invalid property assignment for ''Comp_TuckExtra''. Value must be an integer.')
                            end
                            
                        case 'refmaxiter'
                            if ~isa(U,'double') | length(U) ~= 1
                                error('Invalid property assignment for ''Comp_RefMaxIter''. Value must be a scalar')
                            elseif rem(U,1)
                                error('Invalid property assignment for ''Comp_RefMaxIter''. Value must be an integer')
                            end
                            
                        otherwise
                            error('Invalid property name')
                            
                    end
                end
                
            case 'convcrit'
                CCOp = fieldnames(V);
                for j = 1:length(CCOp)
                    
                    U = V.(CCOp{j});
                    switch lower(CCOp{j})                     
                        case 'fit'
                            if ~isa(U,'double') | length(U) ~= 1
                                error('Invalid property assignment for ''CC_Fit''. Value must be a scalar.')
                            end
                            
                        case 'grad'
                            if ~isa(U,'double') | length(U) ~= 1
                                error('Invalid property assignment for ''CC_Grad''. Value must be a scalar.')
                            end
                            
                        case 'maxiter'
                            if ~isa(U,'double') | length(U) ~= 1
                                error('Invalid property assignment for ''CC_MaxIter''. Value must be an scalar.')
                            elseif rem(U,1)
                                error('Invalid property assignment for ''CC_MaxIter''. Value must be an integer.')
                            end         
                            
                        case 'par'
                            if ~isa(U,'double') | length(U) ~= 1
                                error('Invalid property assignment for ''CC_Par''. Value must be a scalar.')
                            end
                            
                        case 'relfit'
                            if ~isa(U,'double') | length(U) ~= 1
                                error('Invalid property assignment for ''CC_RelFit''. Value must be a scalar.')
                            end
                            
                        otherwise
                            error('Invalid property name')
                            
                    end
                    
                end
                
            case 'diagnostics'
                if ~any(strcmp({'on','off'},lower(V)))
                    error('Invalid assignment for ''Diagnostics''. It can either be ''on'' or ''off''.')
                end
                
            case 'display'
                if ~(isa(V,'char') | isa(V,'double'))
                    error('Invalid property assignment for ''Display''. Value must be an integer or ''None''.')
                end    
                if (isa(V,'char') & ~strcmp('none',lower(V))) | (isa(V,'double') & (logical(rem(V,1)) | length(V) ~= 1))
                    error('Invalid property assignment for ''Display''. Value must be an integer or ''None''.')
                end
                
            case 'initialisation'
                InOp = fieldnames(V);
                for j = 1:length(InOp)
                    
                    U = V.(InOp{j});
                    switch lower(InOp{j})
                        case 'method'
                            if ~any(strcmp({'random','orth','svd','dtld','swatld','best'},lower(U)))
                        
                                error('Invalid property assignment for ''Init_Method''. Value must be either ''random'',''orth'',''svd'',''dtld'' or ''swatld''')
                            end
                            
                        case 'maxiter'
                            if ~isa(U,'double') | length(U) ~= 1
                                error('Invalid property assignment for ''Init_MaxIter''. Value must be a scalar.')
                            elseif rem(U,1)
                                error('Invalid property assignment for ''Init_MaxIter''. Value must be an integer.')
                            end         
                            
                        case 'tol'
                            if ~isa(U,'double')
                                error('Invalid property assignment for ''Init_Tol''. Value must be a scalar')
                            end
                            
                        otherwise
                            error('Invalid property name')
                            
                    end
                    
                end
                
            case 'lambdainit'
                if ~isa(V,'double') & length(V) ~= 1
                    error('Invalid property assignment for ''LambdaInit''. Value must be a scalar.')
                end
                
            case 'lambdaudpar'
                if ~isa(V,'double') | length(V) ~= 2
                    error('Invalid property assignment for ''LambdaUDPar''. Value must be a vector of 2 scalars.')
                end
                if V(1) <= 1
                    error('Invalid property assignment for ''LabdaUDPar'', Value(1) must be strictly greater than 1.')
                end    
                if V(2) >= 1
                    error('Invalid property assignment for ''LabdaUDPar'', Value(2) must be strictly less than 1.')
                end         
                
            case 'regularisation'
                RegOp = fieldnames(V);
                for j = 1:length(RegOp)
                    
                    U = V.(RegOp{j});
                    switch lower(RegOp{j})
                        case 'gamma'
                            if ~isa(U,'double') | length(U) ~= 1
                                error('Invalid assignment for ''Reg_Gamma''. It must be a scalar.')
                            end
                            
                        case 'gammaupdate'
                            if ~isa(U,'double') | length(U) ~= 1
                                error('Invalid assignment for ''Reg_GammaUpdate''. It must be a scalar.')
                            end
                            
                        case 'iter'
                            if ~isa(U,'double') | length(U) ~= 1
                                error('Invalid property assignment for ''Reg_Iter''. Value must be a scalar.')
                            elseif rem(U,1)
                                error('Invalid property assignment for ''Reg_Iter''. Value must be an integer.')
                            end         
                            
                        case 'n_updates'
                            if ~isa(U,'double') | length(U) ~= 1 
                                error('Invalid property assignment for ''Reg_Update''. Value must be an scalar.')
                            elseif rem(U,1)
                                error('Invalid property assignment for ''Reg_Update''. Value must be an integer.')
                            end         
                            
                        otherwise
                            error('Invalid property name')
                            
                    end
                    
                end
                
            otherwise
                error('Invalid property name')
                
        end
    end
    Options = varargin{1};
    
else
    
    CC      = struct('fit',10 * eps,'grad',1e-9,'maxiter',10000,'par',1e-8,'relfit',1e-6);   % Convergence criteria
    Comp    = struct('refmaxiter',10000,'tuckerextra',2,'tuckeriter',3);                    % Compression options
    Init    = struct('maxiter',5,'method','random','tol',1e-5);                             % Initialisation options
    Reg     = struct('gamma',1,'gammaupdate',5,'iter',3,'n_updates',5);                     % Regularisation options for pmf3
    Options = struct('algorithm','dgn','compress','on','compression',Comp,...               % Default values for Options
        'convcrit',CC,'diagnostics','off','display',1,'initialisation',Init,...
        'lambdainit',1,'lambdaudpar',[2,1/3],'regularisation',Reg);
    
    if rem(nargin,2)
        
        if isa(varargin{1},'struct')
            if ~all(strcmp(sort(fieldnames(varargin{1})),sort(fieldnames(Options))))
                error('Incompatible structures')
            else 
                Options_old = varargin{1};
                Options     = varargin{1};
                varargin(1) = [];
                EM          = 'The input Options structure is returned';
            end
        else
            error('Properties and values must come in pairs')
        end
        
    else
        Options_old = Options;
        EM          = 'The default Options structure is returned';
    end
    
    for i = 1:2:length(varargin)
        
        if isequal(lower(varargin{i}(1:min(5,length(varargin{i})))),'init_') 
            Options.initialisation.(lower(varargin{i}(6:end))) = lower(varargin{i+1});
        elseif isequal(lower(varargin{i}(1:min(5,length(varargin{i})))),'comp_')
            Options.compression.(lower(varargin{i}(6:end))) = lower(varargin{i+1});
        elseif isequal(lower(varargin{i}(1:min(3,length(varargin{i})))),'cc_')
            Options.convcrit.(lower(varargin{i}(4:end))) = lower(varargin{i+1});
        elseif isequal(lower(varargin{i}(1:min(4,length(varargin{i})))),'reg_')
            Options.regularisation.(lower(varargin{i}(5:end))) = lower(varargin{i+1});
        else
            Options.(lower(varargin{i})) = lower(varargin{i+1});
        end
        
    end
    
    try
        Options = paroptions(Options);
    catch
        warning('')
        fprintf(['The required update generated an error:\n',lasterr,'\n'])
        fprintf(EM)
        Options = paroptions(Options_old);
    end
    
end
