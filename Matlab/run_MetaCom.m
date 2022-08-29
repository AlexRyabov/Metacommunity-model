
%% Run the metacommunity model
%% The function returns T(Save) - vector of time, N_Dyn (L, L, sp, TSpan) -- the dynamics
%% of species number. R_Dyn(L, L, res, TSpan) the dynamics of the resources

function [TSpan, N_Dyn, R_Dyn] = run_MetaCom...
    (N,  c,   R_st,    r,    m, d, ...
    R,     S,    D,    TSpan, dx, Lx, Ly, Params)
% init sp distr, sp cons rates,  R*,growth rate, mort,diffusiv,
% init res distr, suppl points,  Suppl rates,
% TimeSpan, flux, dx, L



%% Find an appropriate time step
Diff = max(d(:));
switch Params.Solver
    case 'Euler'  
        dtMax = 0.01;
    otherwise
        dtMax = 0.1*abs(TSpan(1)-TSpan(end));
end

if Diff == 0
    dt = dtMax;
else
    dt = min(dtMax, 0.5 * (dx^2/(4*Diff))); %time step should satisfy the stability criterion
end

Params.nSp = size(N, 3);  %number of species
Params.nRes = size(R, 3); %number of reosurces

% second edition

K = MC_GetHalfSaturation(R_st, r,  squeeze(m(1, 1, :)));

%% Initialyze the output aurgments
N_Dyn = zeros([size(N, 1), size(N, 2), size(N, 3), length(TSpan)]);
R_Dyn = zeros([size(R, 1), size(R, 2), size(R, 3), length(TSpan)]);


switch Params.Solver
    case 'Euler'  %explicit euler scheme. Very slow
        %% Save initial conditions
        CurrentTime = 0;
        N_Dyn(:, :, :, 1) = N;
        R_Dyn(:, :, :, 1) = R;
        t = 0;
        for ti=2:length(TSpan)
            SaveStep = 0;
            TimeStep = dt;
            while t < TSpan(ti)
                t = t + TimeStep;
                if t > TSpan(ti)
                    SaveStep = 1;
                    t = TSpan(ti);
                    TimeStep = dt - (t - TSpan(ti));
                end
                if (Params.FertTime > 0) && (t > Params.FertTime) && (Params.TreatmentPlot == 1) && (Params.Fertilized == 0)          % Fertilization
                    S(:, :, 1) = S(:, :, 1) + Params.SAddition(1);
                    S(:, :, 2) = S(:, :, 2) + Params.SAddition(2);
                    Params.Fertilized = 1;
                end
                
                
                %% Run the model
                [f_N, f_R] = SpatLocalChange(R, D, S, N, r, K, m, c);
                [N_rr] = SpatBoundary(Params.BoundaryCond, N); % boundary conditions preparation for N
                [N_Laplace] = Laplace2D_2x2(N_rr, dx); % Laplace operator
                for si=1:size(N, 3)
                    dN(:, :, si) = f_N(:, :, si) + d(si) * N_Laplace(2:end-1, 2:end-1, si);
                end
                
                %Use to have diffusion of the resources
                %         [R_rr] = SpatBoundary(Params.BoundaryCond, R); % optionally also for R
                %         [R_Laplace] = Laplace2D_2x2(R_rr,dx);
                %         for ri=1:res
                %             dR(:, :, ri) = f_R(:, :, ri) + DispRes(ri) * R_Laplace(2:end-1,2:end-1, ri);
                %         end
                dR = f_R; % differential equations for R1 and R2 (in this case both
                
                % Euler step for every iteration
                N = N + TimeStep*dN;
                R = R + TimeStep*dR;
                R(R<0)=0;
                CurrentTime = CurrentTime + TimeStep;
                if SaveStep
                    N_Dyn(:, :, :, ti) = N;
                    R_Dyn(:, :, :, ti) = R;
                    SaveStep = 0;
                    break;
                end
                %Reshuffle
%                 if (ReshuffleStep > 0) && (t > ReshuffleTime)
%                     %S = resourceshuffle(S);
%                     R = R +  5 * rand(size(R));
%                     ReshuffleTime = ReshuffleTime + rand() * ReshuffleStep;
%                 end
            end
        end
    case {'ode45', 'cvode'}
        N_Dyn(:, :, :, 1) = N;
        R_Dyn(:, :, :, 1) = R;
        
        rhsode =@(t, X) rhs_MetaCom(t, X, c, R_st,  r, K, m, d,  S,  D, dx, Lx, Ly, Params);
        switch Params.Solver
            case 'ode45'
                opts = odeset('MaxStep',dt);
                solver = @ode45;
            case 'cvode'
                % Options for integration
                opts = CVodeSetOptions('RelTol', 0, 'AbsTol', 1.e-5,'LinearSolver','Dense');
                solver = @CVodeSolver;
        end
        
        X0 = cat(3, N, R);
        X0 = X0(:);
        Tsol = [];
        Xsol = [];
        if Params.Env.ReshuffleTimes(1) > 0
            t0 = 0;
            for iRT = 1:length(Params.Env.ReshuffleTimes)
                %select TSpan from ReshTime0 to Params.Env.ReshuffleTimes(iRT)
                TSpanSub = TSpan(TSpan >= t0 & TSpan < Params.Env.ReshuffleTimes(iRT));
                if isempty(TSpanSub)
                    [TsolSub, XsolSub] = solver(rhsode, ([t0 Params.Env.ReshuffleTimes(iRT)]), X0, opts);
                else
                    [TsolSub, XsolSub] = solver(rhsode, ([TSpanSub Params.Env.ReshuffleTimes(iRT)]), X0, opts);
                end
                %Run model in the interval from ReshTime0   Params.Env.ReshuffleTimes(iRT)
                t0 = Params.Env.ReshuffleTimes(iRT);
                X0 = XsolSub(end, :);
                %save this solutions
                [~,~, ib]  = intersect(TSpanSub, TsolSub);
                Tsol = [Tsol; TsolSub(ib)];
                Xsol = [Xsol; XsolSub(ib, :)];
                %reshuffle resources
                S = Params.SSupplResourceShuffle(t0, S, Params);
                if sum(abs(Params.RAddition)) 
                     X0_shape = reshape(X0, Ly, Lx, Params.nSp + Params.nRes);
                      R_0(:, :, :) = X0_shape(:, :, Params.nSp + (1:Params.nRes));
                      for ir = 1:Params.nRes
                          R_0(:, :, ir) =  R_0(:, :, ir) + rand(size(R_0(:, :, ir))) * Params.RAddition(1); 
                      end
                      X0_shape(:, :, Params.nSp + (1:Params.nRes)) = R_0(:, :, :);
                      X0 = X0_shape(:)';
                end                
                rhsode =@(t, X) rhs_MetaCom(t, X, c, R_st,  r, K, m, d, S,  D, dx, Lx, Ly, Params);
            end
            %run the model for the rest of time
            if TSpan(end) > t0
                TSpanSub = TSpan(TSpan >= t0);
                if ~isempty(TSpanSub)
                    [TsolSub, XsolSub] = solver(rhsode, unique([t0, TSpanSub]), X0, opts);
                    [~,~, ib]  = intersect(TSpanSub, TsolSub);
                    %save this solutions
                    Tsol = [Tsol; TsolSub(ib)];
                    Xsol = [Xsol; XsolSub(ib, :)];
                end
            elseif TSpan(end) == t0
                Tsol = [Tsol; t0];  %%add the last solution obtained in the loop
                Xsol = [Xsol; X0];
            end
        elseif Params.PlotSolutions
            %open figures
            fg_res = figure();
            fg_spe = figure();
            fg_tr_biom = figure();
            %initial moment of time
            Tsol = TSpan(1);
            Xsol = X0';
            for iRT = 2:length(TSpan)
                [TsolSub, XsolSub] = solver(rhsode, TSpan(iRT-1:iRT), X0, opts);
                X0 = XsolSub(end, :);
                %save this solutions
                Tsol = [Tsol; TsolSub(end)];
                Xsol = [Xsol; XsolSub(end, :)];
                MetaCom_plotResults(Params, Tsol, Xsol, {fg_res, fg_spe, fg_tr_biom});
                rhsode =@(t, X) rhs_MetaCom(t, X, c, R_st,  r, K, m, d, S,  D, dx, Lx, Ly, Params);
            end
            
        else
            [Tsol, Xsol] = solver(rhsode, TSpan, X0, opts);
        end
        
        Xsol = reshape(Xsol, length(Tsol), Ly, Lx, Params.nSp + Params.nRes);
        Xsol = permute(Xsol, [2, 3, 4, 1]); %set time as the last dimension in the array;
        N_Dyn(:, :, :, :) = Xsol(:, :, 1:Params.nSp, :);
        R_Dyn(:, :, :, :) = Xsol(:, :, Params.nSp + (1:Params.nRes), :);
end

if Params.GetFinalSolultion
    X0 = cat(3, N_Dyn(:, :, :, end), R_Dyn(:, :, :, end));
    rhsfsolve = @(X) rhsode(Inf, X);
    Xfin = fsolve(rhsfsolve, X0); toc
    N_Dyn(:, :, :, end + 1) = Xfin(:, :, 1:Params.nSp);
    R_Dyn(:, :, :, end + 1) = Xsol(:, :, Params.nSp + (1:Params.nRes));
    TSpan(end + 1) = Inf;
end

end


function dXdt = rhs_MetaCom(t, X, c, R_st,  r, K, m, d,  S,  D, dx, Lx, Ly, Params)
%t time
%nSp number of species
%nRes number of resources
%X matrix Ly by Lx by nSp + nRes

X = reshape(X, Ly, Lx, Params.nSp + Params.nRes);
N = X(:, :, 1:Params.nSp);  %species
R = X(:, :, Params.nSp+(1:Params.nRes)); %resources


%% Run the model
[f_N, f_R] = MC_LocalChange(t, R, D, S, N, r, K, m, c, Params);
[N_rr] = SpatBoundary(Params.BoundaryCond, N); % boundary conditions preparation for N
[N_Laplace] = Laplace2D_2x2(N_rr, dx); % Laplace operator
for si=1:size(N, 3)
    dN(:, :, si) = f_N(:, :, si) + d(si) * N_Laplace(2:end-1, 2:end-1, si);
end

%Uncomment this part to have diffusion of the resources
if sum(Params.DispRes)>0
    [R_rr] = SpatBoundary(Params.BoundaryCond, R); % optionally also for R
    [R_Laplace] = Laplace2D_2x2(R_rr,dx);
    for ri=1:size(R, 3)
        dR(:, :, ri) = f_R(:, :, ri) + Params.DispRes * R_Laplace(2:end-1,2:end-1, ri);
    end
else
    dR = f_R;
end
dXdt = [dN(:); dR(:)];
end

function  [t,x] = CVodeSolver(f, TSpan, X0, opts)
CVodeInit(f, 'BDF', 'Newton', TSpan(1), X0, opts);
[status, t,x] = CVode(TSpan, 'Normal');
% Free solver memory
CVodeFree;
end
