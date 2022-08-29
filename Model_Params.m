function [dx, Lx, Ly, tmin, tmax,  ...
    res, D, S, R, ResDyn, d, r, m, R_eq, c, N, ModParams] = ...
Model_Params (ModParams, iParam)


Diff = ModParams.Diff;
tmax = ModParams.tmax;
sp   = ModParams.sp;

%% Model parameters
dx = 1; % space resolution

if Diff == 0
    dt = 0.1;
else
    dt = min(0.1, 0.5 * (dx^2/(4*Diff))); %time step should satisfy the stability criterion
end
tmin = 0; % starting time
t = tmin:dt:tmax; % time vector


%rand('seed',1); % optional repetition of the same random values
Lx = ModParams.Lx;
Ly = ModParams.Ly;
%% Resources
res =  size(ModParams.RStar_All, 2); %Number of resources
D = 0.25 * ones(Ly, Lx, res); %  supply rate of the resources (in 1/d) in each point

% % calculate the resource distribution
% if strcmp(ModParams.ReSup, 'eS')   %%equidistantly distributed
%         %generate resource supply points for the equal required range
%         S(:, :, 1) = repmat(linspace(ModParams.ResLow, ModParams.ResHigh, Lx), Ly, 1); % values of S1 and S2 are equidistantly distributed
%         S(:, :, 2) = ModParams.R2toR1SupplRatio * repmat(linspace(ModParams.ResHigh, ModParams.ResLow, Lx), Ly, 1); % values of S1 and S2 are randomly generated
% elseif strcmp(ModParams.ReSup, 'rS')   %%uniformly randomly distributed
%         % generate resource supply point for the unequal required range
%         S_1_pre = [ModParams.ResLow + (ModParams.ResHigh-ModParams.ResLow)*rand(Lx-2,1)];
%         S_1 = sort([ModParams.ResLow; ModParams.ResHigh; S_1_pre]);
%         S_2 = ModParams.R2toR1SupplRatio * (ModParams.ResHigh - S_1 + ModParams.ResLow);
%         %flipud(sort([ModParams.ResLow; ModParams.ResHigh;(ModParams.ResHigh-S_x_pre)+1]));
%         S(:, :, 1) = repmat(S_1', Ly, 1); % values of S1 and S2 are randomly generated
%         S(:, :, 2) = repmat(S_2', Ly, 1); % values of S1 and S2 are randomly generated
% end
% 
% switch ModParams.ResDistr
%     case 'l'
%         %use this for scenario: Random Locations
%         S = resourceshuffle(S, res);
%     case 's'
%         %use this for scenario: Random Supplies
%         S = randresshuffle(S, res);
%     case 'g'
%         %use this for scenario: gradient Supplies
%     otherwise 
%         error('The resource type is undefined');
% end


% Setup resource supply points, if they are not defined
if isempty(ModParams.S)
    switch ModParams.ReSup
        case 'eS'   %generate resource supply points for the equal required range
            S(:, :, 1) = repmat(linspace(ModParams.ResLow, ModParams.ResHigh, Lx), Ly, 1); % values of S1 and S2 are equidistantly distributed
            S(:, :, 2) = ModParams.R2toR1SupplRatio * (ModParams.ResHigh -  S(:, :, 1) + ModParams.ResLow); %
        case 'rSh'   %random generated and shifted from  the diagonal ResLow_i<R_i<ResHigh_i
            %requires definition of  ResLow_i<R_i<ResHigh_i  for both resources
            for ir = 1:res
                S(:, :, ir)  = ModParams.ResLow(ir) + (ModParams.ResHigh(ir) -ModParams.ResLow(ir)) * rand(Ly, Lx);
            end
        case 'rS'   % generate resource supply point for the unequal required range
            S1_pre = [ModParams.ResLow + (ModParams.ResHigh-ModParams.ResLow)*rand(Lx-2,1)];
            S1 = sort([ModParams.ResLow; ModParams.ResHigh; S1_pre]);
            S2 = ModParams.R2toR1SupplRatio * (ModParams.ResHigh - S1 + ModParams.ResLow);
            S(:, :, 1) = repmat(S1', Ly, 1); % values of S1 and S2 are randomly generated
            S(:, :, 2) = repmat(S2', Ly, 1); % values of S1 and S2 are randomly generated
        case {'rLogNorm', 'rLogNorm_g'}   % generate a lognormal distribution of resources with mean = (ResLow + ResHigh)/2 and var = (ResHigh- ResLow)^2;
            %rLogNorm_g sort from left to right
            for ir = 1:res
                if length(ModParams.ResHigh) == res
                    m2 = ((ModParams.ResHigh(ir) + ModParams.ResLow(ir))/2).^2; 
                    v2 = (ModParams.ResHigh(ir) - ModParams.ResLow(ir)).^2;
                elseif length(ModParams.ResHigh) == 1
                    m2 = ((ModParams.ResHigh + ModParams.ResLow)/2).^2; 
                    v2 = (ModParams.ResHigh - ModParams.ResLow).^2;
                else
                    error('length(ModParams.ResHigh) does not equal res');
                end
                mu = log(m2./sqrt(v2 + m2));                       %%see matlab help for lognrnd
                sigma = sqrt(log(v2./m2 + 1));                     %%see matlab help for lognrnd
                S(:, :, ir) = lognrnd(mu, sigma, Ly, Lx); % values of S1 and S2 are randomly generated
            end
            switch ModParams.ReSup

                    %%
            end
        case 'ranSupOff' % generate randomly points along a resource gradient
            if iParam == 1
                slopeFunc = (ModParams.ResHigh-ModParams.ResLow)/(ModParams.ResLow-ModParams.ResHigh);
                slopeIntersect = ModParams.ResHigh-(slopeFunc*ModParams.ResLow);
                randS1_xP_Dummy = sort((ModParams.ResHigh-ModParams.ResLow).*rand(Lx,1)+ModParams.ResLow);
                randS2_yP_Dummy = (randS1_xP_Dummy.*slopeFunc)+slopeIntersect;
                S(:, :, 1) = repmat(randS1_xP_Dummy', Ly, 1); % values of S1 and S2 are randomly generated
                S(:, :, 2) = ModParams.R2toR1SupplRatio * repmat(randS2_yP_Dummy', Ly, 1); % values of S1 and S2 are randomly generated
            else
                S(:, :, 1) = ModParams.S(:, :, 1);
                S(:, :, 2) = ModParams.S(:, :, 2);
            end
        case 'convexS'
            xS1 = linspace(ModParams.ResLow, ModParams.ResHigh+2, Lx/2);
            Sfunc = @(x) (-0.0006*x^3+0.0638*x^2-2.6255*x+43.0737);
            yS1 = arrayfun(Sfunc,xS1);
            S_pre = sort([xS1,yS1]);
            yS1_pre = arrayfun(Sfunc,S_pre);
            S(:, :, 1) = repmat(S_pre, Ly, 1);
            S(:, :, 2) = repmat(yS1_pre, Ly, 1);
    end
    switch ModParams.ResDistr
        case 'l'
            %use this for scenario: Random Locations
            S = resourceshuffle(S, res);
        case 's'
            %use this for scenario: Random Supplies
            S = randresshuffle(S, res);
        case {'g', 'n'}
            %use this for scenario: gradient Supplies
            Ratio = S(:, :, 1)./S(:, :, 2);
            [~, ind] = sort(Ratio(:));
            %%
            for ir = 1:res
                Stmp = S(:, :, ir);
                Stmp = Stmp(ind);
                Stmp = vec2mat(Stmp,size(S, 2));
                S(:, :, ir) = Stmp;
            end
        otherwise
            error('The resource type is undefined');
    end
    ModParams.S = S;
else
    S = ModParams.S;
end

%Initial resource values
R = S; % initial values of R1 and R2 equal the S1 and S2 values 
ResDyn = []; %The dynamics of the average resource values

%% Species parameters
%sp = 20; %Species number
d = Diff * ones(sp, 1); % diffusivity for each species

%maximal growth rate
if isempty(ModParams.MaxGrowthRate) 
    ModParams.MaxGrowthRate = 1 * ones(sp, 1); %maximal growth rate
end
r = ModParams.MaxGrowthRate; 


%ParamDefault.m is empty, then setup mortality
if isempty(ModParams.m)
    m = 0.25* ones(Ly, Lx, sp); %mortality of each species in each point
end

if isempty(ModParams.RStar)
    if strcmp(ModParams.TraitDist, 'eT')  %%equidistantly distributed along a trade-off line
            % critical resource requirements R_eq(species#, resource#)
            ModParams.RStar = [linspace(ModParams.TraitLow, ModParams.TraitHigh, sp)' linspace(ModParams.TraitHigh, ModParams.TraitLow, sp)'];  %Trade-off in resource requirements       
            R_eq = ModParams.RStar;
    elseif strcmp(ModParams.TraitDist, 'rT')  %%uniformly randomly distributed along a trade-off line
            % critical resource requirements unequally distributed
            R_x = sort([ModParams.TraitLow;ModParams.TraitHigh;ModParams.TraitLow+(ModParams.TraitHigh-ModParams.TraitLow)*rand(sp-2,1)]);
            R_y = flipud(R_x);
            ModParams.RStar = [R_x,R_y];
            R_eq = ModParams.RStar;
    end
else
   R_eq = ModParams.RStar;
end


%consumtion of resource j by rate species i. c(species#, resource#)
%optimal foraging theory --> origin to R_q values defines slope of consumption vectors
switch ModParams.comptype
    case 'C'
        alpha = 0.05;   % consumption rate
        ModParams.consumprate = alpha.*R_eq;
        c = ModParams.consumprate;
        %         ModParams.consumprate = [linspace(ModParams.TraitLow/ModParams.consump, ModParams.TraitHigh/ModParams.consump, sp)' ...
        %             linspace(ModParams.TraitHigh/ModParams.consump, ModParams.TraitLow/ModParams.consump, sp)'];
        %         c = ModParams.consumprate
    case 'B'
        %bistability %works only for 2 reosurces
        alpha = 0.05;   % consumption rate
        ModParams.consumprate = alpha.*([R_eq(:,2),R_eq(:,1)]);
        c = ModParams.consumprate;
end

%Initial species number for randomly chosen set of species
%N = 0.01 * rand(L,L,sp);
N = 0.1 * mean(S(:)) * rand(Ly,Lx,sp);
All_spec = 1:sp;

