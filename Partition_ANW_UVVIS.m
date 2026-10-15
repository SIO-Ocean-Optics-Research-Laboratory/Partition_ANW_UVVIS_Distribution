function [Lambda_out, adgopt, aphopt, N_feasible] = Partition_ANW_UVVIS(Lambda_anw,anw,PT)
%Description:                                                                   
%   Partitions a spectrum of the non-water absorption coefficient, anw,
%   into two constituents: the non-phytoplankton absorption coefficient,adg, 
%   and the phytoplankton absorption coefficient, aph. The model
%   operates in the near ultraviolet (UV) and visible (VIS) spectral regions
%   covering the light wavelength range from 350 to 700 nm.
%
%References:
%Stramski, D., and R. A. Reynolds. 2026. Model for partitioning the light
%	absorption coefficient of seawater within the near-UV and visible 
%	spectral region into phytoplankton and non-phytoplankton components.
%	Applied Optics, 65(30), 10467-10482. https://doi.org/10.1364/AO.610292
%	Abbreviated as [SR26] in the comments used to document the code
%
%Inputs:
%   Lambda_anw [n-by-1 numeric]: Values of light wavelength [nm] for the input 
%       spectrum of anw. Must be in the range of 350 to 700 nm. (REQUIRED)
%
%   anw [n-by-1 numeric]: Input values of the non-water absorption
%       coefficient anw [m^-1] corresponding to Lambda_anw. As a minimum
%       requirement, values at 7 specific wavelengths (within 5 nm) of 412, 
%       443, 470, 490, 510, 555, and 670 nm are required for the model to 
%       operate. (REQUIRED)
%       Input anw at additional wavelengths (e.g., hyperspectral) are used
%       by the model if available. The wavelengths of input anw are
%       critically important because the model derives aph only at these
%       wavelengths.
%
%   PT [1-by-m numeric]: Additional percentile values from the pool of feasible
%       solutions to return in addition to the optimal solutions. If not
%       specified, the default is to return the 50th and 99th percentile
%       values, i.e. PT = [50 99]. This default range of percentile
%       solutions is expected to have a high probability (generally >88%
%	    in the VIS, see [SR26]) of containing the true solution. (OPTIONAL)
%
%Additional inputs:
%   shapelib: A stored library of characteristic spectral shapes for adg 
%       provided in a separate file called 'ADG_shapelibrary_UVVIS_498.mat',
%       which must exist in the user path. This file contains a structure
%       named shapelib with the following two variables; 
%       shapelib.lambda (351-by-1 array of wavelengths from 350 to 700 nm at 
%       1-nm interval), and
%       shapelib.adg (351-by-N array of normalized spectral shapes, where N
%       represents the number of adg shapes in the library.
%       The current version of library has 498 spectral shapes. (REQUIRED)
%
%   Constraints [8-by-2 numeric]: An array specifying the lower (column 1) and 
%       upper (column 2) bounds for each of the 8 inequality constraints used
%       by the model. Default values from [SR26] are used but can be modified
%       by the user in the function code below.
%
%Outputs: 
%   Lambda_out [351-by-1 double]: Light wavelengths [nm] corresponding to
%       the model output values of the spectral absorption coefficient adg.
%       For the default adg shape library, this will always be 350 to 700 nm at
%       1 nm intervals. The wavelengths for model output values of aph will
%       correspond to the input values of Lambda_anw.
%
%   adgopt [351-by-(1+length(PT) double): Array containing the optimal
%       solutions for spectral values of adg [m^-1] determined as the 85th
%       percentile value of all feasible solutions at each wavelength
%       (column 1), and additional requested percentile solutions from the
%       pool of feasible solutions (remaining output columns in order of
%       the specified input PT).
%
%   aphopt [n-by-(1+length(PT)) double]: Array containing the optimal
%       solutions for spectral values of aph [m^-1] calculated as the
%       difference between input anw and adgopt (column 1) and additional
%       desired percentile solutions of adg from the pool of feasible
%       solutions (remaining output columns in order of PT). n corresponds
%       to the number of light wavelengths for input Lambda_anw.
%
%   N_feasible [1-by-1 double]: The number of feasible solutions found by the
%       model.
% 
%Version: 1.1
%
%Contacts: R. A. Reynolds, D. Stramski
%
%Version History:
%2024-12-02: Matlab version v1.0 released      
%2026-09-15: Matlab version v1.1 released
	%Added modification to remove from pool of feasible solutions all
		%occurrences where modeled adg exceeds input anw before
		%determination of percentiles and optimal solutions
	%Added citations and equation references to Stramski & Reynolds 2026 [SR26]
                                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%====CHECK INPUT ARGUMENTS====
    arguments
        Lambda_anw (:,1) {mustBeNumeric}
        anw (:,1) {mustBeNumeric}
        %if not specified, default to 50th and 99th percentile values for
        %output in addition to optimal solution
        PT (1,:) {mustBeNumeric} = [50 99]
    end
    nPT = size(PT,2); %number of percentile values to include in output

    %Check for existence of adg shape library file in path
    if ~isfile('ADG_shapelibrary_UVVIS_498.mat')
        error('Error. File ADG_shapelibrary_UVVIS_498.mat not found in path.');
    end

%=====DEFINE CONSTRAINTS=====
%Default array of 8 inequality constraints used by the model, [LB UB]
%where LB and UB refer to Lower and Upper Bound of values, respectively
%See [SR26] Tables 1 & 2
    Constraints = [  0.65	1.15;...  %Constraint 1: aph412/aph443 
                     0.55	0.81;...  %Constraint 2: aph490/aph443
                     0.67	1.41;...  %Constraint 3: aph470/aph412
                     0.07	0.52;...  %Constraint 4: aph555/aph490
                     1.60	6.71;...  %Constraint 5: aph510/aph555
                     1.20	4.60;...  %Constraint 6: aph443/aph670
                     0.33	0.76;...  %Constraint 7: (adg412*anw443) / anw412^2
                     0.004	0.010];   %Constraint 8: (aphn510-anph555) / (555-510)

%=====LOAD LIBRARY OF adg SPECTRAL SHAPES=====
%See [SR26] section 3B & Fig. 7    
    load('ADG_shapelibrary_UVVIS_498.mat','shapelib'); %must exist in path
    nshapes = size(shapelib.adg,2); %number of shape functions for adg_hat

%=====GET WAVELENGTH INFORMATION FOR EVALUATION OF CONSTRAINTS=====
% Get indices for the required wavelengths used to evaluate constraints
%
% FOR INPUT anw (SAMPLE): If an exact wavelength-based match is missing,
% the nearest wavelength within 5 nm is chosen
    [dist412,idx412S] = min(abs(Lambda_anw-412)); if dist412>5, idx412S = []; end
    [dist443,idx443S] = min(abs(Lambda_anw-443)); if dist443>5, idx443S = []; end
    [dist470,idx470S] = min(abs(Lambda_anw-470)); if dist470>5, idx470S = []; end
    [dist490,idx490S] = min(abs(Lambda_anw-490)); if dist490>5, idx490S = []; end
    [dist510,idx510S] = min(abs(Lambda_anw-510)); if dist510>5, idx510S = []; end
    [dist555,idx555S] = min(abs(Lambda_anw-555)); if dist555>5, idx555S = []; end
    [dist670,idx670S] = min(abs(Lambda_anw-670)); if dist670>5, idx670S = []; end
    clear dist*;

% FOR THE ADG SHAPE LIBRARY
    idx412L = find(shapelib.lambda == 412);
    idx443L = find(shapelib.lambda == 443);
    idx470L = find(shapelib.lambda == 470);
    idx490L = find(shapelib.lambda == 490);
    idx510L = find(shapelib.lambda == 510);
    idx555L = find(shapelib.lambda == 555);
    idx670L = find(shapelib.lambda == 670);

%=====GENERATE POOL OF SPECULATIVE SOLUTIONS=====                               
%Basic equations: [SR26] Eqs. 1 & 2
% anw = aph_prime + adg_prime
% adg_prime = A * adg_hat ; adg_prime is a linear function of the shape factors

%Create range of valid values to use for Constraints 1 (Cxi) and 2 (Cyj)
%[SR26] Eqs. 7a & 7b
%Both span range of LB to UB at 0.01 intervals.
%For the current default values, Constraint 1 has 51 discrete values
%and Constraint 2 has 27 discrete values
    Cxi = floor(Constraints(1,1)*100)/100 : 0.01 : ceil(Constraints(1,2)*100)/100; 
    Cyj = floor(Constraints(2,1)*100)/100 : 0.01 : ceil(Constraints(2,2)*100)/100;
        Cyj = transpose(Cyj); %convert to column vector 

%[SR26] Eq. 8
%Solve for A for every combination of Constraints 1 and 2 and adg shapes

%Initialize array for holding values of A for each pair of Cxi,Cyj 
%applied to each of the adg shapes
    A = NaN(numel(Cyj),numel(Cxi),nshapes);

    for i = 1:numel(Cxi)
         for j = 1:numel(Cyj)
             x = Cxi(i); y = Cyj(j);
             for k = 1:nshapes
                 numA = ( (1-x)*(anw(idx490S)-y*anw(idx443S)) )  - ...
                        ( (1-y)*(anw(idx412S)-x*anw(idx443S)) ); %numerator of A in Eq. 8
                 denA = ( (1-x)*(shapelib.adg(idx490L,k)-y*shapelib.adg(idx443L,k)) )...
                        - ( (1-y)*(shapelib.adg(idx412L,k)-x*shapelib.adg(idx443L,k)) ); %denominator of A in Eq. 8
                 A(j,i,k) = numA/denA; %store in array
             end
         end
    end
    clear i j k x y numA denA; %discard unneeded intermediate variables
%Matrix of coefficient A is a 3-D array, Cyj -by- Cxi -by- nshapes

%Require that A is zero or positive; set any negative values of A to NaN
    A(A<0) = NaN;
    
%Generate the speculative solutions for each combination at the specific 7
%wavelengths required for evaluation of constraints. The full spectra will be
%generated only after infeasible solutions are eliminated by application of the
%remaining inequality Constraints 3 through 8.

%Calculate adgprime at specific wavelengths for each shape 
%Preallocate arrays for values at specific constraint wavelengths
    adgprime412 = NaN(size(A));adgprime443 = NaN(size(A));adgprime470 = NaN(size(A));
    adgprime490 = NaN(size(A));adgprime510 = NaN(size(A));adgprime555 = NaN(size(A));
    adgprime670 = NaN(size(A));
    
    for k=1:nshapes
        adgprime412(:,:,k) = A(:,:,k) .* repmat(shapelib.adg(idx412L,k),[size(A,1) size(A,2)]);
        adgprime443(:,:,k) = A(:,:,k) .* repmat(shapelib.adg(idx443L,k),[size(A,1) size(A,2)]);
        adgprime470(:,:,k) = A(:,:,k) .* repmat(shapelib.adg(idx470L,k),[size(A,1) size(A,2)]);
        adgprime490(:,:,k) = A(:,:,k) .* repmat(shapelib.adg(idx490L,k),[size(A,1) size(A,2)]);
        adgprime510(:,:,k) = A(:,:,k) .* repmat(shapelib.adg(idx510L,k),[size(A,1) size(A,2)]);
        adgprime555(:,:,k) = A(:,:,k) .* repmat(shapelib.adg(idx555L,k),[size(A,1) size(A,2)]);
        adgprime670(:,:,k) = A(:,:,k) .* repmat(shapelib.adg(idx670L,k),[size(A,1) size(A,2)]);
    end

%Calculate aphprime at each wavelength as difference between anw and adgprime
    aphprime412 = anw(idx412S) - adgprime412;
    aphprime443 = anw(idx443S) - adgprime443;
    aphprime470 = anw(idx470S) - adgprime470;
    aphprime490 = anw(idx490S) - adgprime490;
    aphprime510 = anw(idx510S) - adgprime510;
    aphprime555 = anw(idx555S) - adgprime555;
    aphprime670 = anw(idx670S) - adgprime670;

%These are all speculative solutions at specific wavelengths.
%With default constraints and adg shape library there is a maximum of
%51*27*498 = 685,746 speculative solutions (if all values of A are positive).
%Next step is to apply remaining 6 constraints to remove infeasible solutions

%==APPLY INEQUALITY CONSTRAINTS TO IDENTIFY AND ELIMINATE INFEASIBLE SOLUTIONS==
%[SR26] section 3E, Tables 1 & 2
%Each remaining constraint is tested on the speculative solutions sequentially,
%identifying and eliminating infeasible solutions through the process by setting
%the values of A to NaN.
%Once indices of feasible solutions (A and shape combinations) which satisfy all
%constraints are identified, full spectra of feasible solutions are then generated.

%Constraint 3 
    C3 = aphprime470./aphprime412;
    A(C3 < Constraints(3,1) | C3 > Constraints(3,2)) = NaN;               
    clear C3;

%Constraint 4
    C4 = aphprime555./aphprime490;
    A(C4 < Constraints(4,1) | C4 > Constraints(4,2)) = NaN; 
    clear C4;

%Constraint 5
    C5 = aphprime510./aphprime555;
    A(C5 < Constraints(5,1) | C5 > Constraints(5,2)) = NaN;
    clear C5;

%Constraint 6
    C6 = aphprime443./aphprime670;
    A(C6 < Constraints(6,1) | C6 > Constraints(6,2)) = NaN;
    clear C6;

%Constraint 7
    C7 = ( adgprime412 .* (adgprime443+aphprime443) ) ./ (adgprime412+aphprime412).^2 ;
    A(C7 < Constraints(7,1) | C7 > Constraints(7,2)) = NaN; 
    clear C7;

%Constraint 8
    aph510n = aphprime510./aphprime443;  aph555n = aphprime555./aphprime443;
    C8 = (aph510n - aph555n) ./ (555-510);
    A(C8 < Constraints(8,1) | C8 > Constraints(8,2)) = NaN; 
    clear C8 aph510n aph555n;

%Values of A that do not satisfy all 8 constraints are now set to NaN
%Get indices of all the feasible combinations of [A,shape#]

%Find ~isnan values and return indices of i, j, k for valid solutions
    idx_feasible = find(~isnan(A)); %index of all feasible values of A
    A_feasible = A(idx_feasible);   %the actual feasible values of A
 
    [~,~,A_feasible_k] = ind2sub(size(A),idx_feasible); 
    %A_feasible_k is the shape number corresponding to each feasible value of A
    %note: numel(unique(A_feasible_k)) will return the number of adg shapes 
    %contributing to the feasible pool

%Now have a list of all feasible values of A (A_feasible) and the corresponding
%adg shape number from the library (A_feasible_k) for each value of A

%Clear some variables to conserve memory before generating all feasible spectra
    clear aphprime* adgprime* A;

%======GENERATE THE POOL OF FEASIBLE SOLUTIONS AT FULL SPECTRAL RESOLUTION=====
%Create a matrix with all the feasible solutions of adg (adgprime)
%at full spectral resolution
    nlambdaL = size(shapelib.lambda,1); %number of wavelengths in the adg shape library

%Multiply the feasible A values by the corresponding adg shape number to get
%each feasible solution of adg at full spectral resolution of library
    adg_feasible = repmat(transpose(A_feasible),nlambdaL,1) .* shapelib.adg(:,A_feasible_k);
    %all feasible adg solutions; rows = lambda, cols = feasible spectra

%Determine the feasible solutions of aph as difference between anw and all
%feasible solutions of adg; i.e. aph = anw - adg
%Calculations of aph can only be done at the wavelengths of input anw spectrum
%aphidxL is the list of shapelib lambda indices which match input anw
    [~,aphidxL] = ismember(Lambda_anw,shapelib.lambda); 
    aph_feasible = repmat(anw,1,numel(A_feasible)) - adg_feasible(aphidxL,:); 

%Include a constraint that absorption cannot be negative.
%Slightly negative values of absorption coefficient sometimes occur, mostly 
%for wavelengths >690 nm (measurement noise region). Rather than eliminate
%entire spectrum, these negative values are set to zero.
%Will not strictly satisfy closure of anw = adg + aph for these individual
%spectra and wavelengths, but should have little influence in 
%determination of optimal solutions
    adg_feasible(adg_feasible<0)= 0;
    aph_feasible(aph_feasible<0)= 0;

%Include a constraint that modeled adg cannot be greater than input anw
%at any wavelength, these values are removed from the feasible solution
%pool
    adg_feasible(adg_feasible > repmat(anw,1,numel(A_feasible))) = NaN;

%=====DETERMINE OPTIMAL SOLUTIONS AND PERCENTILE VALUES FOR OUTPUT=====
%The optimal solution is calculated based on the pool of feasible solutions
%of adg, and the corresponding optimal solution of aph is then determined
%as the difference between measured anw and optimal adg. This means that
%output adg will always be 350 to 700 nm at 1-nm resolution (same as adg
%shape library), but optimal aph spectrum will only include wavelengths of
%input anw.

%[SR26] section E & Fig. 8
%Based on analysis, the optimal solution for adg is chosen as the 85th percentile
%value from the pool of all feasible solutions. By default, the output also 
%includes the 50th and 99th percentile values of adg which provides a range of 
%feasible solutions that is expected to contain the true solution with high probability
    adgopt(:,1) = prctile(adg_feasible,85,2); %85th percentile as first column of output
        adgopt(:,2:nPT+1) = prctile(adg_feasible, PT, 2); %additional percentiles
        aphopt = repmat(anw, 1, nPT+1) - adgopt(aphidxL,:);
 
%Define other function outputs 
    Lambda_out = shapelib.lambda; %Wavelengths corresponding to the adg shape library
    N_feasible = size(adg_feasible,2); %return the Number of feasible solutions
