function effective_conductivity(Keff_file, h_range, user_data, visualise)
% ------------------------------------------------------------------------------
% effective_conductivity.m
%
% Written by Dr Elliot Carr (2013-2014)
% Ecole Centrale Paris and Queensland University of Technology
% 
% This code is part of TwoScalRich.
%
% Solves periodic cell problems computing the effective conductivity (Keff)
% for each discrete value of h in the vector 'h_range'.
%
% ------------------------------------------------------------------------------

% Soil A
KsatA  = user_data.soilA.Ksat;
alphaA = user_data.soilA.alpha;
nA     = user_data.soilA.n;
mA     = user_data.soilA.m;
Se     = (1 + (-alphaA*h_range).^nA).^(-mA);
KA     = KsatA * sqrt(Se) .* (1.0 - (1.0 - (Se.^(1.0/mA))).^mA).^2.0;
% Soil B
KsatB  = user_data.soilB.Ksat;
alphaB = user_data.soilB.alpha;
nB     = user_data.soilB.n;
mB     = user_data.soilB.m;
Se     = (1 + (-alphaB*h_range).^nB).^(-mB);
KB     = KsatB * sqrt(Se) .* (1.0 - (1.0 - (Se.^(1.0/mB))).^mB).^2.0;

% Solve periodic cell problems
no_hvalues = length(h_range); Keff_mat = zeros(no_hvalues,4);
fprintf('%% Calculating the effective hydraulic conductivity\n');
count = 1;
for i = 1:no_hvalues 
    fprintf('Solving periodic cell problems for h = %1.2e (%i of %i)\n',...
            h_range(i),count,no_hvalues);
    for j = 1:2 % columns of Keff
        Keff = solve_cell_problem(j,user_data,KA(i),KB(i));
        Keff_mat(i,2*j-1:2*j) = Keff';
    end 
    count = count + 1;
end

% Save values for linear interpolation of Keff to file
% (see main.cpp and read_effective_conductivity.cpp)
logh_values = log10(-h_range);
logr = log10(exp((log(-h_range(no_hvalues))-log(-h_range(1))) /(no_hvalues-1)));
logh_first = log10(-h_range(1));

fid = fopen(Keff_file, 'w');
fprintf(fid,'%2.15f\n',user_data.micro_mesh.epsilonA);
fprintf(fid,'%i\n',no_hvalues);
fprintf(fid,'%2.15f\n',logr);
fprintf(fid,'%2.15f\n',logh_first);
fprintf(fid,'%2.15f\n',logh_values);
for i = 1:no_hvalues
    fprintf(fid,'%2.15e %2.15e %2.15e %2.15e\n',Keff_mat(i,:));
end
fclose(fid);

% Plot results
if visualise   
    
    figure;
    title('Keff')
       
    h_range_plot = h_range * 100; % Convert from m to cm
    Keff_mat = Keff_mat*100*3600; % Convert to cm/per h
    min_K = min(min(Keff_mat));
    max_K = max(max(Keff_mat));
    
    for k = 1:4
        subplot(2,2,k)
        semilogx(h_range_plot,Keff_mat(:,k),'k.-','MarkerSize',12)
        legend('Keff','Location','NorthWest')
        switch k
            case 1, title('(1,1) entry of Keff') 
            case 2, title('(1,2) entry of Keff') 
            case 3, title('(2,1) entry of Keff')
            case 4, title('(2,2) entry of Keff')
        end
        axis([h_range_plot(1), h_range_plot(end), min_K, max_K]);
        xlabel('h [cm]')
    end
    
end



end

