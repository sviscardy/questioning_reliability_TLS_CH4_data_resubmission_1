% splab = SS_splab(selsp)
%
% This function converts the name of the selected species 'selsp' into a
% label in TeX style 'splab'. It can thus be used to make figures (e.g.
% label, title, legend).
% 
% Note that 'selsp' can be either in upper case or lower case, but must
% follow the way the names of species are defined in the chemical schemes
% (i.e. kpp_standard_chemistry, kpp_HDO_chemistry, etc.).
% 
% Exemples: - h2o and H2O are allowed
%           - h2og is NOT allowed.
% 
% Author: sebastien.viscardy@aeronomie.be (November 17, 2020)
% 
% Last update: sebastien.viscardy@aeronomie.be (April 19, 2022)
%%
function splab = SS_splab(selsp)

%% Main conversions
splab = char(upper(selsp));
splab = regexprep(splab,'CL','Cl');
splab = regexprep(splab,'AR','Ar');
splab = regexprep(splab,'1D','(^1D)');
splab = regexprep(splab,'4','_4');
splab = regexprep(splab,'5','_5');
splab = regexprep(splab,'6','_6');
splab = regexprep(splab,'PLUS','^+');

%% O2(1Delta)
if contains(splab,'2')
    splab = regexprep(splab,'2','_2');
    if contains(splab,'2D')
        splab = regexprep(splab,'_2D','_2(^1\\Delta)');
    end
end

%% O(3P)
if contains(splab,'3')
    splab = regexprep(splab,'3','_3');
    if contains(splab,'3P')
        splab = regexprep(splab,'_3P','(^3P)');
    end
end

%% H2Oice and HDOice --> H_2O_{(ice)} and HDO_{(ice)}
if contains(splab,'ICE')
    splab = regexprep(splab,'ICE','_{(s)}');
end

%% H2O --> H_2O_{(g)}
if contains(splab,'H_2O') && ~contains(splab,'CH_2O')
    kstr = strfind(splab,'H_2O');
    nk   = length(kstr);
    %     disp(splab)
    for ik = 1:nk
        kstr = strfind(splab,'H_2O');
        if length(splab)>=kstr(ik)+4
            if ~strcmpi(splab(kstr(ik)+4),'_')
                splab = [splab(1:kstr(ik)+3),'_{(g)}',splab(kstr(ik)+4:end)];
            end
        elseif length(splab)==kstr(ik)+3
            splab = [splab(1:end),'_{(g)}'];
        end
    end
end

%% HDO --> HDO_{(g)}
if contains(splab,'HDO')
    kstr = strfind(splab,'HDO');
    nk   = length(kstr);
    for ik = 1:nk
        kstr = strfind(splab,'HDO');
        if length(splab)>=kstr(ik)+3
            if ~strcmpi(splab(kstr(ik)+3),'_')
                splab = [splab(1:kstr(ik)+2),'_{(g)}',splab(kstr(ik)+3:end)];
            end
        elseif length(splab)==kstr(ik)+2
            splab = [splab(1:end),'_{(g)}'];
        end
    end
end