function gamma_eff = findEffectiveGamma(r_parent, r_daughters)
    % Solves r_parent^gamma = sum(r_daughters.^gamma) to find effective gamma
    % Inputs:
    %   r_parent    - scalar, radius of the parent vessel
    %   r_daughters - vector, radii of daughter vessels
    % Output:
    %   gamma_eff   - effective radial scaling exponent

    % Define the function to solve
    equation = @(gamma) r_parent.^gamma - sum(r_daughters.^gamma);

    % Use fsolve to find the root, with initial guess gamma = 2.5
    gamma_guess = 2.5;
    options = optimset('Display','off'); % suppress output
    gamma_eff = fsolve(equation, gamma_guess, options);
end