function [x, fval, exitflag, output, grad, hessian] =...
    fit_glm_nested(dt, rt, B, eta, shrink_idx, x0) 
%     [x, fval, exitflag, output, grad, hessian] =...
%         fit_glm_nested(dt, rt, B, eta, shrink_idx)

    % Options    
    options = optimoptions(@fminunc,'GradObj', 'on',...
        'Hessian', 'on', 'Display', 'iter');

    % Initialize
    if nargin < 6
        x0 = zeros(size(B, 2), 1);
    end

    % Optimize
    [x, fval, exitflag, output, grad, hessian] =...
        fminunc(@nestedfun, x0, options);
    
    % Define nested function for passing parameters
    function [f, g, H] = nestedfun(x)
        [f, g, H] = glm_neg_log_post(x, dt, rt, B, eta, shrink_idx);
    end
end