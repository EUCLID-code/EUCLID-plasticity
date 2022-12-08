classdef lib_VEVPHardMixLin_v2
    % lib_VEVPHardMixLin_v2
    % The material behavior depends on the material parameters (e.g., `G0`,
    % `sigma_0`). All material parameters are collected in the parameter
    % vector `theta`. Note, however, that some parameters are stored in
    % `theta` with their reciprocals.
    %

    properties
        name
        n_Maxwell
        n_Maxwell_G
        n_Maxwell_K
        theta_string
        use_reciprocal
        theta_min
        theta_threshold
    end

    methods
        %% Constructor
        function obj = lib_VEVPHardMixLin_v2(n_Maxwell)
            obj.name = 'lib_VEVPHardMixLin_v2';
            obj.n_Maxwell = n_Maxwell;
            obj.n_Maxwell_G = n_Maxwell(1);
            obj.n_Maxwell_K = n_Maxwell(2);
            obj = define_theta_string(obj);
            obj = define_use_reciprocal(obj);
            obj = define_theta_bounds(obj);
            obj.theta_threshold = 1e-4;
        end

        %% Definitions
        function obj = define_theta_string(obj)
            obj.theta_string = {};
            obj.theta_string{end+1} = 'G0';
            obj.theta_string{end+1} = 'K0';
            for i = 1:obj.n_Maxwell_G
                obj.theta_string{end+1} = ['G' num2str(i)];
                obj.theta_string{end+1} = ['1 / g' num2str(i)];
            end
            for i = 1:obj.n_Maxwell_K
                obj.theta_string{end+1} = ['K' num2str(i)];
                obj.theta_string{end+1} = ['1 / k' num2str(i)];
            end
            obj.theta_string{end+1} = 'H_iso';
            obj.theta_string{end+1} = 'H_kin';
            obj.theta_string{end+1} = 'eta';
            obj.theta_string{end+1} = '1 / sigma_0';
        end

        function obj = define_use_reciprocal(obj)
            obj.use_reciprocal = zeros(size(obj.theta_string),'logical');
            for idx = 1:length(obj.use_reciprocal)
                if contains(obj.theta_string(idx),'1 /') || contains(obj.theta_string(idx),'1/')
                    obj.use_reciprocal(idx) = true;
                end
            end
        end

        function obj = define_theta_bounds(obj)
            theta_min_value = 1e-6;
            obj.theta_min = zeros(size(obj.theta_string));
            obj.theta_min(1) = theta_min_value; % (G0) must not be zero
            obj.theta_min(2) = theta_min_value; % (K0) must not be zero

            obj.theta_min(2+(2:2:2*obj.n_Maxwell_G)) = theta_min_value; % (1 / gi) must not be zero
            obj.theta_min(2+2*obj.n_Maxwell_G+(2:2:2*obj.n_Maxwell_K)) = theta_min_value; % (1 / ki) must not be zero

            obj.theta_min(end-1) = theta_min_value; % (eta) must not be zero
            obj.theta_min(end) = theta_min_value; % (1 / sigma_0) must not be zero
        end

        %% Functions
        function n = n_theta(obj)
            % n_theta counts the number of parameters
            n = length(obj.theta_string);
        end

        function theta = get_theta(obj,G0,K0,Gi,gi,Ki,ki,H_iso,H_kin,eta,sigma_0)
            % get_theta for given parameters
            theta_min_value = 1e-6;
            gi(gi < theta_min_value) = theta_min_value;
            ki(ki < theta_min_value) = theta_min_value;
            gi_reciprocal = 1./gi;
            ki_reciprocal = 1./ki;
            theta = [G0,K0,zip(Gi,gi_reciprocal),zip(Ki,ki_reciprocal),H_iso,H_kin,eta,1/sigma_0];
        end

        function theta = get_theta_extended(obj,n_Maxwell,G0,K0,Gi,gi,Ki,ki,H_iso,H_kin,eta,sigma_0)
            % get_theta for given parameters

            n_Maxwell_given = [length(Gi), length(Ki)];
            n_Maxwell_diff = n_Maxwell - n_Maxwell_given;
            Gi = [Gi, zeros(1,n_Maxwell_diff(1))];
            gi = [gi, zeros(1,n_Maxwell_diff(1))];
            Ki = [Ki, zeros(1,n_Maxwell_diff(2))];
            ki = [ki, zeros(1,n_Maxwell_diff(2))];
            theta = obj.get_theta(G0,K0,Gi,gi,Ki,ki,H_iso,H_kin,eta,sigma_0);
        end

        function [G0,K0,Gi,gi,Ki,ki,H_iso,H_kin,eta,sigma_0] = get_parameters(obj,theta)
            % get_parameters for given theta
            counter = 1;
            G0 = theta(counter); counter = counter + 1;
            K0 = theta(counter); counter = counter + 1;
            Gi = zeros(1,obj.n_Maxwell_G);
            gi = zeros(1,obj.n_Maxwell_G);
            for i = 1:obj.n_Maxwell_G
                Gi(i) = theta(counter); counter = counter + 1;
                gi(i) = 1 / theta(counter); counter = counter + 1;
            end
            Ki = zeros(1,obj.n_Maxwell_K);
            ki = zeros(1,obj.n_Maxwell_K);
            for i = 1:obj.n_Maxwell_K
                Ki(i) = theta(counter); counter = counter + 1;
                ki(i) = 1 / theta(counter); counter = counter + 1;
            end
            H_iso = theta(counter); counter = counter + 1;
            H_kin = theta(counter); counter = counter + 1;
            eta = theta(counter); counter = counter + 1;
            sigma_0 = 1 / theta(counter);
        end

        function material = get_material(obj,theta)
            % get_material for given theta
            [G0,K0,Gi,gi,Ki,ki,H_iso,H_kin,eta,sigma_0] = obj.get_parameters(theta);
            material.name = obj.name;
            material.G0 = G0;
            material.K0 = K0;
            material.Gi = Gi;
            material.gi = gi;
            material.Ki = Ki;
            material.ki = ki;
            material.H_iso = H_iso;
            material.H_kin = H_kin;
            material.eta = eta;
            material.sigma_0 = sigma_0;
            material.n_Maxwell_G = length(Gi);
            material.n_Maxwell_K = length(Ki);
            material.Ginf = material.G0 - sum(material.Gi);
            material.Kinf = material.K0 - sum(material.Ki);
            material.eta_Gi = material.Gi.*material.gi; % viscosities
            material.eta_Ki = material.Ki.*material.ki; % viscosities
        end

        function disp_theta(obj,theta)
            % disp_theta displays all parameters
            for idx = 1:obj.n_theta()
                disp([obj.theta_string{idx} ' = ' num2str(theta(idx))])
            end
        end

        function str = get_str_material_type(obj,theta)
            if obj.n_Maxwell_G ~= 1 || obj.n_Maxwell_K ~= 1
                error('Not yet implemented.')
            end
            str = ''; % string
            ans_no = '& \cellcolor{red!15} off ';
            ans_yes = '& \cellcolor{green!15} on ';
            ans_0 = '& \cellcolor{red!15} 0 ';
            ans_1 = '& \cellcolor{green!15} 1 ';
            % elasticity
            if theta(1) == 0 && theta(2) == 0
                str = [str , ans_no];
            else
                str = [str , ans_yes];
            end
            % viscoelasticity
            if theta(3) == 0 && theta(5) == 0
                str = [str , ans_no];
            else
                str = [str , ans_yes];
            end
            % number of Maxwell elements (shear)
            if theta(3) == 0
                str = [str , ans_0];
            else
                str = [str , ans_1];
            end
            % number of Maxwell elements (bulk)
            if theta(5) == 0
                str = [str , ans_0];
            else
                str = [str , ans_1];
            end
            % plasticity
            if theta(10) == 0
                str = [str , ans_no];
            else
                str = [str , ans_yes];
            end
            % viscoplasticity
            if theta(9) == 0
                str = [str , ans_no];
            else
                str = [str , ans_yes];
            end
            % isotropic hardening
            if theta(7) == 0
                str = [str , ans_no];
            else
                str = [str , ans_yes];
            end
            % kinematic hardening
            if theta(8) == 0
                str = [str , ans_no];
            else
                str = [str , ans_yes];
            end
        end

        function str = num2str0(obj,num)
            if num == 0
                str = '0';
            else
                str = num2str(num,'%.4f');
            end
        end

        function str = get_str_material_parameters(obj,theta)
            if obj.n_Maxwell_G ~= 1 || obj.n_Maxwell_K ~= 1
                error('Not yet implemented.')
            end
            str = ''; % string

            counter = 1;
            G0 = theta(counter); counter = counter + 1;
            str = [str , '& ' , obj.num2str0(G0), ' '];
            K0 = theta(counter); counter = counter + 1;
            str = [str , '& ' , obj.num2str0(K0), ' '];
            Gi = zeros(1,obj.n_Maxwell_G);
            gi = zeros(1,obj.n_Maxwell_G);
            for i = 1:obj.n_Maxwell_G
                Gi(i) = theta(counter); counter = counter + 1;
                str = [str , '& ' , obj.num2str0(Gi(i)), ' '];
                if theta(counter) == 0
                    gi(i) = Inf; counter = counter + 1;
                    str = [str , '& ' , '$\rightarrow\infty$', ' '];
                else
                    gi(i) = 1 / theta(counter); counter = counter + 1;
                    str = [str , '& ' , obj.num2str0(gi(i)), ' '];
                end
            end
            Ki = zeros(1,obj.n_Maxwell_K);
            ki = zeros(1,obj.n_Maxwell_K);
            for i = 1:obj.n_Maxwell_K
                Ki(i) = theta(counter); counter = counter + 1;
                str = [str , '& ' , obj.num2str0(Ki(i)), ' '];
                if theta(counter) == 0
                    ki(i) = Inf; counter = counter + 1;
                    str = [str , '& ' , '$\rightarrow\infty$', ' '];
                else
                    ki(i) = 1 / theta(counter); counter = counter + 1;
                    str = [str , '& ' , obj.num2str0(ki(i)), ' '];
                end
            end
            if theta(end) == 0
                sigma_0 = Inf;
                str = [str , '& ' , '$\rightarrow\infty$', ' '];
            else
                sigma_0 = 1 / theta(end);
                str = [str , '& ' , obj.num2str0(sigma_0), ' '];
            end
            eta = theta(end-1);
            str = [str , '& ' , obj.num2str0(eta), ' '];
            H_iso = theta(end-3);
            str = [str , '& ' , obj.num2str0(H_iso), ' '];
            H_kin = theta(end-2);
            str = [str , '& ' , obj.num2str0(H_kin), ' '];
        end

        function theta = apply_threshold(obj,theta)
            if obj.n_Maxwell_G ~= 1 || obj.n_Maxwell_K ~= 1
                error('Not yet implemented.')
            end
            theta(abs(theta) < obj.theta_threshold) = 0;
            % first Maxwell element
            if theta(3)*theta(4) < obj.theta_threshold
                theta(3) = 0;
                theta(4) = 0;
            end
            if theta(5)*theta(6) < obj.theta_threshold
                theta(5) = 0;
                theta(6) = 0;
            end
            % viscoplasticity
            if theta(end) < obj.theta_threshold
                theta(end-1) = 0;
            end
        end

    end
end













