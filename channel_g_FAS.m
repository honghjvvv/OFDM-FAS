classdef channel_g_FAS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % hanjiang.hong
    % H. Hong, K. K. Wong, X. Zhu, C. Rao, D. He, and H. Shin, ``Phased
    % ultra massive array (PUMA),'' submit to {\em IEEE Trans. Wireless
    % Commun.}, 2026.
    % 2024-02
    % 2024-03-01 update channel model 3!!
    % 2024-05-03 update channel model 4, the 2D-FA channel.
    % 2024-05-07 update channel model 5, the finite-scatterer channel
    % 2024-05-09 channel model 6, co-located finite-scatterer channel
    % 2025-03-19 6不对哦, channel model 5 N = 1, dk inf-->
    % 2025-03-05 N = [1,1] for channel nodel 4
    % 2025-07-17 update 5: when N_1 ~= N_2
    % 2026-01-15 J_0 -> j_0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        type            % channel type: 1 - variant mu
                        %               2 - constant mu
                        %               3 - accurate analytical model
                        %               4 - accurate channel for 2D-FA
                        %               5 - finite-scatterer channel
        N               % number of port
        W               % width
        sigma           % 
        g               % the channel
        % parameters for type 1
        mu_k
        % parameters for type 2
        mu
        % parameters for type 3 & 4
        U
        sqrt_lamda
        SSS
        % parameters of type 5 & 6
        K
        Np
        dev_theta
        dev_phi
        dk1      % k x d
        dk2
        theta_0
        phi_0
        theta
        phi
    end
    methods
        function obj = channel_g_FAS(varargin)
            if numel(varargin) == 3
                obj.type = 1;
                obj.N = varargin{1};
                obj.W = varargin{2};
                obj.sigma = sqrt(varargin{3});
            elseif numel(varargin) == 4
                obj.type = varargin{4};
                obj.N = varargin{1};
                obj.W = varargin{2};
                obj.sigma = sqrt(varargin{3});
            elseif numel(varargin) == 6
                obj.type = varargin{4};
                obj.N = varargin{1};
                obj.W = varargin{2};
                obj.sigma = sqrt(varargin{3});
                obj.K = varargin{5};
                obj.Np = varargin{6};
                if 5 ~= obj.type && 6 ~= obj.type
                    error('Input error')
                end
            else
                error('Input Error: input {N, W, sigma^2} or {N, W, sigma^2, type}, or {N, W, Sigma^2, 5/6, K, Np}');
            end
            if size(obj.N, 2) == 2
                if obj.N(1) == 1 && obj.N(2) ~= 1
                    obj.N(1) = [];
                    obj.W(1) = [];
                    obj.type = 3;
                elseif obj.N(1) ~= 1 && obj.N(2) == 1
                    obj.N(2) = [];
                    obj.W(2) = [];
                    obj.type = 3;
                end
            end
            if 1 == obj.type
                obj.mu_k = (bessej(0,(0:obj.N-1)*2*pi/(obj.N-1)*obj.W)).';
            elseif 2 == obj.type
                obj.mu = sqrt(2)*...
                    sqrt(...
                    hypergeom(1/2,[1,3/2],-pi^2*obj.W^2)-besselj(1,2*pi*obj.W)/2/pi/obj.W...
                    );
            elseif 3 == obj.type
                % SSS = besselj(0,(repmat((1:obj.N).',1,obj.N)-repmat(1:obj.N,obj.N,1))*2*pi/(obj.N-1)*obj.W);
                obj.SSS = toeplitz(besselj(0,(0:obj.N-1)*2*pi/(obj.N-1)*obj.W));
                % SSS = toeplitz(sinc((0:obj.N-1)*2*pi/(obj.N-1)*obj.W));
                [obj.U,S,V] = svd(obj.SSS);
                obj.sqrt_lamda = sqrt(diag(S));
            elseif 4 == obj.type
                if length(obj.N) ~= 2 || length(obj.W) ~= 2
                    error('The input of channel model 4 is wrong, N = [N1, N2], W = [W1, W2]')
                end
                N1 = obj.N(1);
                N2 = obj.N(2);
                NN = N1 * N2;
                W1 = obj.W(1);
                W2 = obj.W(2);
                SS = zeros(N2,N2,N1);
                % if N1==1
                    % SS = 1;
                % else
                for idx = 1:N1
                    zzz = 2*pi*sqrt(((idx-1)/(N1-1)*W1)^2+((0:N2-1)/(N2-1)*W2).^2);
                    SS(:,:,idx) = toeplitz(sphbesselj(0,zzz));
                    % SS(:,:,idx) = toeplitz(besselj(0,2*pi*sqrt(((idx-1)/(N1-1)*W1)^2+((0:N2-1)/(N2-1)*W2).^2)));
                end
                % end
                SSS2 = reshape(SS(:,:,toeplitz(1:N1)),N2,NN*N1,1);
                obj.SSS = zeros(NN,NN);

                for idx = 1:N2
                    obj.SSS(:,idx:N2:end) = reshape(SSS2(:,idx:N2:end),NN,N1);
                end
                [obj.U,S,V] = svd(obj.SSS);
                obj.sqrt_lamda = sqrt(diag(S));
            elseif 5 == obj.type
                obj.dev_theta = 5/180*pi;
                obj.dev_phi = 5/180*pi;
                if obj.N(1)==1
                    obj.dk1 = 2 * pi * obj.W(1);
                elseif obj.N(1)>1
                    obj.dk1 = 2 * pi * obj.W(1) / (obj.N(1) - 1);
                end
                if obj.N(2)==1
                    obj.dk2 = 2 * pi * obj.W(1);
                elseif obj.N(1)>1
                    obj.dk2 = 2 * pi * obj.W(2) / (obj.N(2) - 1);
                end
            else
                error('Channel type Error: must be 1, 2, 3, 4, or 5')
            end
        end
        function obj = update(obj)
            if 1 == obj.type
                if isempty(obj.mu_k)
                    error('Error')
                end
                g1 = obj.sigma*(randn+1j*randn)/sqrt(2);
                obj.g = obj.sigma * (sqrt(1-obj.mu_k.^2).*...
                    (randn(obj.N,1)+1j*randn(obj.N,1))/sqrt(2)+obj.mu_k*g1);
            elseif 2 == obj.type
                if isempty(obj.mu)
                    error('Error')
                end
                g1 = obj.sigma*(randn+1j*randn)/sqrt(2);
                obj.g = obj.sigma * (sqrt(1-obj.mu.^2).*...
                    (randn(obj.N,1)+1j*randn(obj.N,1))/sqrt(2)+obj.mu*g1);
            elseif 3 == obj.type
                if (isempty(obj.U) || isempty(obj.sqrt_lamda))
                    error('Error')
                end
                % obj.g = obj.sigma*(randn(obj.N,1)+1j*randn(obj.N,1))/sqrt(2).*(obj.U*obj.sqrt_lamda);
                obj.g = obj.sigma*obj.U*diag(obj.sqrt_lamda)*((randn(obj.N,1)+1j*randn(obj.N,1))/sqrt(2));
            elseif 4 == obj.type
                if (isempty(obj.U) || isempty(obj.sqrt_lamda))
                    error('Error')
                end
                NN = obj.N(1)*obj.N(2);
%                 if NN > 1
                    obj.g = obj.sigma*obj.U*diag(obj.sqrt_lamda)*((randn(NN,1)+1j*randn(NN,1))/sqrt(2));
%                 elseif NN == 1
%                     obj.g = (randn(NN,1)+1j*randn(NN,1))/sqrt(2);
%                 end
            elseif 5 == obj.type
                m = (0:(obj.N(1)-1)).';
                n = 0:(obj.N(2)-1);
                sigma_2 = obj.sigma^2;
                a_variance_LOS = (obj.K*sigma_2/(obj.K+1));
                a_LOS = sqrt(a_variance_LOS)*exp(1j*rand*2*pi);
                obj.theta_0 = (rand-1/2)*pi;% + randlap([1,1],obj.dev_theta);
                obj.phi_0 = rand*2*pi;% + randlap([1,1],obj.dev_phi);
                obj.g = a_LOS * reshape(exp(-1j*(obj.dk1*m*sin(obj.phi_0)*sin(obj.theta_0)+obj.dk2*n*cos(obj.theta_0))),[],1);
                
                a_variance_NLOS = (sigma_2/(obj.K+1))/obj.Np;
                a_NLOS = (randn(1,obj.Np)+1j*randn(1,obj.Np))*sqrt(a_variance_NLOS/2);
                obj.theta = (rand(1,obj.Np)-1/2)*pi;% +randlap([1,obj.Np],obj.dev_theta);
                obj.phi = rand(1,obj.Np)*2*pi;% + randlap([1,obj.Np],obj.dev_phi);
                
                ar = zeros(obj.N(1)*obj.N(2),obj.Np);
                for idx = 1:obj.Np
                    ar(:,idx) = reshape(exp(-1j*(obj.dk1*m*sin(obj.phi(idx))*sin(obj.theta(idx))+obj.dk2*n*cos(obj.theta(idx)))),[],1);
                    obj.g = obj.g + a_NLOS(idx) *ar(:,idx);
                end
            end
        end
    end
end

function x = randlap(SS,dev_sta)
% SS: size = [size1, size2]
% dev_sta: standard deviation
    d = sqrt(2)/dev_sta;
    b = dev_sta/sqrt(2); % 1/d
    exp_d_pi = exp(-d*pi);
    d_div_c = 2*(1-exp_d_pi);
    z = rand(SS);
    x = zeros(SS);
    pp1 = (z<=0.5);
    pp2 = (z>0.5);
    x(pp1) =b*log(d_div_c*z(pp1)+exp_d_pi);
    x(pp2) =-b*log(1-d_div_c*(z(pp2)-0.5)); 
end 

function y = sphbesselj(n,x)
    y = zeros(size(x));
    idx_zero = (x==0);
    if any(idx_zero)
        y(idx_zero & (n==0)) = 1;
        y(idx_zero & (n~=0)) = 0;
    end
    idx_nonzero = ~idx_zero;
    if any(idx_nonzero)
        y(idx_nonzero) = sqrt(pi./(2*x(idx_nonzero))).*besselj(n+0.5,x(idx_nonzero));
    end
end