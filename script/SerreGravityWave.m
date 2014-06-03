%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% -------------------------------------------------- %%%
%%% Dynamic gravity wave simulations to the one-layer  %%%
%%%           Serre-Green-Naghdi equations             %%%
%%% -------------------------------------------------- %%%
%%% As an example, we simulate thr head-on collision   %%%
%%%         of two different solitary waves            %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SerreGravityWave ()
    
    clear all, close all, format longE
    
    %%% Physical parameters:
    g = 1.0;        % gravity acceleration
    d = 1.0;        % undisturbed water depth
    d2 = d^2;
    
    %%% Numerical parameters:
    l = 75.0;               % half-length of the domain [-l, l]
    N = 2048;               % number of Fourier modes in discrete solution
    dx = 2*l/N;             % distance between two physical points
    x = (1-N/2:N/2)'*dx;    % physical space discretization
    
    %%% Definition of various pseudo-diff operators:
    k = [0:N/2 1-N/2:-1]'*pi/l;     % vector of wavenumbers
    kmax = max(k);
    Nvec = 1:N/2+1;
    
    % antialising treatment (2/3 rule)
    j = (N/4+2:N/4*3);              % the frequencies we make zero
    k(j) = 0;
    
    K2 = k.^2;                      % the 2nd derivative symbol x (-1)
    kinv = zeros(size(k));          % declaration of 1/k
    kinv(k ~= 0) = 1./k(k ~= 0);    % needed to pass to dimensional variable \eta
    
    % some constants appearing in the Serre equations:
    C13 = 1/3;                  % a coefficient in the Serre equations
    theta = 0.5;                % relaxation parameter
    Linv = 1./(1+C13*d2*K2);    % operator related to the linear dispersion relation
    
    omega = sqrt(g*d*K2.*Linv); % the dispersion relation of the Serre equations
    ominv = zeros(size(omega)); % declaration of 1/omega
    % this variable is needed to recover the dimensional velocity field
    ominv(omega ~= 0) = 1./omega(omega ~= 0);

    %%% Initial conditions specification:
    a1 = 0.2; a2 = 0.3; ampl = max(a1, a2);
    [eta, u] = InitCond (a1, a2);
    
    eta_hat = fft(eta);                 % Fourier transform of the free surface
    u_hat = fft(u);                     % Fourier transform of the depth-averaged velocity
    eta0 = eta_hat(1);                  % save the average value
    eta_hat(j) = 0;
    etax = real(ifft(1i*k.*eta_hat));   % x-derivative to compute q below
    
    u0 = u_hat(1);                      % save the average value
    u_hat(j) = 0;                       % dealiasing
    ux = real(ifft(1i*k.*u_hat));       % we compute the derivative u_x
    uxx = real(ifft(-K2.*u_hat));       % and u_xx to initialize the q :
    q = u - (d+eta).*alias(C13*(d+eta).*uxx + etax.*ux);
    q0 = mean(q);                       % save the average value
    
    % Finally the initial condition in Fourier space and dynamic variables:
    v = zeros(N,2);                     % vector of dynamic variables
    v(:,1) = 1i*k.*eta_hat;             % dimensionless free surface
    v(:,2) = 1i*omega.*fft(q)/g;        % and the potential vorticity
    
    %%% Parameters for the time stepping:
    t = 0.0;        % the discrete time variable
    T = 40.0;       % final simulation time
    M = 1000;       % number of time steps
    m = 10;         % we plot the solution every m steps
    dt = (T-t)/M;   % individual time step
    
    tol = 1e-14;                % tolerance for the fixed point iterations
    Times = (t:m*dt:T)';        % vector of time moments where we store the solution
    Energy = zeros(M/m+1,1);    % vector of the total energy values
    Eta = zeros(M/m+1,N);       % Snapshots of the free surface to make a space-time plot

    Energy(1) = Hamilt(eta, u); % store the initial energy of the system
    Eta(1,:) = eta;             % store the free surface shape at the initial time
    
    % set-up the size of the graphical window:
    scrsz = get(0,'ScreenSize');
    figure('Position', [scrsz(3)/4 0 scrsz(3)/2 3/4*scrsz(4)])

    % and let's plot the initial condition:
    Plot (eta, u, t);
    
    % main loop in time:
    for jj=1:M
        % make one time step of the Serre's solver:
        [v, u_hat] = SerreOneStep(v, u_hat, dt);
        
        % if we made already m steps, we plot the current numerical results:
        if (mod(jj,m) == 0)
          t = jj*dt; % current moment of time
          
          % we reconstruct the physical variables:
          eta_hat = -1i*kinv.*v(:,1); eta_hat(1) = eta0;
          eta = real(ifft(eta_hat)); h = d + eta;
          u = real(ifft(u_hat));
          
          frame = jj/m + 1;
          Energy(frame) = Hamilt(eta, u);   % store the total energy to check the accuracy
          Eta(frame,:) = eta;               % store the free surface shape
          
          Plot(eta, u, t);
        end % if (draw)
        
    end % for jj

    %%% At the end of the simulation we show the space-time plot:
    figure;
    % we take only a subset of saved snapshots:
    waterfall(x, Times(1:2:end), Eta(1:2:end,:)), grid on, hold off
    colormap('gray'), view([-13 58]);
    xlabel('$x$', 'interpreter', 'latex', 'fontsize', 16);
    ylabel('$t$', 'interpreter', 'latex', 'fontsize', 16);
    zlabel('$\eta(x,t)$', 'interpreter', 'latex', 'fontsize', 16);
    set(gcf, 'Color', 'w')

    %%% At the end of the simulation we show the energy conservation:
    figure;
    plot(Times, Energy - mean(Energy), 'bo-', 'LineWidth', 1.0), grid off
    xlabel('$t$', 'interpreter', 'latex', 'fontsize', 14);
    ylabel('$\mathcal{E}(t) - \langle\mathcal{E}\rangle$', 'interpreter', 'latex', 'fontsize', 14);
    title('Energy conservation');
    set(gcf, 'color', 'w');

    %%% END OF THE MAIN SCRIPT %%%
    %%% ---------------------------------------------------- %%%

    %%% We set-up an initial condition consisting of two 
    % counter-propagating solitary waves.
    %%% Input parameters:
    % a1 : left solitary wave amplitude
    % a2 : right solitary wave amplitude
    function [eta, u] = InitCond (a1, a2)
        c1 = sqrt (g*(d+a1));               % SW1 speed
        k1 = sqrt(3*a1/(a1+d))/d;           % SW1 decay rate
        c2 = sqrt (g*(d+a2));               % SW2 speed
        k2 = sqrt(3*a2/(a2+d))/d;           % SW2 decay rate

        eta1 = a1*sech(0.5*k1*(x+25)).^2;   % SW1
        eta2 = a2*sech(0.5*k2*(x-25)).^2;   % SW2

        eta = eta1 + eta2;
        h = d + eta;                        % total water depth
        u = (c1*eta1 - c2*eta2)./h;         % velocity field
    end % InitCond ()

    %%% ---------------------------------------------------- %%%

    % Plotting function of the current solution:
    function Plot (eta, u, t)
        % Plot the current snapshot:
        subplot(3,1,1)
        plot (x, eta, 'b-', 'LineWidth', 1.4), grid off, hold on
        plot (x, ampl*ones(size(x)), 'k--'), hold on
        plot (x, (a1 + a2)*ones(size(x)), 'r--'), hold off
        axis([-l l -0.2*ampl 2.2*ampl])
        legend('Free surface', 'Heighest wave ampl.', 'Sum of two ampl.',...
            'location', 'NorthEast');
        xlabel('$x$', 'interpreter', 'latex', 'fontsize', 14);
        ylabel('$\eta(x,t)$', 'interpreter', 'latex', 'fontsize', 14);
        title(['Free surface \eta(x,t) at t = ',num2str(t,'%5.2f'),...
            '; dt = ',num2str(dt,'%4.3f')])
          
        subplot(3,1,2)
        plot (x, u, 'k-', 'LineWidth', 1.4), grid off
        axis([-l l -1.2*ampl 1.2*ampl])
        xlabel('$x$', 'interpreter', 'latex', 'fontsize', 14);
        ylabel('$u(x,t)$', 'interpreter', 'latex', 'fontsize', 14);
        title('Horizontal velocity')
        
        subplot(3,1,3)
        eta_hat = fft(eta);
        loglog(k(Nvec), abs(eta_hat(Nvec,1)).^2, 'r-', 'LineWidth', 1.4), grid on
        axis([1 0.9*kmax 1e-35 1e5])
        xlabel('$k$', 'interpreter', 'latex', 'fontsize', 14);
        ylabel ('$|\hat{\eta}(k,t)|^2$', 'interpreter', 'latex', 'fontsize', 14);
        title('Fourier spectrum');

        set(gcf, 'color', 'w');
        drawnow
    end % Plot ()
    
    %%% Function which computes the total energy in the Serre system:
    % Input parameters (both in the real space):
    %  eta: free surface elevation
    %    u: horizontal velocity
    % Output parameters:
    %    H: value of the total energy (= Hamiltonian in this case)
    function H = Hamilt(eta, u)
        u_hat = fft(u); ux = real(ifft(1i*k.*u_hat));
        h = d + eta;
        etax = real(ifft(1i*k.*fft(eta)));
        
        H = 0.5*dx*sum(h.*u.^2 + C13*(h.^3).*ux.^2 + g*eta.^2);
    end % Hamilt()    
    
    %%% ---------------------------------------------------- %%%
    
    %%% Dealiasing function:
    % Input parameters:
    %    fr: nonlinear signal in the real space
    % Output parameters:
    %  filt: the filtered signal
    function filt = alias(fr)
      f_hat = fft(fr);
      f_hat(j) = 0;
      filt = real(ifft(f_hat));
    end % alias()
    
    %%% ---------------------------------------------------- %%%
    
    function w = turn(v, dt)
        w = zeros(N,2);
        
        cosw = cos(omega*dt);
        sinw = sin(omega*dt);
        
        w(:,1) =    cosw.*v(:,1) + 1i*sinw.*v(:,2);
        w(:,2) = 1i*sinw.*v(:,1) +    cosw.*v(:,2);
        w(j,:) = 0; % antialiasing
    end % turn()

    %%% ---------------------------------------------------- %%%
    
    %%% Right-hand side of the Serre equations:
    % Input parameters:
    %    w : vector of dynamic variables (i*k*eta, i*omega*q/g) in Fourier space
    %   u0 : velocity field on the previous time step (as the initial guess)
    % Output parameters:
    %   rhs: the right-hand side
    % u_hat: Fourier transform of the reconstructed horizontal velocity
    function [rhs, u_hat] = RHS(w, u0_hat, dt)
        
        % declaration of the result and memory preallocation
        rh = zeros(N,2);    % RHS without rotation
        v = turn(w, -dt);   % we turn back the linear integration
        
        % return to dimensional variables
        eta_hat = -1i*kinv.*v(:,1);     % * 1/(i*k)
        eta_hat(1) = eta0;
        q_hat = -1i*g*ominv.*v(:,2);    % * g/(i*omega)
        q_hat(1) = q0;
        
        qr = real(ifft(q_hat));             % potential vorticity in real space
        
        eta = real(ifft(eta_hat));          % we recover \eta in the real space
        etave = sum(eta)/N;                 % average value of \eta
        eprim = eta - etave;                % \eta' = \eta - \eta_ave
        etax = real(ifft(1i*k.*eta_hat));   % \eta_x in real space
        etaxx = real(ifft(-K2.*eta_hat));   % \eta_x in real space
        h = d + eta;                        % the total water depth
        
        err = inf;
        while (err > tol) % in this loop we recover u from q
            u0x = real(ifft(1i*k.*u0_hat)); % u0_x
            u0xx = real(ifft(-K2.*u0_hat)); % u0_xx
            u_hat = theta*(q_hat + fft(h.*alias(etax.*u0x) +...
                C13*(alias(eprim.^2) + 2*(d+etave)*eprim).*u0xx))./(1+C13*(d+etave)^2*K2) +...
                (1-theta)*u0_hat;           % relaxation term
            u_hat(j) = 0;                   % dealiasing
            % compute difference between iterations:
            err = norm(real(ifft(u_hat-u0_hat)), inf);
            u0_hat = u_hat;                 % initialisation for the next step
        end % while (err > tol)
        
        u_hat(1) = u0;
        u = real(ifft(u_hat));          % recovered velocity
        ux = real(ifft(1i*k.*u_hat));   % u_x in phys. space
        uxx = real(ifft(-K2.*u_hat));   % u_xx in phys. space
        
        hux = alias(h.*ux); % combination appearing in q_t equation
        
        % Difference between q and u_linear
        FNL = fft(C13*(alias(eta.^2) + 2*d*eta).*uxx + hux.*etax);
        
        % now we assemble the RHS
        rh(:,1) = K2.*fft(eta.*u) + d*K2.*FNL./(1+C13*d^2*K2);
        rh(:,2) = -omega.*k.*fft(0.5*u.^2 + 0.5*hux.^2 - u.*qr)/g;

        rhs = turn(rh, dt); % we rotate the RHS and return the result
    end % RHS ()
    
    %%% ---------------------------------------------------- %%%
    
    %%% This function makes one adaptive time step:
    % Input parameters:
    %   v: solution in Fourier space at the previous time step
    %   u: depth-averaged velocity at the previous time step
    %  dt: the time step value
    % Output parameters:
    %   w: solution in Fourier space at the following time step
    %   u: depth-averaged velocity at the next time step
    %%% Time-stepping: This Runge-Kutta-Prince-Dormand method 4-5
    function [v, u] = SerreOneStep(v, u, dt)
        [k1, u] = RHS(v, u, 0);
        [k2, u] = RHS(v + dt*k1/5, u, dt/5);
        [k3, u] = RHS(v + dt/40*(3*k1 + 9*k2), u, 3*dt/10);
        [k4, u] = RHS(v + dt/10*(3*k1 - 9*k2 + 12*k3), u, 3*dt/5);
        [k5, u] = RHS(v + dt/729*(226*k1 - 675*k2 + 880*k3 + 55*k4), u, 2*dt/3);
        [k6, u] = RHS(v + dt/2970*(-1991*k1 + 7425*k2 - 2660*k3 - 10010*k4 + 10206*k5), u, dt);
        w = v + dt*(31/540*k1 + 190/297*k3 - 145/108*k4 + 351/220*k5 + 1/20*k6);
        v = turn(w, -dt);
    end % SerreOneStep()

    %%% ---------------------------------------------------- %%%

end % SerreGravityWave ()