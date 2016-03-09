%% Setup
% Problem Data
a = 1;
c = -1;
L = 1;

% Inputs
all_num_elems = [4 2 2 4];
all_elem_orders = [1 2 2 2];
all_alpha = [0 1/2 1/3 1/2];
names = {'4L', '2Q0.5', '2Q0.3', '4Q'};

for i = 1:length(all_num_elems)

    %% Meat of the code
    % Fetch data
    num_elems = all_num_elems(i);
    elem_order = all_elem_orders(i);
    alpha = all_alpha(i);
    % Make connectivity matrix
    num_nodes = num_elems*(elem_order+1) - (num_elems - 1);
    B = zeros(num_elems, elem_order+1);
    if (elem_order == 1)
        B(:,1) = 1:num_nodes-1;
        B(:,2) = 2:num_nodes;
    else
        B(:,1) = 1:2:num_nodes-2;
        B(:,2) = 2:2:num_nodes-1;
        B(:,3) = 3:2:num_nodes;
    end
    % Generate grid
    h = L/num_elems;  % Constant
    if (elem_order == 1)
        node_locs = 0:h:L;
    else
        node_locs = zeros(1,num_nodes);
        node_locs(1:2:end) = 0:h:L;
        node_locs(2:2:end) = node_locs(1:2:end-1) + alpha*h;
    end

    % Construct K, F
    K = zeros(num_nodes, num_nodes);
    F = zeros(num_nodes, 1);
    for elem = 1:num_elems
        % Construct stiffness matrix and source term
        sctr = B(elem,:);
        xa = node_locs(sctr(1));
        if (elem_order == 1)
            Kelem = a/h*[1 -1; -1 1] + c*h/6*[2 1; 1 2];
            Felem = zeros(2,1);
            Felem(1) = -h/12*(h^2 + 4*h*xa + 6*xa^2);
            Felem(2) = -h/12*(3*h^2 + 8*h*xa + 6*xa^2);
        else
            Kelem(1,1) = a/h + a/(3*alpha^2*h) + c*h*(10*alpha^2 - 5*alpha + 1)/(30*alpha^2);
            Kelem(1,2) = a/(3*alpha^2*h*(alpha - 1)) - c*h*(5*alpha - 2)/(60*alpha^2*(alpha - 1));
            Kelem(1,3) = -a*(3*alpha^2 - 3*alpha + 1)/(3*alpha*h*(alpha - 1)) + c*h*(10*alpha^2 - 10*alpha + 3)/(60*alpha*(alpha - 1));
            Kelem(2,2) = a/(3*alpha^2*h*(alpha^2 - 2*alpha + 1)) + c*h/(30*alpha^2*(alpha^2 - 2*alpha + 1));
            Kelem(2,3) = -a/(3*alpha*h*(alpha^2 - 2*alpha + 1)) - c*h*(5*alpha - 3)/(60*alpha*(alpha^2 - 2*alpha + 1));
            Kelem(3,3) = a*(3*alpha^2 - 6*alpha + 4)/(3*h*(alpha^2 - 2*alpha + 1)) + c*h*(10*alpha^2 - 15*alpha + 6)/(30*(alpha^2 - 2*alpha + 1));
            Felem(1) = -h*(5*alpha*h^2 + 20*alpha*h*xa + 30*alpha*xa^2 - 3*h^2 - 10*h*xa - 10*xa^2)/(60*alpha);
            Felem(2) = h*(3*h^2 + 10*h*xa + 10*xa^2)/(60*alpha*(alpha - 1));
            Felem(3) = -h*(30*alpha*xa^2 - 12*h^2 + 15*h*(alpha*h - 2*xa) + 20*xa*(2*alpha*h - xa))/(60*alpha - 60);
            Kelem(2,1) = Kelem(1,2); Kelem(3,1) = Kelem(1,3); Kelem(3,2) = Kelem(2,3);
        end
        % Scatter to global
        K(sctr,sctr) = K(sctr,sctr) + Kelem;
        F(sctr) = F(sctr) + Felem;
    end
    
    % Apply BC
    % Always fixed at ends, so we just take the inner part
    Kc = K(2:end-1,2:end-1);
    Fc = F(2:end-1);

    % Solve for U
    U = zeros(num_nodes,1);
    U(2:end-1) = Kc\Fc;
    
    % Finally, construct solution field. We use many points so we can 
    % accurately visualize the quadratic element solutions
    Ux = [];
    xx = [];
    for elem = 1:num_elems
        sctr = B(elem,:);
        xa = node_locs(sctr(1));
        xb = node_locs(sctr(end));
        x = linspace(xa, xb, 21);
        if (xb == 0.5 && xa == 0) || xb == 0.25
            x = [x 0.1666];
        end
        if (xb == 1.0 && xa == 0.5) || xb == 0.75
            x = [x 0.666];
        end
        x = sort(x);
        xbar = x - xa;
        if (elem_order == 1)
            phis = [1 - xbar/h; xbar/h]';
        else
            phis = [(1 - xbar/h).*(1 - (1/alpha)*xbar/h);
                    1/(alpha*(1-alpha))*xbar/h.*(1 - xbar/h);
                    -alpha/(1-alpha)*xbar/h.*(1 - 1/alpha*xbar/h)]';
        end
        uelem = phis*U(sctr);
        xx = [xx x];
        Ux = [Ux uelem'];
    end
    xx
    Ux
    dat = [xx' Ux'];
    csvwrite(names{i}, dat);
    dat = [node_locs' U];
    csvwrite(strcat(names{i}, 'x'), dat);
end
