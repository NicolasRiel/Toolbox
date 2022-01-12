clear; set(0,'defaultfigurecolor',[1 1 1]);
% 19-03-2021 N. Riel
% generalized simplex discretization in n-Dimensions with lower and upper bounds

%% Define the parameters %%
D           = 4;                        % number of dimensions
ABS         = 100;                      % max number of point
step        = 1/ABS;                    % compositional step
istep(1:D)  = 20;

% istep(1:D)  = [8,4,8,16];         % lower bound

MAXB(1:D)   = [60,60,20,60];    % upper bound
MINB(1:D)   = [0,0,0,10];         % lower bound

% MINB(1:D) = [0, 0, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0];
% MAXB(1:D) = [3, 1, 1, 1, 1, 5, 1, 1, 1, 7, 1, 3];

% MAXB(1:D)   = ABS;                      % upper bound
% MINB(1:D)   = 0;                        % lower bound

IP(1:D)     = 0;                        % indexes
ID(1:D)     = 0;                        % point coordinates
LB(1:D)     = 0;                        % active lower bound
UB(1:D)     = 0;                        % active upper bound

v           = [];                       % vector of compositions
w           = [];                       % vector of coordinates

%% Run the generalized simplex algorithm

UB = update_ub(D,ABS,MAXB,UB,IP);
LB = update_lb(D,MINB,LB);
IP = get_ip(D,MINB,LB);

p = D;
while IP(1) == 0
    
    ID(D) = ABS;
    for k=1:D-1
        ID(k) = IP(k+1);
        ID(D) = ID(D) - IP(k+1);
    end
    
    % save the coordinates and composition vectors
    if (ID(D) >= MINB(D) && ID(D) <= MAXB(D) && mod(ID(D),istep(1)) == 0)
        w = [w; ID];    
        v = [v; step*ID];
    end

    IP(p) = IP(p) + istep(p);
    while IP(p) > UB(p)
        IP(p)   = LB(p);
        p       = p - 1;
        IP(p)   = IP(p) + istep(p);
        UB      = update_ub(D,ABS,MAXB,UB,IP);

        if IP(p) <= UB(p)
            p = D;
        end
        
    end
end


% check results here %
info = check_sum(v);
disp(info)
if info == 1
   disp('Calculation is successful'); 
end

%% Plot results
figure(1)

subplot(1,3,1)
scatter3(v(:,1),v(:,2),v(:,3),'MarkerFaceColor','b')
daspect([1 1 1])
xlabel('1') 
ylabel('2')
zlabel('3')

subplot(1,3,2)
scatter3(v(:,2),v(:,3),v(:,4),'MarkerFaceColor','b')
daspect([1 1 1])
xlabel('2') 
ylabel('3')
zlabel('4')

subplot(1,3,3)
scatter3(v(:,3),v(:,4),v(:,1),'MarkerFaceColor','b')
daspect([1 1 1])
xlabel('3') 
ylabel('4')
zlabel('1')





%% Function list
function info = check_sum(v)
    R = sum(v,2);
    minSum = min(R);
    maxSum = min(R);
    
    if abs(maxSum - minSum) < 1e-8
        info = 1;
    else
        info = 0;
    end
end


function UB = update_ub(D,ABS,MAXB,UB,IP)
    UB(1) = 1e10;
    UB(2) = MAXB(1);
    for k=3:D
        UB(k) = ABS;
        for l=2:(k-1)
           UB(k) = UB(k) - IP(l);
        end
        UB(k) = min(UB(k),MAXB(k-1));
    end
end

function LB = update_lb(D,MINB,LB)
    for k=2:D
        LB(k) = MINB(k-1);
    end
end

function IP = get_ip(D,MINB,IP)
    for k=2:D
        IP(k) = MINB(k-1);
    end
end
