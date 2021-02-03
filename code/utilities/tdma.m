function x = tdma(A, b)

% Data and initial definitions
N = size(A,1);  % Linear system dimension
P = zeros(N,1); % Vector of Ps
Q = P;          % Vector of Qs

% Compute Ps and Qs which is better to not compute algorithmically
P(1) = A(1,3)/A(1,2);
Q(1) = b(1)/A(1,2);
P(N) = 0;
% Compute algorithmically (for)
for i = 2:N-1
    P(i) = A(i,3)/(A(i,2) - A(i,1)*P(i-1));
    Q(i) = (b(i) + A(i,1)*Q(i-1))/(A(i,2) - A(i,1)*P(i-1));
end
% Compute Q(N) not algorithmically again because I'm stupid
Q(N) = (b(N) + A(N,1)*Q(N-1))/(A(N,2) - A(N,1)*P(N-1));

% Compute solution
x = zeros(N,1);
x(N) = Q(N);
for i = N-1:-1:1
    x(i) = P(i)*x(i+1) + Q(i);
end


end