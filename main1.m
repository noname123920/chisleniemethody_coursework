graphics_toolkit("gnuplot");

a = 0; b = 1;
k = @(x) x.^2 + 1;
q = @(x) x + 1;
f = @(x) x.^5 - 21*x.^4 + 23*x.^3 - 14*x.^2 + 9*x - 2;
u_exact = @(x) x.^4 - 2*x.^3 + x.^2 + 3*x;

tol = 1e-4;
N = 5;
err_est = tol + 1;
first = true;

while err_est > tol
    h = (b - a)/(N - 1);
    x = linspace(a, b, N)';
    xmid = (x(1:end-1) + x(2:end))/2;
    kmid = k(xmid);
    qv = q(x);
    fv = f(x);

    a_diag = zeros(N,1);
    b_diag = zeros(N,1);
    c_diag = zeros(N,1);
    d_vec = zeros(N,1);

    for i = 2:N-1
        a_diag(i-1) = -kmid(i-1)/h^2;
        c_diag(i) = -kmid(i)/h^2;
        b_diag(i) = (kmid(i-1) + kmid(i))/h^2 + qv(i);
        d_vec(i) = fv(i);
    end

    b_diag(1) = 1; c_diag(1) = 0; d_vec(1) = u_exact(0);
    a_diag(N-1) = 0; b_diag(N) = 1; d_vec(N) = u_exact(1);

    n = length(d_vec);
    aa = a_diag(1:n-1);
    bb = b_diag;
    cc = c_diag(1:n-1);
    dd = d_vec;

    for i = 2:n
        w = aa(i-1) / bb(i-1);
        bb(i) = bb(i) - w * cc(i-1);
        dd(i) = dd(i) - w * dd(i-1);
    end

    U = zeros(n,1);
    U(n) = dd(n) / bb(n);
    for i = n-1:-1:1
        U(i) = (dd(i) - cc(i) * U(i+1)) / bb(i);
    end

    if ~first
        U_coarse = U(1:2:end);
        err_est = max(abs(U_coarse - U_prev)) / 3;
        fprintf('N=%d, h=%.6f, err_est=%.6e\n', N, h, err_est);
        if err_est < tol
            x_final1 = x;
            U_final1 = U;
            break;
        end
    else
        first = false;
    end

    U_prev = U;
    N = 2*N - 1;
end

x_final1 = x;
U_final1 = U;
u_exact_vals = u_exact(x_final1);
err_fact = max(abs(U_final1 - u_exact_vals));
fprintf('Actual_error: %.6e\n', err_fact);
fprintf('Runge_error_estimate: %.6e\n', err_est);

figure;
plot(x_final1, U_final1, 'ro', x_final1, u_exact_vals, 'b-');
xlabel('x'); ylabel('u(x)');
legend('Numerical', 'Exact', 'Location', 'northwest');
grid on;

figure;
plot(x_final1, abs(U_final1 - u_exact_vals), 'k-o');
xlabel('x'); ylabel('Error');
grid on;

