graphics_toolkit("gnuplot");

a = 0; b = 1;
k = @(x) x.^2 + 1;
q = @(x) 10*(x + 1);
f = @(x) 10*x.^5 - 30*x.^4 + 14*x.^3 + 22*x.^2 + 36*x - 2;
u_exact = @(x) x.^4 - 2*x.^3 + x.^2 + 3*x;
gamma = 3;
delta = 3;

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

    k_minus_half = k(a - h/2);
    k_plus_half = k(a + h/2);
    b_diag(1) = (k_minus_half + k_plus_half)/h^2 + qv(1);
    c_diag(1) = -(k_minus_half + k_plus_half)/h^2;
    d_vec(1) = fv(1) + 2*k_minus_half/h * gamma;

    k_minus_half = k(b - h/2);
    k_plus_half = k(b + h/2);
    a_diag(N-1) = -(k_minus_half + k_plus_half)/h^2;
    b_diag(N) = (k_minus_half + k_plus_half)/h^2 + qv(N);
    d_vec(N) = fv(N) + 2*k_plus_half/h * delta;

    aa = a_diag(1:N-1);
    cc = c_diag(1:N-1);
    n = length(d_vec);
    bb = b_diag;
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
            x_final2 = x;
            U_final2 = U;
            break;
        end
    else
        first = false;
    end

    U_prev = U;
    N = 2*N - 1;
end

x_final2 = x;
U_final2 = U;
fprintf('Runge_error_estimate: %.6e\n', err_est);
u_exact_vals = u_exact(x_final2);
err_fact2 = max(abs(U_final2 - u_exact_vals));
fprintf('Actual_error: %.6e\n', err_fact2);

figure;
plot(x_final2, U_final2, 'ro', x_final2, u_exact_vals, 'b-');
xlabel('x'); ylabel('u(x)');
legend('Numerical', 'Exact', 'Location', 'northwest');
grid on;

figure;
plot(x_final2, abs(U_final2 - u_exact_vals), 'k-o');
xlabel('x'); ylabel('Error');
grid on;
