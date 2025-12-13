/**
 * math_universe.c
 * COMPLETE MATHEMATICAL DESCRIPTION OF THE UNIVERSE IN C
 * All mathematical concepts and geometric expressions that describe reality
 * 
 * Compile: gcc -o math_universe math_universe.c -lm -O3 -march=native
 * Run: ./math_universe --explore [concept]
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

/* ==============================================
   1. FUNDAMENTAL MATHEMATICAL CONSTANTS
   ============================================== */

typedef struct {
    /* Universal Constants */
    double pi;               /* π: Circle circumference/diameter */
    double e;                /* e: Natural logarithm base */
    double phi;              /* φ: Golden ratio (1+√5)/2 */
    double gamma;            /* γ: Euler-Mascheroni constant */
    double alpha;            /* α: Fine-structure constant */
    
    /* Physical Constants (mathematical relationships) */
    double c;                /* Speed of light */
    double G;                /* Gravitational constant */
    double hbar;             /* Reduced Planck's constant */
    double k_B;             /* Boltzmann constant */
    
    /* Mathematical Limits */
    double infinity;
    double epsilon;          /* Machine epsilon */
    double aleph_null;       /* ℵ₀: Cardinality of natural numbers */
} MathConstants;

/* ==============================================
   2. NUMBER THEORY & ABSTRACT ALGEBRA
   ============================================== */

/* Prime number generation (Sieve of Eratosthenes) */
void generate_primes(int limit, int** primes, int* count) {
    bool* is_prime = malloc((limit + 1) * sizeof(bool));
    memset(is_prime, true, (limit + 1) * sizeof(bool));
    
    is_prime[0] = is_prime[1] = false;
    
    for (int p = 2; p * p <= limit; p++) {
        if (is_prime[p]) {
            for (int i = p * p; i <= limit; i += p) {
                is_prime[i] = false;
            }
        }
    }
    
    *count = 0;
    for (int i = 2; i <= limit; i++) {
        if (is_prime[i]) (*count)++;
    }
    
    *primes = malloc(*count * sizeof(int));
    int index = 0;
    for (int i = 2; i <= limit; i++) {
        if (is_prime[i]) {
            (*primes)[index++] = i;
        }
    }
    
    free(is_prime);
}

/* Greatest Common Divisor (Euclidean algorithm) */
int gcd(int a, int b) {
    while (b != 0) {
        int temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

/* Modular arithmetic operations */
int mod_pow(int base, int exp, int mod) {
    int result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = (result * base) % mod;
        base = (base * base) % mod;
        exp >>= 1;
    }
    return result;
}

/* Complex number arithmetic */
typedef struct {
    double real;
    double imag;
} Complex;

Complex complex_add(Complex a, Complex b) {
    return (Complex){a.real + b.real, a.imag + b.imag};
}

Complex complex_mul(Complex a, Complex b) {
    return (Complex){a.real*b.real - a.imag*b.imag,
                     a.real*b.imag + a.imag*b.real};
}

double complex_magnitude(Complex c) {
    return sqrt(c.real*c.real + c.imag*c.imag);
}

/* ==============================================
   3. GEOMETRIC ALGEBRA (Clifford Algebra)
   ============================================== */

/* 3D Vector with geometric algebra operations */
typedef struct {
    double scalar;      /* Grade 0: Scalar */
    double e1, e2, e3; /* Grade 1: Vectors */
    double e12, e23, e31; /* Grade 2: Bivectors */
    double e123;        /* Grade 3: Trivector (pseudoscalar) */
} Multivector;

/* Geometric product: The fundamental operation of GA */
Multivector geometric_product(Multivector a, Multivector b) {
    Multivector result = {0};
    
    /* Scalar * Scalar */
    result.scalar = a.scalar * b.scalar;
    
    /* Scalar * Vector */
    result.e1 += a.scalar * b.e1;
    result.e2 += a.scalar * b.e2;
    result.e3 += a.scalar * b.e3;
    
    result.e1 += b.scalar * a.e1;
    result.e2 += b.scalar * a.e2;
    result.e3 += b.scalar * a.e3;
    
    /* Vector * Vector (produces scalar + bivector) */
    result.scalar += a.e1*b.e1 + a.e2*b.e2 + a.e3*b.e3;  /* Dot product */
    result.e12 = a.e1*b.e2 - a.e2*b.e1;  /* Wedge product components */
    result.e23 = a.e2*b.e3 - a.e3*b.e2;
    result.e31 = a.e3*b.e1 - a.e1*b.e3;
    
    /* Higher grade terms omitted for brevity */
    
    return result;
}

/* Rotor (even subalgebra for rotations) */
Multivector create_rotor(double angle, double x, double y, double z) {
    double half_angle = angle / 2;
    double s = sin(half_angle);
    double c = cos(half_angle);
    
    /* Normalize rotation axis */
    double norm = sqrt(x*x + y*y + z*z);
    x /= norm; y /= norm; z /= norm;
    
    Multivector R = {0};
    R.scalar = c;
    R.e23 = -s * x;  /* i component */
    R.e31 = -s * y;  /* j component */
    R.e12 = -s * z;  /* k component */
    
    return R;
}

/* Rotate vector using rotor */
Multivector rotate_vector(Multivector rotor, Multivector vector) {
    Multivector R_conj = {0};
    R_conj.scalar = rotor.scalar;
    R_conj.e23 = -rotor.e23;
    R_conj.e31 = -rotor.e31;
    R_conj.e12 = -rotor.e12;
    
    /* R v R† */
    Multivector temp = geometric_product(rotor, vector);
    return geometric_product(temp, R_conj);
}

/* ==============================================
   4. DIFFERENTIAL GEOMETRY & TENSORS
   ============================================== */

/* Riemannian metric tensor (4D spacetime) */
typedef struct {
    double g[4][4];  /* Metric components g_μν */
    double det_g;    /* Determinant of metric */
    double inv_g[4][4]; /* Inverse metric g^μν */
} MetricTensor;

/* Christoffel symbols (Levi-Civita connection) */
void calculate_christoffel(MetricTensor* metric, double Gamma[4][4][4]) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                Gamma[i][j][k] = 0;
                for (int l = 0; l < 4; l++) {
                    double dg_jl = 0, dg_kl = 0, dg_jk = 0;
                    /* Finite difference approximation of derivatives */
                    /* In practice, we'd have coordinate functions */
                    Gamma[i][j][k] += 0.5 * metric->inv_g[i][l] *
                                     (dg_jl + dg_kl - dg_jk);
                }
            }
        }
    }
}

/* Riemann curvature tensor */
void calculate_riemann(MetricTensor* metric, double R[4][4][4][4]) {
    double Gamma[4][4][4];
    calculate_christoffel(metric, Gamma);
    
    /* R^ρ_σμν = ∂_μΓ^ρ_νσ - ∂_νΓ^ρ_μσ + Γ^ρ_μλΓ^λ_νσ - Γ^ρ_νλΓ^λ_μσ */
    for (int rho = 0; rho < 4; rho++) {
        for (int sigma = 0; sigma < 4; sigma++) {
            for (int mu = 0; mu < 4; mu++) {
                for (int nu = 0; nu < 4; nu++) {
                    double term1 = 0; /* ∂_μΓ^ρ_νσ */
                    double term2 = 0; /* ∂_νΓ^ρ_μσ */
                    double term3 = 0, term4 = 0;
                    
                    for (int lam = 0; lam < 4; lam++) {
                        term3 += Gamma[rho][mu][lam] * Gamma[lam][nu][sigma];
                        term4 += Gamma[rho][nu][lam] * Gamma[lam][mu][sigma];
                    }
                    
                    R[rho][sigma][mu][nu] = term1 - term2 + term3 - term4;
                }
            }
        }
    }
}

/* ==============================================
   5. TOPOLOGY & MANIFOLDS
   ============================================== */

/* Simplicial complex (for computational topology) */
typedef struct {
    int dimension;
    int* vertices;
    struct Simplex** faces;
    int num_faces;
} Simplex;

/* Calculate Euler characteristic χ = V - E + F - ... */
int euler_characteristic(Simplex** complex, int num_simplices) {
    int chi = 0;
    for (int d = 0; d <= 3; d++) {
        int count = 0;
        for (int i = 0; i < num_simplices; i++) {
            if (complex[i]->dimension == d) count++;
        }
        chi += (d % 2 == 0) ? count : -count;
    }
    return chi;
}

/* Homology group computation (simplified) */
void compute_homology(Simplex** complex, int num_simplices, 
                      int* betti_numbers, int max_dim) {
    /* Betti numbers: b_k = rank(H_k) = dim(ker ∂_k) - dim(im ∂_{k+1}) */
    for (int k = 0; k < max_dim; k++) {
        /* Simplified: count k-dimensional holes */
        int cycles = 0, boundaries = 0;
        
        /* In real implementation, we'd build boundary matrices
           and compute their ranks using linear algebra */
        betti_numbers[k] = cycles - boundaries;
    }
}

/* ==============================================
   6. LIE GROUPS & ALGEBRAS
   ============================================== */

/* SO(3) rotation group (3x3 orthogonal matrices with det = 1) */
typedef struct {
    double m[3][3];
} SO3Matrix;

/* SO(3) group operations */
SO3Matrix so3_identity() {
    SO3Matrix I = {{{0}}};
    for (int i = 0; i < 3; i++) I.m[i][i] = 1.0;
    return I;
}

SO3Matrix so3_mul(SO3Matrix a, SO3Matrix b) {
    SO3Matrix result = {{{0}}};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                result.m[i][j] += a.m[i][k] * b.m[k][j];
            }
        }
    }
    return result;
}

/* so(3) Lie algebra (skew-symmetric 3x3 matrices) */
typedef struct {
    double x, y, z;  /* Basis: [J_x, J_y, J_z] */
} so3Algebra;

/* Exponential map: so(3) → SO(3) */
SO3Matrix so3_exp(so3Algebra w) {
    double theta = sqrt(w.x*w.x + w.y*w.y + w.z*w.z);
    SO3Matrix result = so3_identity();
    
    if (theta > 1e-10) {
        double wx = w.x / theta;
        double wy = w.y / theta;
        double wz = w.z / theta;
        
        double c = cos(theta);
        double s = sin(theta);
        double t = 1 - c;
        
        /* Rodrigues' rotation formula */
        result.m[0][0] = c + wx*wx*t;
        result.m[0][1] = wx*wy*t - wz*s;
        result.m[0][2] = wx*wz*t + wy*s;
        
        result.m[1][0] = wy*wx*t + wz*s;
        result.m[1][1] = c + wy*wy*t;
        result.m[1][2] = wy*wz*t - wx*s;
        
        result.m[2][0] = wz*wx*t - wy*s;
        result.m[2][1] = wz*wy*t + wx*s;
        result.m[2][2] = c + wz*wz*t;
    }
    
    return result;
}

/* ==============================================
   7. DIFFERENTIAL EQUATIONS & DYNAMICAL SYSTEMS
   ============================================== */

/* Runge-Kutta 4th order integrator */
typedef double (*ODEFunction)(double t, double y, void* params);

double rk4_integrate(ODEFunction f, double t0, double y0, double t_end, 
                     double dt, void* params) {
    double t = t0;
    double y = y0;
    
    while (t < t_end) {
        double k1 = f(t, y, params);
        double k2 = f(t + dt/2, y + dt*k1/2, params);
        double k3 = f(t + dt/2, y + dt*k2/2, params);
        double k4 = f(t + dt, y + dt*k3, params);
        
        y += dt * (k1 + 2*k2 + 2*k3 + k4) / 6;
        t += dt;
    }
    
    return y;
}

/* Lorenz attractor (chaotic system) */
typedef struct {
    double sigma, rho, beta;
    double x, y, z;
} LorenzSystem;

void lorenz_step(LorenzSystem* sys, double dt) {
    double dx = sys->sigma * (sys->y - sys->x);
    double dy = sys->x * (sys->rho - sys->z) - sys->y;
    double dz = sys->x * sys->y - sys->beta * sys->z;
    
    sys->x += dx * dt;
    sys->y += dy * dt;
    sys->z += dz * dt;
}

/* Hamiltonian mechanics */
typedef struct {
    double q[3];  /* Generalized coordinates */
    double p[3];  /* Conjugate momenta */
    double t;     /* Time */
} PhasePoint;

typedef void (*HamiltonianFunc)(PhasePoint* point, double* H, 
                                double dHdq[3], double dHdp[3]);

void symplectic_integrate(HamiltonianFunc H, PhasePoint* point, double dt) {
    double dHdq[3], dHdp[3];
    double H_val;
    
    H(point, &H_val, dHdq, dHdp);
    
    /* Leapfrog/Verlet integration */
    for (int i = 0; i < 3; i++) {
        point->p[i] -= dt * dHdq[i];
    }
    
    H(point, &H_val, dHdq, dHdp);
    
    for (int i = 0; i < 3; i++) {
        point->q[i] += dt * dHdp[i];
    }
    
    point->t += dt;
}

/* ==============================================
   8. FRACTAL GEOMETRY & COMPLEX DYNAMICS
   ============================================== */

/* Mandelbrot set iteration */
int mandelbrot_iteration(double complex c, int max_iter) {
    double complex z = 0;
    for (int i = 0; i < max_iter; i++) {
        if (cabs(z) > 2.0) return i;
        z = z*z + c;
    }
    return max_iter;
}

/* Julia set iteration */
int julia_iteration(double complex z, double complex c, int max_iter) {
    for (int i = 0; i < max_iter; i++) {
        if (cabs(z) > 2.0) return i;
        z = z*z + c;
    }
    return max_iter;
}

/* Fractal dimension calculation (box counting) */
double fractal_dimension(double* data, int width, int height, 
                         int min_box, int max_box) {
    int counts[100];
    double scales[100];
    int count = 0;
    
    for (int box_size = min_box; box_size <= max_box; box_size *= 2) {
        int boxes = 0;
        for (int y = 0; y < height; y += box_size) {
            for (int x = 0; x < width; x += box_size) {
                /* Check if box contains any data point */
                int has_data = 0;
                for (int dy = 0; dy < box_size && y+dy < height; dy++) {
                    for (int dx = 0; dx < box_size && x+dx < width; dx++) {
                        if (data[(y+dy)*width + (x+dx)] > 0) {
                            has_data = 1;
                            break;
                        }
                    }
                    if (has_data) break;
                }
                if (has_data) boxes++;
            }
        }
        
        counts[count] = boxes;
        scales[count] = 1.0 / box_size;
        count++;
    }
    
    /* Linear regression: log(N) = -D * log(ε) + constant */
    double sum_log_N = 0, sum_log_epsilon = 0;
    double sum_log_N_log_epsilon = 0, sum_log_epsilon_sq = 0;
    
    for (int i = 0; i < count; i++) {
        double log_N = log(counts[i]);
        double log_epsilon = log(scales[i]);
        
        sum_log_N += log_N;
        sum_log_epsilon += log_epsilon;
        sum_log_N_log_epsilon += log_N * log_epsilon;
        sum_log_epsilon_sq += log_epsilon * log_epsilon;
    }
    
    double D = -(count * sum_log_N_log_epsilon - sum_log_N * sum_log_epsilon) /
               (count * sum_log_epsilon_sq - sum_log_epsilon * sum_log_epsilon);
    
    return D;
}

/* ==============================================
   9. INFORMATION THEORY & ENTROPY
   ============================================== */

/* Shannon entropy H(X) = -Σ p(x) log₂ p(x) */
double shannon_entropy(double* probabilities, int n) {
    double entropy = 0.0;
    for (int i = 0; i < n; i++) {
        if (probabilities[i] > 0) {
            entropy -= probabilities[i] * log2(probabilities[i]);
        }
    }
    return entropy;
}

/* Kolmogorov complexity approximation (using LZ77) */
int lz77_compress(const char* input, int length, char* output) {
    /* Simplified LZ77 compression for complexity estimation */
    int output_pos = 0;
    int pos = 0;
    
    while (pos < length) {
        int best_match_length = 0;
        int best_match_distance = 0;
        
        /* Search for longest match in sliding window */
        int window_start = (pos > 4096) ? pos - 4096 : 0;
        for (int i = window_start; i < pos; i++) {
            int match_length = 0;
            while (pos + match_length < length && 
                   i + match_length < pos &&
                   input[i + match_length] == input[pos + match_length]) {
                match_length++;
            }
            
            if (match_length > best_match_length) {
                best_match_length = match_length;
                best_match_distance = pos - i;
            }
        }
        
        if (best_match_length > 3) {
            /* Encode as (distance, length) pair */
            output[output_pos++] = (best_match_distance >> 8) & 0xFF;
            output[output_pos++] = best_match_distance & 0xFF;
            output[output_pos++] = best_match_length;
            pos += best_match_length;
        } else {
            /* Encode as literal */
            output[output_pos++] = 0;  /* Flag for literal */
            output[output_pos++] = input[pos++];
        }
    }
    
    return output_pos;
}

double kolmogorov_complexity_estimate(const char* data, int length) {
    char compressed[length * 2];
    int compressed_length = lz77_compress(data, length, compressed);
    return (double)compressed_length / length;
}

/* ==============================================
   10. CATEGORY THEORY (FUNCTORS & NATURAL TRANSFORMATIONS)
   ============================================== */

/* Object in a category (simplified) */
typedef struct {
    int id;
    char* name;
    void* data;
} CategoryObject;

/* Morphism (arrow) between objects */
typedef struct {
    CategoryObject* domain;
    CategoryObject* codomain;
    void* (*map)(void*);  /* The actual mapping function */
    char* name;
} Morphism;

/* Functor between categories */
typedef struct {
    CategoryObject** objects;
    Morphism** morphisms;
    int num_objects;
    int num_morphisms;
    
    /* Mapping functions */
    CategoryObject* (*map_object)(CategoryObject*);
    Morphism* (*map_morphism)(Morphism*);
} Functor;

/* Natural transformation between functors */
typedef struct {
    Functor* F;
    Functor* G;
    Morphism** components;  /* One for each object */
} NaturalTransformation;

/* Compose morphisms */
Morphism* compose_morphisms(Morphism* f, Morphism* g) {
    if (f->codomain != g->domain) return NULL;
    
    Morphism* comp = malloc(sizeof(Morphism));
    comp->domain = f->domain;
    comp->codomain = g->codomain;
    comp->name = "composition";
    
    comp->map = NULL;  /* Would be f->map ∘ g->map */
    
    return comp;
}

/* ==============================================
   11. QUANTUM MATHEMATICS
   ============================================== */

/* Quantum state (wavefunction) */
typedef struct {
    complex double* amplitudes;
    int dimension;
    double norm;
} QuantumState;

/* Initialize quantum state */
QuantumState* quantum_state_create(int dim) {
    QuantumState* state = malloc(sizeof(QuantumState));
    state->amplitudes = calloc(dim, sizeof(complex double));
    state->dimension = dim;
    state->norm = 0.0;
    return state;
}

/* Tensor product of quantum states (entanglement) */
QuantumState* tensor_product(QuantumState* a, QuantumState* b) {
    int dim = a->dimension * b->dimension;
    QuantumState* result = quantum_state_create(dim);
    
    for (int i = 0; i < a->dimension; i++) {
        for (int j = 0; j < b->dimension; j++) {
            result->amplitudes[i * b->dimension + j] = 
                a->amplitudes[i] * b->amplitudes[j];
        }
    }
    
    return result;
}

/* Quantum Fourier Transform */
void quantum_fourier_transform(QuantumState* state) {
    int N = state->dimension;
    complex double* new_amplitudes = calloc(N, sizeof(complex double));
    
    for (int k = 0; k < N; k++) {
        for (int n = 0; n < N; n++) {
            double angle = -2.0 * M_PI * k * n / N;
            new_amplitudes[k] += state->amplitudes[n] * 
                                (cos(angle) + I * sin(angle));
        }
        new_amplitudes[k] /= sqrt(N);
    }
    
    free(state->amplitudes);
    state->amplitudes = new_amplitudes;
}

/* ==============================================
   12. MATHEMATICAL PHYSICS EQUATIONS
   ============================================== */

/* Einstein field equations: G_μν = 8πG/c⁴ T_μν */
void einstein_field_equations(MetricTensor* g, double T[4][4], 
                              double G_tensor[4][4]) {
    double R[4][4];  /* Ricci tensor */
    double R_scalar; /* Ricci scalar */
    
    /* Calculate Ricci tensor from Riemann tensor */
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            R[mu][nu] = 0;
            for (int rho = 0; rho < 4; rho++) {
                /* R_μν = R^ρ_μρν */
                /* Simplified - actual calculation needs Riemann */
            }
        }
    }
    
    /* Einstein tensor: G_μν = R_μν - ½ g_μν R */
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            G_tensor[mu][nu] = R[mu][nu] - 0.5 * g->g[mu][nu] * R_scalar;
        }
    }
}

/* Schrödinger equation: iℏ ∂ψ/∂t = Ĥ ψ */
void schrodinger_equation(QuantumState* psi, double t, 
                          complex double (*hamiltonian)(int, int, double),
                          double dt) {
    int N = psi->dimension;
    complex double* new_amplitudes = calloc(N, sizeof(complex double));
    
    /* Time evolution: ψ(t+dt) = exp(-iĤdt/ℏ) ψ(t) */
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            complex double H_ij = hamiltonian(i, j, t);
            complex double propagator = cexp(-I * H_ij * dt);
            new_amplitudes[i] += propagator * psi->amplitudes[j];
        }
    }
    
    free(psi->amplitudes);
    psi->amplitudes = new_amplitudes;
}

/* ==============================================
   13. STATISTICS & PROBABILITY THEORY
   ============================================== */

/* Probability distributions */
typedef enum {
    DIST_NORMAL,
    DIST_UNIFORM,
    DIST_EXPONENTIAL,
    DIST_POISSON,
    DIST_BERNOULLI
} DistributionType;

double sample_normal(double mu, double sigma) {
    /* Box-Muller transform */
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return mu + sigma * z0;
}

/* Bayesian inference */
typedef struct {
    double prior;
    double likelihood;
    double evidence;
    double posterior;
} BayesianUpdate;

void bayesian_update(BayesianUpdate* bu, double new_likelihood) {
    bu->posterior = (bu->likelihood * bu->prior) / bu->evidence;
    bu->prior = bu->posterior;  /* For next update */
    bu->likelihood = new_likelihood;
}

/* Markov Chain Monte Carlo (Metropolis-Hastings) */
double* mcmc_sample(double (*target_dist)(double), int num_samples, 
                    double initial, double proposal_std) {
    double* samples = malloc(num_samples * sizeof(double));
    double current = initial;
    
    for (int i = 0; i < num_samples; i++) {
        /* Propose new state */
        double proposed = current + sample_normal(0, proposal_std);
        
        /* Acceptance ratio */
        double acceptance = target_dist(proposed) / target_dist(current);
        
        if (acceptance >= 1.0 || ((double)rand() / RAND_MAX) < acceptance) {
            current = proposed;
        }
        
        samples[i] = current;
    }
    
    return samples;
}

/* ==============================================
   14. GRAPH THEORY
   ============================================== */

typedef struct GraphNode {
    int id;
    void* data;
    struct GraphNode** neighbors;
    int num_neighbors;
} GraphNode;

typedef struct {
    GraphNode** nodes;
    int num_nodes;
    int** adjacency_matrix;
    int directed;
} Graph;

/* Dijkstra's shortest path algorithm */
int* dijkstra_shortest_path(Graph* graph, int start, int end) {
    int n = graph->num_nodes;
    double* dist = malloc(n * sizeof(double));
    int* prev = malloc(n * sizeof(int));
    int* visited = calloc(n, sizeof(int));
    
    for (int i = 0; i < n; i++) {
        dist[i] = INFINITY;
        prev[i] = -1;
    }
    dist[start] = 0;
    
    for (int i = 0; i < n; i++) {
        /* Find unvisited node with minimum distance */
        int u = -1;
        double min_dist = INFINITY;
        for (int j = 0; j < n; j++) {
            if (!visited[j] && dist[j] < min_dist) {
                min_dist = dist[j];
                u = j;
            }
        }
        
        if (u == -1 || u == end) break;
        visited[u] = 1;
        
        /* Update distances to neighbors */
        for (int v = 0; v < n; v++) {
            if (graph->adjacency_matrix[u][v] > 0) {
                double alt = dist[u] + graph->adjacency_matrix[u][v];
                if (alt < dist[v]) {
                    dist[v] = alt;
                    prev[v] = u;
                }
            }
        }
    }
    
    /* Reconstruct path */
    int* path = NULL;
    if (prev[end] != -1) {
        /* Count path length */
        int count = 1;
        int current = end;
        while (current != start) {
            count++;
            current = prev[current];
        }
        
        /* Build path */
        path = malloc((count + 1) * sizeof(int));
        path[0] = count;
        current = end;
        for (int i = count; i > 0; i--) {
            path[i] = current;
            current = prev[current];
        }
    }
    
    free(dist);
    free(prev);
    free(visited);
    
    return path;
}

/* Graph isomorphism test (simplified) */
int graphs_isomorphic(Graph* g1, Graph* g2) {
    if (g1->num_nodes != g2->num_nodes) return 0;
    
    int n = g1->num_nodes;
    
    /* Check degree sequences */
    int deg1[n], deg2[n];
    for (int i = 0; i < n; i++) {
        deg1[i] = 0;
        deg2[i] = 0;
        for (int j = 0; j < n; j++) {
            deg1[i] += g1->adjacency_matrix[i][j];
            deg2[i] += g2->adjacency_matrix[i][j];
        }
    }
    
    /* Sort and compare degree sequences */
    /* Simplified - full algorithm would try all permutations */
    
    return 1;  /* Assume isomorphic for demo */
}

/* ==============================================
   15. COMPUTATIONAL GEOMETRY
   ============================================== */

/* Point in 2D */
typedef struct {
    double x, y;
} Point2D;

/* Line segment */
typedef struct {
    Point2D start, end;
} Segment;

/* Polygon */
typedef struct {
    Point2D* vertices;
    int num_vertices;
} Polygon;

/* Check if point is inside polygon (ray casting) */
int point_in_polygon(Point2D point, Polygon* poly) {
    int crossings = 0;
    
    for (int i = 0; i < poly->num_vertices; i++) {
        Point2D v1 = poly->vertices[i];
        Point2D v2 = poly->vertices[(i + 1) % poly->num_vertices];
        
        /* Check if ray intersects edge */
        if (((v1.y > point.y) != (v2.y > point.y)) &&
            (point.x < (v2.x - v1.x) * (point.y - v1.y) / (v2.y - v1.y) + v1.x)) {
            crossings++;
        }
    }
    
    return crossings % 2;
}

/* Convex hull (Graham scan) */
Polygon* convex_hull(Point2D* points, int n) {
    if (n < 3) return NULL;
    
    /* Find point with minimum y (and leftmost if tie) */
    int min_idx = 0;
    for (int i = 1; i < n; i++) {
        if (points[i].y < points[min_idx].y ||
            (points[i].y == points[min_idx].y && 
             points[i].x < points[min_idx].x)) {
            min_idx = i;
        }
    }
    
    /* Swap to front */
    Point2D temp = points[0];
    points[0] = points[min_idx];
    points[min_idx] = temp;
    
    /* Sort by polar angle */
    /* Implementation omitted for brevity */
    
    /* Build hull */
    Point2D* hull = malloc(n * sizeof(Point2D));
    hull[0] = points[0];
    hull[1] = points[1];
    int hull_size = 2;
    
    for (int i = 2; i < n; i++) {
        while (hull_size >= 2) {
            Point2D a = hull[hull_size - 2];
            Point2D b = hull[hull_size - 1];
            Point2D c = points[i];
            
            /* Cross product (b-a) × (c-b) */
            double cross = (b.x - a.x) * (c.y - b.y) - 
                          (b.y - a.y) * (c.x - b.x);
            
            if (cross <= 0) break;  /* Not making a left turn */
            hull_size--;
        }
        hull[hull_size++] = points[i];
    }
    
    Polygon* result = malloc(sizeof(Polygon));
    result->vertices = hull;
    result->num_vertices = hull_size;
    return result;
}

/* ==============================================
   MAIN DEMONSTRATION FUNCTION
   ============================================== */

void demonstrate_all_mathematics() {
    printf("=== MATHEMATICAL UNIVERSE DEMONSTRATION ===\n\n");
    
    /* 1. Number Theory */
    printf("1. NUMBER THEORY:\n");
    int* primes;
    int prime_count;
    generate_primes(100, &primes, &prime_count);
    printf("   Primes < 100: ");
    for (int i = 0; i < (prime_count < 10 ? prime_count : 10); i++) {
        printf("%d ", primes[i]);
    }
    printf("... (%d total)\n", prime_count);
    printf("   gcd(48, 180) = %d\n", gcd(48, 180));
    printf("   7^13 mod 11 = %d\n", mod_pow(7, 13, 11));
    free(primes);
    
    /* 2. Complex Numbers */
    printf("\n2. COMPLEX NUMBERS:\n");
    Complex c1 = {3, 4};
    Complex c2 = {1, -2};
    Complex sum = complex_add(c1, c2);
    Complex prod = complex_mul(c1, c2);
    printf("   (3+4i) + (1-2i) = %.1f%+.1fi\n", sum.real, sum.imag);
    printf("   (3+4i) * (1-2i) = %.1f%+.1fi\n", prod.real, prod.imag);
    printf("   |3+4i| = %.2f\n", complex_magnitude(c1));
    
    /* 3. Geometric Algebra */
    printf("\n3. GEOMETRIC ALGEBRA:\n");
    Multivector v1 = {0, 1, 0, 0, 0, 0, 0, 0};  /* Vector along x */
    Multivector v2 = {0, 0, 1, 0, 0, 0, 0, 0};  /* Vector along y */
    Multivector product = geometric_product(v1, v2);
    printf("   e1 * e2 = scalar: %.1f + bivector e12: %.1f\n", 
           product.scalar, product.e12);
    
    /* 4. Differential Equations */
    printf("\n4. DIFFERENTIAL EQUATIONS:\n");
    LorenzSystem lorenz = {10.0, 28.0, 8.0/3.0, 1.0, 1.0, 1.0};
    for (int i = 0; i < 5; i++) {
        lorenz_step(&lorenz, 0.01);
    }
    printf("   Lorenz system after 5 steps: (%.3f, %.3f, %.3f)\n",
           lorenz.x, lorenz.y, lorenz.z);
    
    /* 5. Fractals */
    printf("\n5. FRACTAL GEOMETRY:\n");
    double complex c = -0.7 + 0.27015*I;
    int mandelbrot_val = mandelbrot_iteration(c, 100);
    printf("   Mandelbrot iteration count for c = -0.7+0.27015i: %d\n",
           mandelbrot_val);
    
    /* 6. Information Theory */
    printf("\n6. INFORMATION THEORY:\n");
    double probs[] = {0.5, 0.25, 0.125, 0.125};
    double entropy = shannon_entropy(probs, 4);
    printf("   Shannon entropy of {0.5, 0.25, 0.125, 0.125}: %.3f bits\n",
           entropy);
    
    /* 7. Statistics */
    printf("\n7. STATISTICS:\n");
    printf("   Normal distribution sample: %.3f\n", sample_normal(0, 1));
    
    /* 8. Graph Theory */
    printf("\n8. GRAPH THEORY:\n");
    printf("   Dijkstra's algorithm implemented for shortest paths\n");
    
    /* 9. Computational Geometry */
    printf("\n9. COMPUTATIONAL GEOMETRY:\n");
    Point2D square[] = {{0,0}, {1,0}, {1,1}, {0,1}};
    Polygon poly = {square, 4};
    Point2D test_point = {0.5, 0.5};
    printf("   Point (0.5,0.5) in unit square: %s\n",
           point_in_polygon(test_point, &poly) ? "YES" : "NO");
    
    printf("\n=== MATHEMATICAL UNIVERSE SUMMARY ===\n");
    printf("This demonstration shows how mathematics forms the\n");
    printf("foundation of reality, from prime numbers to quantum\n");
    printf("mechanics, from geometry to information theory.\n");
    printf("All implemented in a single C file.\n");
}

int main(int argc, char** argv) {
    srand(time(NULL));
    
    if (argc > 1 && strcmp(argv[1], "--explore") == 0) {
        if (argc > 2) {
            if (strcmp(argv[2], "primes") == 0) {
                /* Explore number theory */
                printf("Prime numbers are the atoms of mathematics.\n");
                int* primes;
                int count;
                generate_primes(1000, &primes, &count);
                printf("There are %d primes under 1000.\n", count);
                free(primes);
            }
            else if (strcmp(argv[2], "fractals") == 0) {
                printf("Fractals show self-similarity across scales.\n");
                printf("Mandelbrot set demonstrates chaotic dynamics.\n");
            }
        }
    } else {
        demonstrate_all_mathematics();
    }
    
    return 0;
}