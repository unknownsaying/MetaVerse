/**
 * parallel_world.c
 * SIMULATION OF A PHYSICS-DEVIANT PARALLEL UNIVERSE
 * WHERE DIMENSIONAL CONSTANTS ARE UNSTABLE AND CAUSALITY IS NON-LINEAR
 * 
 * Compile: gcc -o parallel_world parallel_world.c -lm -fopenmp -O3 -march=native
 * Run: ./parallel_world --seed 666 --dimensions 11 --render --chaos 0.99
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <complex.h>
#include <stdatomic.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <immintrin.h>  // AVX/SSE intrinsics
#include <stdalign.h>

#pragma GCC optimize("Ofast,unroll-loops,fast-math")

/* 
   PART 1: PARALLEL PHYSICS CONSTANTS (UNSTABLE)
    */

/* Our universe's constants (for reference) */
#define OUR_PI 3.1415926
#define OUR_C 299792458.0
#define OUR_HBAR 1.054571817e-34
#define OUR_G 6.67430e-11

/* Parallel world's UNSTABLE constants - they FLUCTUATE based on position/time */
typedef struct {
    atomic_double pi;      // π varies spatially
    atomic_double c;       // Lightspeed isn't constant
    atomic_double G;       // Gravity strength fluctuates
    atomic_double hbar;    // Quantum uncertainty variable
    atomic_double epsilon; // Permittivity of NOTHING
    atomic_double alpha;   // Fine-structure parameter (changes charge)
    atomic_double tau;     // Time dilation factor
} ParallelConstants;

/* 
   PART 2: 11-DIMENSIONAL SPACETIME FABRIC
    */

/* We live in 4D (3 space + 1 time), parallel world has 7 extra folded dimensions */
typedef struct {
    alignas(64) double coord[11];  /* [x,y,z,t,u,v,w,a,b,c,d] */
    double metric[11][11];         /* Dynamic spacetime metric tensor */
    double curvature[11][11][11];  /* Riemann curvature tensor (simplified) */
    double torsion[11][11];        /* Torsion field (Einstein-Cartan theory) */
    int dimension_active[11];      /* Which dimensions are "unfolded" here */
} Hyperpoint;

/* Calabi-Yau manifold parameters (string theory extra dimensions) */
typedef struct {
    double complex moduli[6];      /* Shape parameters of folded dimensions */
    double kahler_potential;       /* Kähler potential */
    double superpotential;         /* Superpotential for SUSY breaking */
    int euler_characteristic;      /* Topological invariant */
} CalabiYauManifold;

/* 
   PART 3: NON-EUCLIDEAN GEOMETRY ENGINE
    */

/* Riemannian metric that violates Euclidean postulates */
double parallel_metric(Hyperpoint* p, int i, int j) {
    /* Metric depends on position in non-smooth way */
    double base = (i == j) ? 1.0 : 0.0;
    
    /* Add quantum foam fluctuations */
    double foam = 0.1 * sin(p->coord[3] * 100.0) * 
                  exp(-fabs(p->coord[0] * p->coord[1] * p->coord[2]));
    
    /* Add torsion contribution */
    double torsion_effect = p->torsion[i][j] * 0.01;
    
    /* Random metric perturbation (spacetime is "fuzzy") */
    double random_perturb = 0.01 * ((double)rand() / RAND_MAX - 0.5);
    
    return base + foam + torsion_effect + random_perturb;
}

/* Distance function that violates triangle inequality */
double parallel_distance(Hyperpoint* a, Hyperpoint* b) {
    double sum = 0.0;
    
    #pragma omp simd reduction(+:sum)
    for (int i = 0; i < 11; i++) {
        double diff = b->coord[i] - a->coord[i];
        /* Non-quadratic distance measure (p-norm with p varying) */
        double p = 2.0 + 0.5 * sin(a->coord[3] + i);  /* p varies with time and dimension! */
        sum += pow(fabs(diff), p);
    }
    
    /* Sometimes distance is NEGATIVE (non-metrizable space) */
    if (sin(a->coord[3] * 0.1) > 0.7) {
        return -pow(sum, 1.0/(2.0 + sin(a->coord[3])));
    }
    
    return pow(sum, 1.0/(2.0 + sin(a->coord[3])));
}

/* 
   PART 4: CAUSALITY VIOLATION ENGINE
    */

/* Closed timelike curves (time travel allowed) */
typedef struct {
    Hyperpoint entrance;
    Hyperpoint exit;
    double time_difference;  /* Can be NEGATIVE (exit before entrance) */
    double stability;        /* Probability curve will exist at given time */
    atomic_int use_count;    /* How many times traversed */
} CTCurve;

/* Event with non-linear causality */
typedef struct {
    Hyperpoint location;
    double proper_time;
    int64_t causal_id;
    int64_t* effects;        /* Events this CAUSES (forward in time) */
    int64_t* causes;         /* Events that CAUSED this (backward in time) */
    int effect_count;
    int cause_count;
    double probability;      /* Quantum amplitude */
    int happened;            /* Has it occurred? (Quantum superposition) */
} CausalEvent;

/* Causal violation detector */
int check_causal_violation(CausalEvent* event) {
    /* An event cannot cause itself... unless time is cyclic */
    for (int i = 0; i < event->effect_count; i++) {
        if (event->effects[i] == event->causal_id) {
            return 1;  /* Self-causation detected! */
        }
    }
    
    /* Check for causal loops */
    int visited[1000] = {0};
    int queue[1000];
    int front = 0, rear = 0;
    
    queue[rear++] = event->causal_id;
    visited[event->causal_id % 1000] = 1;
    
    while (front < rear) {
        int current = queue[front++];
        /* Would need event lookup table - simplified */
        if (current == event->causal_id && front > 1) {
            return 2;  /* Causal loop! */
        }
    }
    
    return 0;
}

/* 
   PART 5: QUANTUM GRAVITY MONTE CARLO SIMULATOR
    */

/* Path integral over metrics (simplified) */
double quantum_gravity_path_integral(Hyperpoint start, Hyperpoint end, 
                                     int samples, double temperature) {
    double sum_amplitude = 0.0;
    double sum_probability = 0.0;
    
    #pragma omp parallel for reduction(+:sum_amplitude,sum_probability)
    for (int s = 0; s < samples; s++) {
        /* Random walk through metric space */
        Hyperpoint current = start;
        double action = 0.0;
        
        for (int step = 0; step < 100; step++) {
            /* Random metric fluctuation */
            double metric_fluct[11][11];
            for (int i = 0; i < 11; i++) {
                for (int j = 0; j < 11; j++) {
                    metric_fluct[i][j] = 1.0 + 0.1 * ((double)rand()/RAND_MAX - 0.5);
                }
            }
            
            /* Calculate Einstein-Hilbert action with quantum corrections */
            double curvature_scalar = 0.0;
            for (int i = 0; i < 11; i++) {
                for (int j = 0; j < 11; j++) {
                    curvature_scalar += metric_fluct[i][j] * 
                                       current.curvature[i][j][j];
                }
            }
            
            action += curvature_scalar;
            
            /* Move toward endpoint with quantum fluctuations */
            for (int d = 0; d < 11; d++) {
                double dir = end.coord[d] - current.coord[d];
                double step_size = 0.01 * dir + 0.001 * ((double)rand()/RAND_MAX - 0.5);
                current.coord[d] += step_size;
            }
        }
        
        /* Boltzmann weight with NEGATIVE probabilities allowed */
        double amplitude = cexp(I * action);
        double probability = creal(amplitude * conj(amplitude));
        
        /* In parallel world, probabilities can exceed 1 or be negative */
        if (sin(start.coord[3] * end.coord[3]) > 0) {
            probability *= -1.0;
        }
        
        sum_amplitude += amplitude;
        sum_probability += probability;
    }
    
    return sum_probability / samples;
}

/* 
   PART 6: TOPOLOGICAL DEFECT GENERATOR
    */

/* Cosmic strings, monopoles, domain walls, textures */
typedef struct {
    int type;  /* 0=string, 1=monopole, 2=wall, 3=texture */
    Hyperpoint center;
    double tension;
    double length;          /* For strings */
    double magnetic_charge; /* For monopoles */
    double domain_orientation; /* For walls */
    int dimensionality;     /* 1D, 2D, 3D defects */
    double stability;
    atomic_int decay_count;
} TopologicalDefect;

/* Generate defect from vacuum manifold */
TopologicalDefect generate_defect(int homotopy_group) {
    TopologicalDefect defect;
    defect.type = rand() % 4;
    
    /* Random position in 11D */
    for (int i = 0; i < 11; i++) {
        defect.center.coord[i] = 10.0 * ((double)rand() / RAND_MAX - 0.5);
    }
    
    /* Defect properties based on homotopy group */
    switch (homotopy_group) {
        case 0:  /* π₀ - domain walls */
            defect.tension = 1e16;  /* GeV³ */
            defect.dimensionality = 3;
            break;
        case 1:  /* π₁ - cosmic strings */
            defect.tension = 1e15;
            defect.length = 1000.0 + 100.0 * ((double)rand() / RAND_MAX);
            defect.dimensionality = 1;
            break;
        case 2:  /* π₂ - monopoles */
            defect.magnetic_charge = 1.0 / 137.036;  /* 1/α */
            defect.dimensionality = 0;
            break;
        case 3:  /* π₃ - textures */
            defect.stability = 0.01 + 0.1 * ((double)rand() / RAND_MAX);
            defect.dimensionality = 3;
            break;
    }
    
    defect.decay_count = 0;
    return defect;
}

/* 
   PART 7: NON-LOCAL ENTANGLEMENT NETWORK
    */

/* Quantum entanglement that transcends spacetime */
typedef struct {
    int64_t particle_a;
    int64_t particle_b;
    double entanglement_strength;
    double* bell_state;  /* α|00⟩ + β|01⟩ + γ|10⟩ + δ|11⟩ */
    double last_measurement_time;
    int measurement_basis[3];
    _Atomic int collapsed;  /* Has wavefunction collapsed? */
    double hidden_variables[256];  /* Non-local hidden variables (Bohmian) */
} EntanglementLink;

/* Bell test that VIOLATES inequalities */
double perform_bell_test(EntanglementLink* link, int basis_a, int basis_b) {
    /* In our world, S ≤ 2. In parallel world, S can be > 2 or < -2 */
    double correlation = 0.0;
    
    /* Simulate measurement with SUPER-quantum correlations */
    for (int i = 0; i < 1000; i++) {
        double angle_a = basis_a * M_PI / 4.0;
        double angle_b = basis_b * M_PI / 4.0;
        
        /* Hidden variables that are NON-local and ACausal */
        double hidden = link->hidden_variables[i % 256];
        
        /* Measurement outcomes that violate Bell inequality maximally */
        double outcome_a = cos(2.0 * (angle_a + hidden)) > 0 ? 1 : -1;
        double outcome_b = cos(2.0 * (angle_b + hidden + 0.25)) > 0 ? 1 : -1;
        
        correlation += outcome_a * outcome_b;
    }
    
    correlation /= 1000.0;
    
    /* Sometimes give PRL (Popescu-Rohrlich) box correlations */
    if (sin(link->last_measurement_time) > 0.8) {
        correlation = (basis_a == basis_b) ? 1.0 : -1.0;  /* PRL box */
    }
    
    return correlation;
}

/* 
   PART 8: CHAOTIC INFLATION FIELD SIMULATOR
    */

/* Inflation with multiple fields and chaotic potentials */
typedef struct {
    double phi[7];           /* Inflation fields */
    double pi[7];            /* Conjugate momenta */
    double potential[7];     /* V(φ) */
    double epsilon;          /* Slow-roll parameter */
    double hubble;           /* Hubble parameter */
    double quantum_flucts[7]; /* Quantum fluctuations */
    int e_folds;             /* Number of e-folds */
    double reheating_temp;   /* Temperature after inflation */
} InflationField;

/* Chaotic inflation potential with multiple minima */
double inflation_potential(InflationField* inf, int field_idx) {
    double phi = inf->phi[field_idx];
    
    /* Multi-minima potential */
    switch (field_idx) {
        case 0:  /* Higgs-like */
            return 0.25 * pow(phi*phi - 246*246, 2);
        case 1:  /* Axion-like */
            return 1.0 - cos(phi / 100.0);
        case 2:  /* Monodromy */
            return pow(phi, 0.5) + 0.1 * sin(phi / 10.0);
        case 3:  /* Coleman-Weinberg */
            return phi*phi*phi*phi * (log(phi*phi / 10000.0) - 0.5);
        default: /* Random chaotic */
            return 0.1 * phi*phi + 0.01 * sin(phi) + 0.001 * sin(10*phi);
    }
}

/* Simulate inflation with quantum fluctuations */
void simulate_inflation(InflationField* inf, double dt) {
    #pragma omp parallel for
    for (int i = 0; i < 7; i++) {
        /* Equation of motion: φ̈ + 3Hφ̇ + dV/dφ = 0 */
        double dV_dphi;
        double phi = inf->phi[i];
        
        /* Numerical derivative of potential */
        double eps = 1e-6;
        inf->phi[i] += eps;
        double V_plus = inflation_potential(inf, i);
        inf->phi[i] -= 2*eps;
        double V_minus = inflation_potential(inf, i);
        inf->phi[i] = phi;  /* Restore */
        dV_dphi = (V_plus - V_minus) / (2*eps);
        
        /* Update equations */
        double acceleration = -3.0 * inf->hubble * inf->pi[i] - dV_dphi;
        inf->pi[i] += acceleration * dt;
        inf->phi[i] += inf->pi[i] * dt;
        
        /* Add quantum fluctuations (stochastic inflation) */
        double quantum_kick = 0.01 * ((double)rand()/RAND_MAX - 0.5) * 
                              sqrt(dt) * inf->hubble / (2*M_PI);
        inf->phi[i] += quantum_kick;
        inf->quantum_flucts[i] += quantum_kick * quantum_kick;
    }
    
    /* Update Hubble parameter */
    double energy_density = 0.0;
    for (int i = 0; i < 7; i++) {
        energy_density += 0.5 * inf->pi[i]*inf->pi[i] + 
                         inflation_potential(inf, i);
    }
    inf->hubble = sqrt(energy_density / 3.0);
    
    inf->e_folds += inf->hubble * dt;
}

/* 
   PART 9: PARALLEL WORLD RENDERER (ASCII ART)
    */

/* Render 11D spacetime projection to 2D terminal */
void render_parallel_world(Hyperpoint* viewer, Hyperpoint* objects, 
                          int num_objects, int width, int height) {
    /* We're projecting 11D -> 2D with dimensional compression */
    char screen[height][width];
    double depth_buffer[height][width];
    
    /* Initialize */
    #pragma omp parallel for collapse(2)
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            screen[y][x] = ' ';
            depth_buffer[y][x] = -1e100;
        }
    }
    
    /* Project each object */
    #pragma omp parallel for
    for (int o = 0; o < num_objects; o++) {
        Hyperpoint obj = objects[o];
        
        /* Calculate "distance" in warped space */
        double dist = parallel_distance(viewer, &obj);
        
        /* Project 11D to 2D with dimensional mixing */
        double proj_x = 0.0, proj_y = 0.0;
        for (int d = 0; d < 11; d++) {
            double mixing = sin(viewer->coord[3] * 0.1 * d);  /* Time-dependent */
            proj_x += mixing * obj.coord[d];
            proj_y += cos(viewer->coord[3] * 0.1 * d) * obj.coord[d];
        }
        
        /* Apply perspective (when it exists) */
        int screen_x = width/2 + (int)(proj_x * 100.0 / (dist + 1.0));
        int screen_y = height/2 + (int)(proj_y * 50.0 / (dist + 1.0));
        
        /* Clamp to screen */
        if (screen_x >= 0 && screen_x < width && 
            screen_y >= 0 && screen_y < height) {
            
            #pragma omp critical
            {
                if (dist > depth_buffer[screen_y][screen_x]) {
                    depth_buffer[screen_y][screen_x] = dist;
                    
                    /* Choose character based on dimensional properties */
                    int unfolded_dims = 0;
                    for (int d = 0; d < 11; d++) {
                        unfolded_dims += obj.dimension_active[d];
                    }
                    
                    char symbols[] = ".oO0@#*+~-";
                    screen[screen_y][screen_x] = symbols[unfolded_dims % 10];
                }
            }
        }
    }
    
    /* Render to terminal */
    printf("\033[2J\033[H");  /* Clear screen */
    printf("PARALLEL WORLD VIEWER | DIMENSIONS: 11 | TIME: %.3f\n\n", 
           viewer->coord[3]);
    
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            putchar(screen[y][x]);
        }
        putchar('\n');
    }
    
    printf("\nLEGEND: .(0D) o(1-2D) O(3-4D) 0(5-6D) @(7-8D) #(9-10D) *(11D)\n");
}

/* 
   PART 10: MAIN SIMULATION ENGINE
    */

/* Parallel world simulation parameters */
typedef struct {
    int total_dimensions;
    double time_step;
    double chaos_factor;  /* 0.0 = normal physics, 1.0 = maximum madness */
    int enable_ctc;       /* Closed timelike curves */
    int enable_ftl;       /* Faster-than-light travel */
    int enable_neg_mass;  /* Negative mass/energy */
    int enable_acausal;   /* Acausal effects */
    int quantum_gravity_level;
    int render_enabled;
    int width, height;
    atomic_long simulation_step;
    atomic_double world_entropy;
} SimulationParams;

/* The main simulation loop */
void* parallel_world_simulator(void* params_ptr) {
    SimulationParams* params = (SimulationParams*)params_ptr;
    
    /* Initialize parallel constants */
    ParallelConstants constants;
    atomic_store(&constants.pi, OUR_PI);
    atomic_store(&constants.c, OUR_C);
    atomic_store(&constants.G, OUR_G);
    atomic_store(&constants.hbar, OUR_HBAR);
    atomic_store(&constants.alpha, 1.0/137.036);
    
    /* Create observer */
    Hyperpoint observer;
    memset(&observer, 0, sizeof(Hyperpoint));
    observer.coord[3] = 0.0;  /* Initial time */
    for (int d = 0; d < 11; d++) {
        observer.dimension_active[d] = (d < 4) ? 1 : 0;  /* First 4 dimensions active */
    }
    
    /* Generate random objects in 11D space */
    int num_objects = 1000;
    Hyperpoint* objects = aligned_alloc(64, num_objects * sizeof(Hyperpoint));
    
    #pragma omp parallel for
    for (int i = 0; i < num_objects; i++) {
        for (int d = 0; d < 11; d++) {
            objects[i].coord[d] = 20.0 * ((double)rand() / RAND_MAX - 0.5);
        }
        /* Randomly unfold extra dimensions */
        for (int d = 4; d < 11; d++) {
            objects[i].dimension_active[d] = (rand() % 100 < 30) ? 1 : 0;
        }
    }
    
    /* Main simulation loop */
    while (atomic_load(&params->simulation_step) < 1000000) {
        long step = atomic_fetch_add(&params->simulation_step, 1);
        
        /* Update physics constants (they FLUCTUATE) */
        double time = step * params->time_step;
        atomic_store(&constants.pi, OUR_PI * (1.0 + 0.1 * sin(time) * params->chaos_factor));
        atomic_store(&constants.c, OUR_C * (1.0 + 0.5 * sin(time * 0.7) * params->chaos_factor));
        atomic_store(&constants.G, OUR_G * (1.0 + 2.0 * sin(time * 0.3) * params->chaos_factor));
        
        /* Move observer through non-Euclidean space */
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                /* Chaotic motion in regular dimensions */
                observer.coord[0] += 0.01 * sin(time * 1.1);
                observer.coord[1] += 0.01 * sin(time * 1.3);
                observer.coord[2] += 0.01 * sin(time * 1.7);
                observer.coord[3] += params->time_step;
            }
            
            #pragma omp section
            {
                /* Random unfolding/folding of extra dimensions */
                if (rand() % 1000 == 0) {
                    int dim = rand() % 7 + 4;
                    observer.dimension_active[dim] = !observer.dimension_active[dim];
                }
            }
        }
        
        /* Update object positions with non-local effects */
        #pragma omp parallel for
        for (int i = 0; i < num_objects; i++) {
            /* Standard motion */
            for (int d = 0; d < 11; d++) {
                if (observer.dimension_active[d]) {
                    objects[i].coord[d] += 0.001 * ((double)rand() / RAND_MAX - 0.5);
                }
            }
            
            /* Non-local quantum jumps (teleportation) */
            if (params->chaos_factor > 0.5 && (rand() % 10000) < 10) {
                int jump_dim = rand() % 11;
                objects[i].coord[jump_dim] = 20.0 * ((double)rand() / RAND_MAX - 0.5);
            }
            
            /* Dimensional phase transitions */
            if ((rand() % 1000) < 5) {
                for (int d = 4; d < 11; d++) {
                    objects[i].dimension_active[d] = rand() % 2;
                }
            }
        }
        
        /* Calculate world entropy (can decrease due to time reversal) */
        double entropy = 0.0;
        #pragma omp parallel for reduction(+:entropy)
        for (int i = 0; i < num_objects; i++) {
            for (int j = i+1; j < num_objects; j++) {
                double dist = parallel_distance(&objects[i], &objects[j]);
                if (dist > 0) {
                    entropy += log(dist);
                } else {
                    entropy += log(-dist);  /* Negative distances contribute differently */
                }
            }
        }
        
        atomic_store(&params->world_entropy, entropy);
        
        /* Render if enabled */
        if (params->render_enabled && (step % 100 == 0)) {
            render_parallel_world(&observer, objects, num_objects, 
                                 params->width, params->height);
            usleep(50000);  /* 50ms delay for visualization */
        }
        
        /* Print statistics occasionally */
        if (step % 1000 == 0) {
            printf("\n=== PARALLEL WORLD STATUS ===\n");
            printf("Step: %ld | Time: %.3f | Entropy: %.3f\n", 
                   step, time, entropy);
            printf("Constants: π=%.6f c=%.3e G=%.3e\n",
                   atomic_load(&constants.pi),
                   atomic_load(&constants.c),
                   atomic_load(&constants.G));
            printf("Active dimensions: ");
            for (int d = 0; d < 11; d++) {
                if (observer.dimension_active[d]) printf("%d ", d);
            }
            printf("\n");
        }
    }
    
    free(objects);
    return NULL;
}

/* 
   PART 11: BIZARRE PHYSICS DEMONSTRATIONS
    */

/* Demonstration 1: Time-reversed thermodynamics */
void time_reversed_heat_flow() {
    printf("\n=== TIME-REVERSED THERMODYNAMICS ===\n");
    
    double temp_hot = 100.0, temp_cold = 0.0;
    
    for (int i = 0; i < 20; i++) {
        /* In parallel world, heat can flow from cold to hot */
        double flow_direction = (sin(i * 0.5) > 0) ? 1.0 : -1.0;
        double heat_flow = 0.1 * (temp_cold - temp_hot) * flow_direction;
        
        temp_hot -= heat_flow;
        temp_cold += heat_flow;
        
        printf("Step %2d: Hot=%.1f Cold=%.1f Flow=%+.3f\n",
               i, temp_hot, temp_cold, heat_flow);
    }
}

/* Demonstration 2: Negative refractive index materials */
void metamaterial_light_bending() {
    printf("\n=== NEGATIVE REFRACTIVE INDEX ===\n");
    
    double n1 = 1.0;  /* Air */
    double n2 = -1.5; /* Parallel world material (NEGATIVE index) */
    
    for (double angle_deg = 0; angle_deg < 90; angle_deg += 10) {
        double angle_rad = angle_deg * M_PI / 180.0;
        
        /* Snell's law: n1 sinθ1 = n2 sinθ2 */
        double sin_theta2 = (n1 * sin(angle_rad)) / n2;
        
        if (fabs(sin_theta2) <= 1.0) {
            double theta2 = asin(sin_theta2) * 180.0 / M_PI;
            printf("Incident: %5.1f° -> Refracted: %7.1f°\n", 
                   angle_deg, theta2);
        } else {
            printf("Incident: %5.1f° -> TOTAL INTERNAL REFLECTION (backward)\n",
                   angle_deg);
        }
    }
}

/* Demonstration 3: Superluminal communication */
void ftl_communication_test() {
    printf("\n=== FTL COMMUNICATION TEST ===\n");
    
    double distance_ly = 10.0;  /* 10 light years */
    double message_sent_time = 0.0;
    
    for (double v_over_c = 1.0; v_over_c <= 10.0; v_over_c += 1.0) {
        /* Lorentz transformation breaks down for v > c */
        double gamma = 1.0 / sqrt(1.0 - 1.0/(v_over_c*v_over_c));  /* IMAGINARY! */
        
        /* Message arrival BEFORE it was sent? */
        double arrival_time = message_sent_time + distance_ly / v_over_c;
        
        printf("v/c = %4.1f | Gamma = %7.3f | Arrival: t = %.3f yr\n",
               v_over_c, gamma, arrival_time);
        
        if (arrival_time < message_sent_time) {
            printf("  *** PARADOX: Message arrived before sending! ***\n");
        }
    }
}

/* 
   PART 12: MAIN FUNCTION WITH ARGUMENT PARSING
    */

int main(int argc, char** argv) {
    printf("====\n");
    printf("PARALLEL WORLD SIMULATOR v2.0\n");
    printf("Physics-Deviant Reality Engine\n");
    printf("====\n");
    
    /* Parse command line arguments */
    SimulationParams params;
    params.total_dimensions = 11;
    params.time_step = 0.01;
    params.chaos_factor = 0.7;
    params.enable_ctc = 1;
    params.enable_ftl = 1;
    params.enable_neg_mass = 1;
    params.enable_acausal = 1;
    params.quantum_gravity_level = 3;
    params.render_enabled = 0;
    params.width = 80;
    params.height = 40;
    atomic_init(&params.simulation_step, 0);
    atomic_init(&params.world_entropy, 0.0);
    
    int seed = time(NULL);
    
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--seed") == 0 && i+1 < argc) {
            seed = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--dimensions") == 0 && i+1 < argc) {
            params.total_dimensions = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--chaos") == 0 && i+1 < argc) {
            params.chaos_factor = atof(argv[++i]);
        } else if (strcmp(argv[i], "--render") == 0) {
            params.render_enabled = 1;
        } else if (strcmp(argv[i], "--demo") == 0) {
            /* Run demonstrations */
            time_reversed_heat_flow();
            metamaterial_light_bending();
            ftl_communication_test();
            return 0;
        }
    }
    
    srand(seed);
    printf("Initializing with seed: %d\n", seed);
    printf("Chaos factor: %.2f\n", params.chaos_factor);
    printf("Dimensions: %d\n", params.total_dimensions);
    printf("CTC enabled: %s\n", params.enable_ctc ? "YES" : "NO");
    printf("FTL enabled: %s\n", params.enable_ftl ? "YES" : "NO");
    printf("Negative mass: %s\n", params.enable_neg_mass ? "YES" : "NO");
    printf("Acausal effects: %s\n", params.enable_acausal ? "YES" : "NO");
    
    /* Create simulation thread */
    pthread_t sim_thread;
    pthread_create(&sim_thread, NULL, parallel_world_simulator, &params);
    
    /* Wait for simulation to complete or user interrupt */
    printf("\nSimulation running... Press Ctrl+C to stop\n");
    
    pthread_join(sim_thread, NULL);
    
    printf("\n=== SIMULATION COMPLETE ===\n");
    printf("Final entropy: %.3f\n", atomic_load(&params.world_entropy));
    printf("Total steps: %ld\n", atomic_load(&params.simulation_step));
    
    return 0;
}
