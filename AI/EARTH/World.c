/**
 * earth_physics.c
 * A complete Earth physics simulation system covering multiple physical domains
 * Compile: gcc -o earth_physics earth_physics.c -lm -lGL -lGLU -lglut -O3
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>

/* ==============================================
   PART 1: FOUNDATIONAL PHYSICS ENGINE
   ============================================== */

typedef struct { double x, y, z; } Vector3;
typedef struct { double lat, lon, alt; } GeoCoord;  // Latitude, Longitude, Altitude

/* Core physics constants (SI units) */
const double G = 6.67430e-11;      // Gravitational constant
const double R_EARTH = 6371000.0;  // Earth radius (m)
const double M_EARTH = 5.9722e24;  // Earth mass (kg)
const double OMEGA_EARTH = 7.2921159e-5;  // Earth angular velocity (rad/s)

/* Vector operations */
Vector3 vec_add(Vector3 a, Vector3 b) { return (Vector3){a.x+b.x, a.y+b.y, a.z+b.z}; }
Vector3 vec_sub(Vector3 a, Vector3 b) { return (Vector3){a.x-b.x, a.y-b.y, a.z-b.z}; }
Vector3 vec_mul(Vector3 v, double s) { return (Vector3){v.x*s, v.y*s, v.z*s}; }
double vec_dot(Vector3 a, Vector3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
double vec_length(Vector3 v) { return sqrt(v.x*v.x + v.y*v.y + v.z*v.z); }
Vector3 vec_normalize(Vector3 v) { double l = vec_length(v); return (Vector3){v.x/l, v.y/l, v.z/l}; }

/* Coordinate transformations */
Vector3 geo_to_cartesian(GeoCoord geo) {
    double lat_rad = geo.lat * M_PI/180.0;
    double lon_rad = geo.lon * M_PI/180.0;
    double r = R_EARTH + geo.alt;
    return (Vector3){
        r * cos(lat_rad) * cos(lon_rad),
        r * cos(lat_rad) * sin(lon_rad),
        r * sin(lat_rad)
    };
}

GeoCoord cartesian_to_geo(Vector3 pos) {
    double r = vec_length(pos);
    double lat = asin(pos.z / r) * 180.0/M_PI;
    double lon = atan2(pos.y, pos.x) * 180.0/M_PI;
    return (GeoCoord){lat, lon, r - R_EARTH};
}

/* ==============================================
   PART 2: GEOPHYSICS & EARTH STRUCTURE
   ============================================== */

typedef struct {
    char name[32];
    double depth_start;    // km from surface
    double depth_end;      // km from surface
    double density;        // kg/m³
    double temperature;    // K
    double pressure;       // Pa
    char composition[64];
} EarthLayer;

EarthLayer earth_layers[] = {
    {"Crust", 0, 35, 2700, 293, 101325, "Granite, Basalt"},
    {"Upper Mantle", 35, 410, 3400, 1673, 1.4e9, "Peridotite"},
    {"Transition Zone", 410, 660, 3700, 2073, 2.4e10, "Ringwoodite"},
    {"Lower Mantle", 660, 2891, 4500, 3073, 1.4e11, "Bridgmanite"},
    {"Outer Core", 2891, 5150, 9900, 4073, 3.3e11, "Liquid Iron-Nickel"},
    {"Inner Core", 5150, 6371, 13000, 5700, 3.6e12, "Solid Iron"}
};

/* Calculate pressure at given depth using PREM model */
double pressure_at_depth(double depth_km) {
    double depth_m = depth_km * 1000;
    // Simplified polynomial approximation of PREM model
    double p = 101325;  // Surface pressure
    if (depth_km > 0) {
        p += 3.3e7 * depth_km;  // Crustal gradient
    }
    if (depth_km > 35) {
        p += 1.1e8 * (depth_km - 35);  // Mantle gradient
    }
    if (depth_km > 2891) {
        p += 5.5e7 * (depth_km - 2891);  // Core gradient
    }
    return p;
}

/* Calculate temperature at depth using geotherm */
double temperature_at_depth(double depth_km) {
    double surface_temp = 293;  // 20°C
    double mantle_temp = 1673;  // 1400°C at 35km
    double core_temp = 5700;    // Inner core boundary
    
    if (depth_km <= 35) {
        return surface_temp + (mantle_temp - surface_temp) * (depth_km / 35);
    } else if (depth_km <= 2891) {
        return mantle_temp + (core_temp - mantle_temp) * ((depth_km - 35) / (2891 - 35));
    } else {
        return core_temp;
    }
}

/* ==============================================
   PART 3: ORBITAL MECHANICS & GRAVITY
   ============================================== */

/* Two-body problem: Calculate gravitational force */
Vector3 gravitational_force(Vector3 pos1, double mass1, Vector3 pos2, double mass2) {
    Vector3 r_vec = vec_sub(pos2, pos1);
    double r = vec_length(r_vec);
    if (r < 1e-6) return (Vector3){0,0,0};
    
    double force_mag = G * mass1 * mass2 / (r * r);
    Vector3 force_dir = vec_normalize(r_vec);
    return vec_mul(force_dir, force_mag);
}

/* Calculate orbital velocity for circular orbit at altitude */
double orbital_velocity(double altitude) {
    double r = R_EARTH + altitude;
    return sqrt(G * M_EARTH / r);
}

/* Calculate orbital period (Kepler's third law) */
double orbital_period(double altitude) {
    double r = R_EARTH + altitude;
    return 2 * M_PI * sqrt(r * r * r / (G * M_EARTH));
}

/* 4th-order Runge-Kutta integration for orbital motion */
typedef struct { Vector3 pos, vel; double mass; } Satellite;

Satellite rk4_orbit_step(Satellite sat, double dt, Vector3 (*force_func)(Vector3, Vector3, double)) {
    Vector3 k1_v = force_func(sat.pos, sat.vel, sat.mass);
    Vector3 k1_r = sat.vel;
    
    Vector3 k2_v = force_func(vec_add(sat.pos, vec_mul(k1_r, dt/2)), 
                              vec_add(sat.vel, vec_mul(k1_v, dt/2)), sat.mass);
    Vector3 k2_r = vec_add(sat.vel, vec_mul(k1_v, dt/2));
    
    Vector3 k3_v = force_func(vec_add(sat.pos, vec_mul(k2_r, dt/2)), 
                              vec_add(sat.vel, vec_mul(k2_v, dt/2)), sat.mass);
    Vector3 k3_r = vec_add(sat.vel, vec_mul(k2_v, dt/2));
    
    Vector3 k4_v = force_func(vec_add(sat.pos, vec_mul(k3_r, dt)), 
                              vec_add(sat.vel, vec_mul(k3_v, dt)), sat.mass);
    Vector3 k4_r = vec_add(sat.vel, vec_mul(k3_v, dt));
    
    Satellite result;
    result.pos = vec_add(sat.pos, vec_mul(vec_add(vec_add(k1_r, vec_mul(k2_r, 2)), 
                        vec_add(vec_mul(k3_r, 2), k4_r)), dt/6));
    result.vel = vec_add(sat.vel, vec_mul(vec_add(vec_add(k1_v, vec_mul(k2_v, 2)), 
                        vec_add(vec_mul(k3_v, 2), k4_v)), dt/6));
    result.mass = sat.mass;
    return result;
}

/* Earth's gravitational force with J2 perturbation (oblateness) */
Vector3 earth_gravity_with_j2(Vector3 pos, Vector3 vel, double mass) {
    double r = vec_length(pos);
    Vector3 r_hat = vec_normalize(pos);
    
    // Central gravity
    Vector3 g_central = vec_mul(r_hat, -G * M_EARTH / (r * r));
    
    // J2 perturbation (Earth's oblateness)
    double J2 = 1.08263e-3;
    double z = pos.z;
    double r2 = r * r;
    double r5 = r2 * r2 * r;
    
    Vector3 g_j2 = (Vector3){
        5 * z * z / r2 - 1,
        5 * z * z / r2 - 1,
        5 * z * z / r2 - 3
    };
    g_j2 = vec_mul(g_j2, 1.5 * J2 * G * M_EARTH * R_EARTH * R_EARTH / r5);
    
    // Coriolis and centrifugal forces (rotating frame)
    Vector3 omega = {0, 0, OMEGA_EARTH};
    Vector3 coriolis = vec_mul(vec_cross(vel, omega), 2);
    Vector3 centrifugal = vec_mul(vec_cross(omega, vec_cross(omega, pos)), 1);
    
    return vec_add(vec_add(g_central, g_j2), vec_add(coriolis, centrifugal));
}

Vector3 vec_cross(Vector3 a, Vector3 b) {
    return (Vector3){
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    };
}

/* ==============================================
   PART 4: ATMOSPHERIC PHYSICS & FLUID DYNAMICS
   ============================================== */

/* International Standard Atmosphere model */
typedef struct {
    double altitude;      // m
    double temperature;   // K
    double pressure;      // Pa
    double density;       // kg/m³
    double speed_of_sound; // m/s
} AtmosphereLayer;

/* Calculate atmospheric properties using ISA model */
AtmosphereLayer atmosphere_at_altitude(double altitude_m) {
    AtmosphereLayer layer;
    layer.altitude = altitude_m;
    
    if (altitude_m <= 11000) {  // Troposphere
        layer.temperature = 288.15 - 0.0065 * altitude_m;
        layer.pressure = 101325 * pow(1 - altitude_m/44330, 5.256);
    } else if (altitude_m <= 20000) {  // Lower Stratosphere
        layer.temperature = 216.65;
        layer.pressure = 22632 * exp(-(altitude_m - 11000)/6340);
    } else if (altitude_m <= 32000) {  // Upper Stratosphere
        layer.temperature = 216.65 + 0.001 * (altitude_m - 20000);
        layer.pressure = 5474.9 * pow(216.65/layer.temperature, 34.163);
    } else {
        layer.temperature = 228.65 + 0.0028 * (altitude_m - 32000);
        layer.pressure = 868.02 * pow(228.65/layer.temperature, 12.0);
    }
    
    layer.density = layer.pressure / (287.05 * layer.temperature);
    layer.speed_of_sound = sqrt(1.4 * 287.05 * layer.temperature);
    return layer;
}

/* Calculate drag force on an object */
double drag_force(double density, double velocity, double drag_coeff, double area) {
    return 0.5 * density * velocity * velocity * drag_coeff * area;
}

/* Bernoulli's principle for fluid flow */
double bernoulli_equation(double pressure, double density, double velocity, double height) {
    return pressure + 0.5 * density * velocity * velocity + density * 9.81 * height;
}

/* ==============================================
   PART 5: THERMODYNAMICS & HEAT TRANSFER
   ============================================== */

/* Heat transfer models */
typedef struct {
    double thermal_conductivity;  // W/m·K
    double specific_heat;         // J/kg·K
    double density;               // kg/m³
} Material;

Material materials[] = {
    {"Granite", 2.9, 790, 2700},
    {"Water", 0.6, 4186, 1000},
    {"Air", 0.026, 1005, 1.2},
    {"Iron", 80.2, 449, 7874}
};

/* Fourier's Law of Heat Conduction */
double heat_conduction(double k, double area, double dT, double thickness) {
    return k * area * dT / thickness;
}

/* Calculate temperature distribution over time (1D heat equation) */
void simulate_heat_transfer(double* temp, int n, double dx, double dt, double alpha, int steps) {
    double* new_temp = malloc(n * sizeof(double));
    
    for (int s = 0; s < steps; s++) {
        for (int i = 1; i < n-1; i++) {
            new_temp[i] = temp[i] + alpha * dt/(dx*dx) * (temp[i+1] - 2*temp[i] + temp[i-1]);
        }
        // Boundary conditions (insulated)
        new_temp[0] = new_temp[1];
        new_temp[n-1] = new_temp[n-2];
        
        memcpy(temp, new_temp, n * sizeof(double));
    }
    free(new_temp);
}

/* ==============================================
   PART 6: ELECTROMAGNETISM & EARTH'S MAGNETIC FIELD
   ============================================== */

/* International Geomagnetic Reference Field (simplified) */
Vector3 earth_magnetic_field(GeoCoord geo, double year) {
    double lat_rad = geo.lat * M_PI/180.0;
    double lon_rad = geo.lon * M_PI/180.0;
    
    // Dipole approximation
    double B0 = 3.12e-5;  // Tesla at equator
    double r = (R_EARTH + geo.alt) / R_EARTH;
    
    double Br = -2 * B0 * pow(1/r, 3) * sin(lat_rad);
    double Btheta = -B0 * pow(1/r, 3) * cos(lat_rad);
    
    // Convert to Cartesian (ENU coordinates)
    Vector3 B_enu = {
        -Btheta * cos(lon_rad) - Br * sin(lat_rad) * cos(lon_rad),
        -Btheta * sin(lon_rad) - Br * sin(lat_rad) * sin(lon_rad),
        Br * cos(lat_rad)
    };
    
    return B_enu;
}

/* Lorentz force on charged particle in magnetic field */
Vector3 lorentz_force(double charge, Vector3 velocity, Vector3 B_field) {
    return vec_mul(vec_cross(velocity, B_field), charge);
}

/* ==============================================
   PART 7: SEISMIC WAVE PROPAGATION
   ============================================== */

typedef struct {
    double P_wave_speed;   // Primary wave speed (m/s)
    double S_wave_speed;   // Secondary wave speed (m/s)
    double density;        // Medium density (kg/m³)
} SeismicMedium;

SeismicMedium crust_medium = {5800, 3400, 2700};
SeismicMedium mantle_medium = {8000, 4500, 3400};

/* Calculate seismic wave arrival times */
double seismic_arrival_time(double distance_km, double depth_km, bool is_p_wave) {
    double v = is_p_wave ? crust_medium.P_wave_speed : crust_medium.S_wave_speed;
    
    // Simple direct path calculation
    double distance_m = distance_km * 1000;
    double depth_m = depth_km * 1000;
    
    double path_length = sqrt(distance_m*distance_m + depth_m*depth_m);
    return path_length / v;
}

/* ==============================================
   PART 8: VISUALIZATION SYSTEM (ASCII/OpenGL)
   ============================================== */

#ifdef WITH_OPENGL
#include <GL/glut.h>

void render_earth() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // Draw Earth as sphere
    glColor3f(0.0, 0.0, 1.0);
    glutSolidSphere(1.0, 50, 50);
    
    // Draw atmosphere layers
    glColor4f(0.0, 1.0, 1.0, 0.3);
    for(int i=0; i<5; i++) {
        glPushMatrix();
        glScalef(1.0 + i*0.1, 1.0 + i*0.1, 1.0 + i*0.1);
        glutWireSphere(1.0, 20, 20);
        glPopMatrix();
    }
    
    glutSwapBuffers();
}
#endif

/* ASCII visualization of Earth layers */
void print_earth_cross_section() {
    printf("\n=== EARTH CROSS SECTION ===\n");
    printf("Depth (km)  Layer           Density Temp(K) Pressure\n");
    printf("----------------------------------------------------\n");
    
    for(int i=0; i<6; i++) {
        EarthLayer layer = earth_layers[i];
        printf("%4.0f-%-6.0f %-15s %6.0f   %5.0f   %.1e Pa\n",
               layer.depth_start, layer.depth_end,
               layer.name, layer.density,
               layer.temperature, layer.pressure);
    }
}

/* ==============================================
   PART 9: INTEGRATED SIMULATION EXAMPLES
   ============================================== */

/* Example 1: Satellite orbit simulation */
void simulate_satellite_orbit() {
    printf("\n=== SATELLITE ORBIT SIMULATION ===\n");
    
    Satellite sat;
    sat.pos = (Vector3){R_EARTH + 400000, 0, 0};  // 400km altitude
    sat.vel = (Vector3){0, orbital_velocity(400000), 0};
    sat.mass = 1000;  // 1 ton satellite
    
    double dt = 10;  // 10 second steps
    int steps = 1000;
    
    for(int i=0; i<steps; i++) {
        sat = rk4_orbit_step(sat, dt, earth_gravity_with_j2);
        
        if(i % 100 == 0) {
            GeoCoord geo = cartesian_to_geo(sat.pos);
            printf("Time: %4ds | Altitude: %6.0fkm | Lat: %6.1f° | Lon: %6.1f°\n",
                   i*10, geo.alt/1000, geo.lat, geo.lon);
        }
    }
}

/* Example 2: Atmospheric entry simulation */
void simulate_atmospheric_entry() {
    printf("\n=== ATMOSPHERIC ENTRY SIMULATION ===\n");
    
    double altitude = 120000;  // Start at 120km
    double velocity = 7800;    // Orbital velocity
    double mass = 5000;        // Spacecraft mass
    double drag_area = 15;     // Cross-sectional area
    double drag_coeff = 1.5;
    
    double dt = 0.1;
    double time = 0;
    
    while(altitude > 0 && time < 600) {
        AtmosphereLayer atm = atmosphere_at_altitude(altitude);
        double drag = drag_force(atm.density, velocity, drag_coeff, drag_area);
        double gravity = 9.81 * mass * pow(R_EARTH/(R_EARTH+altitude), 2);
        
        double acceleration = (drag - gravity) / mass;
        velocity += acceleration * dt;
        altitude -= velocity * dt;
        time += dt;
        
        if((int)time % 10 == 0) {
            printf("Time: %4.1fs | Alt: %6.0fm | Vel: %6.0f m/s | Temp: %5.0fK\n",
                   time, altitude, velocity, atm.temperature);
        }
    }
}

/* Example 3: Seismic wave propagation */
void simulate_seismic_event() {
    printf("\n=== SEISMIC EVENT SIMULATION ===\n");
    
    double epicenter_lat = 35.0;
    double epicenter_lon = 140.0;
    double depth = 10.0;  // km
    
    // Calculate arrival times at different distances
    printf("Distance from epicenter (km) | P-wave arrival | S-wave arrival\n");
    printf("--------------------------------------------------------------\n");
    
    for(double distance=10; distance<=1000; distance*=2) {
        double p_time = seismic_arrival_time(distance, depth, true);
        double s_time = seismic_arrival_time(distance, depth, false);
        
        printf("%8.0f km                  | %6.1f s       | %6.1f s\n",
               distance, p_time, s_time);
    }
}

/* Example 4: Geothermal gradient */
void simulate_geothermal() {
    printf("\n=== GEOTHERMAL GRADIENT SIMULATION ===\n");
    
    int n_points = 100;
    double depth_max = 50;  // km
    double dx = depth_max / n_points;
    
    // Initialize temperature array (surface to depth)
    double* temperature = malloc(n_points * sizeof(double));
    for(int i=0; i<n_points; i++) {
        double depth = i * dx;
        temperature[i] = temperature_at_depth(depth);
    }
    
    // Simulate heat conduction
    double thermal_diffusivity = 1e-6;  // m²/s for crust
    double dt = 3600 * 24 * 365;        // 1 year in seconds
    double dx_m = dx * 1000;           // Convert km to m
    
    simulate_heat_transfer(temperature, n_points, dx_m, dt, thermal_diffusivity, 10);
    
    printf("Depth(km) Temperature(K) Pressure(Pa)\n");
    printf("--------------------------------------\n");
    for(int i=0; i<n_points; i+=10) {
        double depth = i * dx;
        printf("%8.1f %13.0f %.2e\n", 
               depth, temperature[i], pressure_at_depth(depth));
    }
    
    free(temperature);
}

/* ==============================================
   PART 10: MAIN CONTROL & USER INTERFACE
   ============================================== */

void print_menu() {
    printf("\n===== EARTH PHYSICS SIMULATION SYSTEM =====\n");
    printf("1. View Earth Structure\n");
    printf("2. Simulate Satellite Orbit\n");
    printf("3. Simulate Atmospheric Entry\n");
    printf("4. Simulate Seismic Event\n");
    printf("5. Simulate Geothermal Gradient\n");
    printf("6. Calculate Orbital Parameters\n");
    printf("7. Atmosphere Profile\n");
    printf("8. Magnetic Field Calculator\n");
    printf("9. Run All Simulations\n");
    printf("0. Exit\n");
    printf("===========================================\n");
    printf("Choice: ");
}

int main() {
    srand(time(NULL));
    
    int choice;
    
    do {
        print_menu();
        scanf("%d", &choice);
        
        switch(choice) {
            case 1:
                print_earth_cross_section();
                break;
            case 2:
                simulate_satellite_orbit();
                break;
            case 3:
                simulate_atmospheric_entry();
                break;
            case 4:
                simulate_seismic_event();
                break;
            case 5:
                simulate_geothermal();
                break;
            case 6: {
                printf("\n=== ORBITAL PARAMETERS ===\n");
                for(double alt=200; alt<=36000; alt*=2) {
                    printf("Altitude: %6.0fkm | Velocity: %6.0f m/s | Period: %6.1f min\n",
                           alt, orbital_velocity(alt*1000), 
                           orbital_period(alt*1000)/60);
                }
                break;
            }
            case 7: {
                printf("\n=== ATMOSPHERE PROFILE ===\n");
                printf("Altitude(m) Temp(K) Pressure(Pa) Density(kg/m³) Speed of Sound(m/s)\n");
                printf("--------------------------------------------------------------------\n");
                for(double alt=0; alt<=100000; alt+=10000) {
                    AtmosphereLayer atm = atmosphere_at_altitude(alt);
                    printf("%9.0f %7.0f %.3e %12.3e %17.0f\n",
                           alt, atm.temperature, atm.pressure, 
                           atm.density, atm.speed_of_sound);
                }
                break;
            }
            case 8: {
                printf("\n=== MAGNETIC FIELD CALCULATOR ===\n");
                GeoCoord test_points[] = {
                    {0, 0, 0},      // Equator
                    {45, 0, 0},     // Mid-latitude
                    {90, 0, 0}      // North pole
                };
                for(int i=0; i<3; i++) {
                    Vector3 B = earth_magnetic_field(test_points[i], 2024);
                    printf("Lat: %3.0f° | Field: (%6.1f, %6.1f, %6.1f) nT\n",
                           test_points[i].lat, 
                           B.x*1e9, B.y*1e9, B.z*1e9);
                }
                break;
            }
            case 9:
                printf("\n=== RUNNING ALL SIMULATIONS ===\n");
                print_earth_cross_section();
                simulate_satellite_orbit();
                simulate_atmospheric_entry();
                simulate_seismic_event();
                simulate_geothermal();
                break;
            case 0:
                printf("Exiting Earth Physics Simulation System.\n");
                break;
            default:
                printf("Invalid choice. Please try again.\n");
        }
        
    } while(choice != 0);
    
    return 0;
}