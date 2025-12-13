/**
 * multiverse.c - A complete MultiVerse simulation with enhanced math,
 * persistent state, networking, and visualization in a single C file.
 * 
 * Architecture: Multiple parallel universes (Multiverse) containing
 * their own Metaverse instances with mathematical, persistent, networked,
 * and visual dimensions.
 * 
 * Compile: gcc -o multiverse multiverse.c -lm -lSDL2 -pthread -O3
 * Run Server: ./multiverse --server --port 8080 --universe-count 3
 * Run Client: ./multiverse --client --host 127.0.0.1 --port 8080
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <stdatomic.h>
#include <signal.h>
#include <stdbool.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <fcntl.h>

/* === 1. ENHANCED MATHEMATICAL LIBRARY (Fractals, Quaternions, Calculus, Chaos) === */

/* ---- 1.1 Vector/Matrix Library (Extended) ---- */
typedef struct { double x, y, z, w; } Vec4;
typedef struct { double x, y, z; } Vec3;
typedef struct { double x, y; } Vec2;
typedef struct { double m[4][4]; } Mat4;
typedef struct { double w, x, y, z; } Quat;

/* Vector operations */
Vec3 vec3_add(Vec3 a, Vec3 b) { return (Vec3){a.x+b.x, a.y+b.y, a.z+b.z}; }
Vec3 vec3_sub(Vec3 a, Vec3 b) { return (Vec3){a.x-b.x, a.y-b.y, a.z-b.z}; }
Vec3 vec3_mul(Vec3 a, double s) { return (Vec3){a.x*s, a.y*s, a.z*s}; }
double vec3_dot(Vec3 a, Vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
Vec3 vec3_cross(Vec3 a, Vec3 b) {
    return (Vec3){a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}
double vec3_length(Vec3 v) { return sqrt(v.x*v.x + v.y*v.y + v.z*v.z); }

/* Quaternion operations */
Quat quat_mul(Quat a, Quat b) {
    return (Quat){
        a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
        a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
        a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x,
        a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w
    };
}
Quat quat_from_axis_angle(Vec3 axis, double angle) {
    double s = sin(angle/2);
    return (Quat){cos(angle/2), axis.x*s, axis.y*s, axis.z*s};
}
Vec3 quat_rotate(Vec3 v, Quat q) {
    Quat p = {0, v.x, v.y, v.z};
    Quat q_conj = {q.w, -q.x, -q.y, -q.z};
    Quat result = quat_mul(quat_mul(q, p), q_conj);
    return (Vec3){result.x, result.y, result.z};
}

/* Matrix operations */
Mat4 mat4_identity() {
    Mat4 m = {{{0}}};
    for(int i=0; i<4; i++) m.m[i][i] = 1.0;
    return m;
}
Mat4 mat4_perspective(double fov, double aspect, double near, double far) {
    double f = 1.0 / tan(fov/2);
    Mat4 m = {{{0}}};
    m.m[0][0] = f / aspect;
    m.m[1][1] = f;
    m.m[2][2] = (far+near)/(near-far);
    m.m[2][3] = -1;
    m.m[3][2] = (2*far*near)/(near-far);
    return m;
}

/* ---- 1.2 Fractal Mathematics (Mandelbrot, Julia, 3D Fractals) ---- */
typedef struct {
    double real, imag;
} Complex;

Complex complex_mul(Complex a, Complex b) {
    return (Complex){a.real*b.real - a.imag*b.imag, a.real*b.imag + a.imag*b.real};
}
Complex complex_add(Complex a, Complex b) {
    return (Complex){a.real+b.real, a.imag+b.imag};
}

/* Mandelbrot set iteration (returns iteration count) */
int mandelbrot_iter(Complex c, int max_iter) {
    Complex z = {0, 0};
    for(int i=0; i<max_iter; i++) {
        if(z.real*z.real + z.imag*z.imag > 4.0) return i;
        z = complex_add(complex_mul(z, z), c);
    }
    return max_iter;
}

/* 3D Fractal noise (Perlin-inspired for terrain generation) */
double fractal_noise_3d(double x, double y, double z, int octaves, double persistence) {
    double total = 0;
    double frequency = 1;
    double amplitude = 1;
    double max_value = 0;
    
    for(int i=0; i<octaves; i++) {
        // Simple pseudo-random hash for each octave
        double nx = sin(x * frequency * 12.9898 + y * 78.233 + z * 37.719) * 43758.5453;
        double ny = sin(y * frequency * 78.233 + z * 12.9898 + x * 37.719) * 43758.5453;
        double nz = sin(z * frequency * 37.719 + x * 78.233 + y * 12.9898) * 43758.5453;
        
        double noise = nx - floor(nx) + ny - floor(ny) + nz - floor(nz);
        noise = noise - floor(noise);
        
        total += noise * amplitude;
        max_value += amplitude;
        amplitude *= persistence;
        frequency *= 2;
    }
    return total / max_value;
}

/* ---- 1.3 Differential Equation Solvers (RK4, Verlet) ---- */
typedef struct { Vec3 pos, vel; } ParticleState;
typedef Vec3 (*ForceFunc)(ParticleState, void*);

/* 4th Order Runge-Kutta integration for physics */
ParticleState rk4_integrate(ParticleState state, double dt, ForceFunc force, void* data) {
    ParticleState k1, k2, k3, k4, temp;
    
    // k1 = f(t, y)
    Vec3 a1 = force(state, data);
    k1.pos = state.vel;
    k1.vel = a1;
    
    // k2 = f(t + dt/2, y + dt*k1/2)
    temp.pos = vec3_add(state.pos, vec3_mul(k1.pos, dt/2));
    temp.vel = vec3_add(state.vel, vec3_mul(k1.vel, dt/2));
    Vec3 a2 = force(temp, data);
    k2.pos = temp.vel;
    k2.vel = a2;
    
    // k3 = f(t + dt/2, y + dt*k2/2)
    temp.pos = vec3_add(state.pos, vec3_mul(k2.pos, dt/2));
    temp.vel = vec3_add(state.vel, vec3_mul(k2.vel, dt/2));
    Vec3 a3 = force(temp, data);
    k3.pos = temp.vel;
    k3.vel = a3;
    
    // k4 = f(t + dt, y + dt*k3)
    temp.pos = vec3_add(state.pos, vec3_mul(k3.pos, dt));
    temp.vel = vec3_add(state.vel, vec3_mul(k3.vel, dt));
    Vec3 a4 = force(temp, data);
    k4.pos = temp.vel;
    k4.vel = a4;
    
    // Final integration
    ParticleState result;
    result.pos = vec3_add(state.pos, vec3_mul(vec3_add(vec3_add(k1.pos, vec3_mul(k2.pos, 2)), 
                       vec3_add(vec3_mul(k3.pos, 2), k4.pos)), dt/6));
    result.vel = vec3_add(state.vel, vec3_mul(vec3_add(vec3_add(k1.vel, vec3_mul(k2.vel, 2)), 
                       vec3_add(vec3_mul(k3.vel, 2), k4.vel)), dt/6));
    return result;
}

/* ---- 1.4 Chaotic Systems (Lorenz Attractor, Double Pendulum) ---- */
typedef struct {
    double x, y, z;
    double sigma, rho, beta;
} LorenzSystem;

void lorenz_step(LorenzSystem* ls, double dt) {
    double dx = ls->sigma * (ls->y - ls->x);
    double dy = ls->x * (ls->rho - ls->z) - ls->y;
    double dz = ls->x * ls->y - ls->beta * ls->z;
    
    ls->x += dx * dt;
    ls->y += dy * dt;
    ls->z += dz * dt;
}

/* === 2. PERSISTENT STATE SYSTEM (Database-like with Versioning) === */

#define MAX_UNIVERSES 16
#define MAX_OBJECTS_PER_UNIVERSE 1024
#define MAX_NAME_LEN 64
#define STATE_VERSION 2

#pragma pack(push, 1)
typedef struct {
    char id[32];
    char name[MAX_NAME_LEN];
    int32_t type;                // 0=avatar, 1=item, 2=portal, 3=field
    Vec3 position;
    Vec3 velocity;
    Quat rotation;
    double scale;
    double energy;
    double creation_time;
    double last_update;
    uint32_t flags;
    char metadata[128];          // JSON-like metadata string
} PersistentObject;

typedef struct {
    int32_t version;
    char universe_id[32];
    char name[MAX_NAME_LEN];
    double creation_timestamp;
    double last_save_timestamp;
    uint32_t object_count;
    PersistentObject objects[MAX_OBJECTS_PER_UNIVERSE];
    double physics_constants[8]; // G, c, etc.
    char save_checksum[32];
} UniverseState;

typedef struct {
    int32_t multiverse_version;
    char multiverse_name[MAX_NAME_LEN];
    uint32_t universe_count;
    char universe_ids[MAX_UNIVERSES][32];
    double global_time;
    uint64_t total_transactions;
} MultiVerseState;
#pragma pack(pop)

/* State Management Functions */
void generate_checksum(UniverseState* state, char* checksum) {
    // Simple checksum for demo (in production use SHA256)
    unsigned long sum = 0;
    unsigned char* bytes = (unsigned char*)state;
    size_t len = sizeof(UniverseState) - 32; // Exclude checksum field itself
    
    for(size_t i=0; i<len; i++) {
        sum = (sum * 31 + bytes[i]) % 1000000007;
    }
    snprintf(checksum, 32, "%08lx", sum);
}

int save_universe_state(const char* filename, UniverseState* state) {
    generate_checksum(state, state->save_checksum);
    state->last_save_timestamp = (double)time(NULL);
    
    FILE* f = fopen(filename, "wb");
    if(!f) return 0;
    
    size_t written = fwrite(state, sizeof(UniverseState), 1, f);
    fclose(f);
    return written == 1;
}

int load_universe_state(const char* filename, UniverseState* state) {
    FILE* f = fopen(filename, "rb");
    if(!f) return 0;
    
    size_t read = fread(state, sizeof(UniverseState), 1, f);
    fclose(f);
    
    if(read != 1) return 0;
    
    // Verify checksum
    char calculated[32];
    generate_checksum(state, calculated);
    return strncmp(state->save_checksum, calculated, 32) == 0;
}

/* === 3. NETWORK LAYER (Multi-threaded Client-Server with Protocol Buffers) === */

#define NET_PORT 8080
#define MAX_CLIENTS 32
#define NET_BUFFER_SIZE 8192
#define PROTOCOL_VERSION 1

typedef enum {
    MSG_CONNECT = 1,
    MSG_DISCONNECT = 2,
    MSG_OBJECT_UPDATE = 3,
    MSG_UNIVERSE_STATE = 4,
    MSG_CHAT = 5,
    MSG_COMMAND = 6
} MessageType;

#pragma pack(push, 1)
typedef struct {
    uint32_t type;
    uint32_t size;
    uint64_t timestamp;
    uint32_t crc;
} MessageHeader;

typedef struct {
    MessageHeader header;
    char client_id[32];
    char universe_id[32];
} ConnectMessage;

typedef struct {
    MessageHeader header;
    PersistentObject object;
    double server_timestamp;
} ObjectUpdateMessage;

typedef struct {
    MessageHeader header;
    UniverseState state;
} UniverseStateMessage;
#pragma pack(pop)

/* Network Server Structure */
typedef struct {
    int socket;
    struct sockaddr_in address;
    pthread_t thread;
    atomic_bool running;
    int client_sockets[MAX_CLIENTS];
    pthread_mutex_t client_mutex;
} NetworkServer;

/* Network Client Structure */
typedef struct {
    int socket;
    struct sockaddr_in server_addr;
    pthread_t receive_thread;
    atomic_bool connected;
    char client_id[32];
} NetworkClient;

uint32_t calculate_crc(const void* data, size_t length) {
    // Simple CRC for demo (use proper CRC32 in production)
    const unsigned char* bytes = (const unsigned char*)data;
    uint32_t crc = 0xFFFFFFFF;
    for(size_t i=0; i<length; i++) {
        crc ^= bytes[i];
        for(int j=0; j<8; j++) {
            crc = (crc >> 1) ^ (0xEDB88320 & -(crc & 1));
        }
    }
    return ~crc;
}

/* Simple Non-blocking Socket Setup */
int set_nonblocking(int sock) {
    int flags = fcntl(sock, F_GETFL, 0);
    if(flags == -1) return -1;
    return fcntl(sock, F_SETFL, flags | O_NONBLOCK);
}

/* === 4. VISUALIZATION ENGINE (SDL-based 3D Renderer) === */

#ifdef WITH_SDL
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

typedef struct {
    SDL_Window* window;
    SDL_GLContext gl_context;
    int width, height;
    atomic_bool should_close;
    pthread_t render_thread;
    
    // Camera
    Vec3 camera_pos;
    Vec3 camera_front;
    Vec3 camera_up;
    float yaw, pitch;
    
    // Projection
    Mat4 projection;
    Mat4 view;
    
    // Render data
    GLuint shader_program;
    GLuint vao, vbo;
    GLuint texture_id;
} Visualizer;

/* Simple OpenGL Shaders (embedded as strings) */
const char* vertex_shader_source = 
    "#version 330 core\n"
    "layout(location = 0) in vec3 aPos;\n"
    "layout(location = 1) in vec3 aColor;\n"
    "uniform mat4 projection;\n"
    "uniform mat4 view;\n"
    "uniform mat4 model;\n"
    "out vec3 ourColor;\n"
    "void main() {\n"
    "    gl_Position = projection * view * model * vec4(aPos, 1.0);\n"
    "    ourColor = aColor;\n"
    "}\n";

const char* fragment_shader_source = 
    "#version 330 core\n"
    "in vec3 ourColor;\n"
    "out vec4 FragColor;\n"
    "void main() {\n"
    "    FragColor = vec4(ourColor, 1.0);\n"
    "}\n";

GLuint compile_shader(GLenum type, const char* source) {
    GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &source, NULL);
    glCompileShader(shader);
    
    GLint success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if(!success) {
        char info_log[512];
        glGetShaderInfoLog(shader, 512, NULL, info_log);
        fprintf(stderr, "Shader compilation failed: %s\n", info_log);
    }
    return shader;
}

GLuint create_shader_program() {
    GLuint vertex_shader = compile_shader(GL_VERTEX_SHADER, vertex_shader_source);
    GLuint fragment_shader = compile_shader(GL_FRAGMENT_SHADER, fragment_shader_source);
    
    GLuint program = glCreateProgram();
    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);
    glLinkProgram(program);
    
    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);
    
    return program;
}

void render_universe(Visualizer* vis, UniverseState* state) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glUseProgram(vis->shader_program);
    
    // Update camera
    vis->camera_front.x = cos(vis->yaw) * cos(vis->pitch);
    vis->camera_front.y = sin(vis->pitch);
    vis->camera_front.z = sin(vis->yaw) * cos(vis->pitch);
    vec3_normalize(&vis->camera_front);
    
    // Calculate view matrix
    Vec3 target = vec3_add(vis->camera_pos, vis->camera_front);
    vis->view = mat4_look_at(vis->camera_pos, target, vis->camera_up);
    
    // Set uniforms
    GLuint proj_loc = glGetUniformLocation(vis->shader_program, "projection");
    GLuint view_loc = glGetUniformLocation(vis->shader_program, "view");
    glUniformMatrix4fv(proj_loc, 1, GL_FALSE, (float*)&vis->projection);
    glUniformMatrix4fv(view_loc, 1, GL_FALSE, (float*)&vis->view);
    
    // Render objects as colored cubes
    for(uint32_t i=0; i<state->object_count; i++) {
        PersistentObject* obj = &state->objects[i];
        
        Mat4 model = mat4_identity();
        model = mat4_translate(model, obj->position);
        
        // Set model matrix uniform
        GLuint model_loc = glGetUniformLocation(vis->shader_program, "model");
        glUniformMatrix4fv(model_loc, 1, GL_FALSE, (float*)&model);
        
        // Draw cube (simplified)
        glDrawArrays(GL_TRIANGLES, 0, 36);
    }
    
    SDL_GL_SwapWindow(vis->window);
}
#endif

/* === 5. MULTIVERSE CORE ARCHITECTURE === */

typedef struct {
    char id[32];
    char name[MAX_NAME_LEN];
    UniverseState state;
    pthread_mutex_t state_mutex;
    
    // Physics thread
    pthread_t physics_thread;
    atomic_bool physics_running;
    double physics_time_step;
    
    // Network
    NetworkServer* net_server;
    
    // Visualization
#ifdef WITH_SDL
    Visualizer* visualizer;
#endif
    
    // Chaotic systems for procedural generation
    LorenzSystem lorenz;
    double fractal_seed;
} UniverseInstance;

typedef struct {
    char name[MAX_NAME_LEN];
    UniverseInstance universes[MAX_UNIVERSES];
    int universe_count;
    
    // Global state
    MultiVerseState global_state;
    pthread_mutex_t global_mutex;
    
    // Thread management
    pthread_t net_thread;
    atomic_bool running;
    
    // Portals between universes (inter-universe connections)
    struct {
        char from_universe[32];
        char to_universe[32];
        Vec3 from_position;
        Vec3 to_position;
        double cooldown;
    } portals[16];
    int portal_count;
} MultiVerse;

/* Universe Physics Thread Function */
void* universe_physics_thread(void* arg) {
    UniverseInstance* universe = (UniverseInstance*)arg;
    
    struct timespec ts;
    double accumulator = 0;
    double current_time = (double)clock() / CLOCKS_PER_SEC;
    
    while(universe->physics_running) {
        double new_time = (double)clock() / CLOCKS_PER_SEC;
        double frame_time = new_time - current_time;
        current_time = new_time;
        
        // Clamp frame time
        if(frame_time > 0.25) frame_time = 0.25;
        accumulator += frame_time;
        
        pthread_mutex_lock(&universe->state_mutex);
        
        while(accumulator >= universe->physics_time_step) {
            // Update chaotic systems
            lorenz_step(&universe->lorenz, universe->physics_time_step);
            
            // Update fractal seed
            universe->fractal_seed = fmod(universe->fractal_seed * 1.6180339887, 1.0);
            
            // Update objects with RK4 physics
            for(uint32_t i=0; i<universe->state.object_count; i++) {
                PersistentObject* obj = &universe->state.objects[i];
                
                // Define force function (gravity toward origin)
                Vec3 gravity_force(ParticleState state, void* data) {
                    double G = -9.81;
                    double distance = vec3_length(state.pos);
                    if(distance > 0.001) {
                        Vec3 direction = vec3_mul(state.pos, -1.0/distance);
                        return vec3_mul(direction, G);
                    }
                    return (Vec3){0,0,0};
                }
                
                // Apply physics
                ParticleState pstate = {obj->position, obj->velocity};
                ParticleState new_state = rk4_integrate(pstate, universe->physics_time_step, 
                                                       gravity_force, NULL);
                
                obj->position = new_state.pos;
                obj->velocity = new_state.vel;
                obj->last_update = current_time;
                
                // Apply fractal terrain if object is near "ground"
                if(obj->position.y < 0) {
                    double noise = fractal_noise_3d(
                        obj->position.x, obj->position.z, universe->fractal_seed, 4, 0.5
                    );
                    obj->position.y = noise * 10 - 5;
                }
            }
            
            accumulator -= universe->physics_time_step;
            universe->state.last_save_timestamp = current_time;
        }
        
        pthread_mutex_unlock(&universe->state_mutex);
        
        // Sleep to prevent CPU hogging
        ts.tv_sec = 0;
        ts.tv_nsec = 1000000; // 1ms
        nanosleep(&ts, NULL);
    }
    
    return NULL;
}

/* Network Thread Function */
void* network_server_thread(void* arg) {
    NetworkServer* server = (NetworkServer*)arg;
    
    // Create socket
    server->socket = socket(AF_INET, SOCK_STREAM, 0);
    if(server->socket < 0) {
        perror("Socket creation failed");
        return NULL;
    }
    
    // Set socket options
    int opt = 1;
    setsockopt(server->socket, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt));
    set_nonblocking(server->socket);
    
    // Bind
    server->address.sin_family = AF_INET;
    server->address.sin_addr.s_addr = INADDR_ANY;
    server->address.sin_port = htons(NET_PORT);
    
    if(bind(server->socket, (struct sockaddr*)&server->address, sizeof(server->address)) < 0) {
        perror("Bind failed");
        return NULL;
    }
    
    // Listen
    listen(server->socket, MAX_CLIENTS);
    printf("Server listening on port %d\n", NET_PORT);
    
    fd_set read_fds;
    struct timeval tv = {0, 100000}; // 100ms timeout
    
    while(server->running) {
        FD_ZERO(&read_fds);
        FD_SET(server->socket, &read_fds);
        int max_fd = server->socket;
        
        pthread_mutex_lock(&server->client_mutex);
        for(int i=0; i<MAX_CLIENTS; i++) {
            if(server->client_sockets[i] > 0) {
                FD_SET(server->client_sockets[i], &read_fds);
                if(server->client_sockets[i] > max_fd) {
                    max_fd = server->client_sockets[i];
                }
            }
        }
        pthread_mutex_unlock(&server->client_mutex);
        
        int activity = select(max_fd+1, &read_fds, NULL, NULL, &tv);
        
        if(activity < 0 && errno != EINTR) {
            perror("Select error");
        }
        
        // New connection
        if(FD_ISSET(server->socket, &read_fds)) {
            struct sockaddr_in client_addr;
            socklen_t addr_len = sizeof(client_addr);
            int new_socket = accept(server->socket, (struct sockaddr*)&client_addr, &addr_len);
            
            if(new_socket >= 0) {
                set_nonblocking(new_socket);
                
                pthread_mutex_lock(&server->client_mutex);
                for(int i=0; i<MAX_CLIENTS; i++) {
                    if(server->client_sockets[i] == 0) {
                        server->client_sockets[i] = new_socket;
                        printf("New client connected: %s:%d\n", 
                               inet_ntoa(client_addr.sin_addr), ntohs(client_addr.sin_port));
                        break;
                    }
                }
                pthread_mutex_unlock(&server->client_mutex);
            }
        }
        
        // Check client sockets for data
        pthread_mutex_lock(&server->client_mutex);
        for(int i=0; i<MAX_CLIENTS; i++) {
            int client_sock = server->client_sockets[i];
            if(client_sock > 0 && FD_ISSET(client_sock, &read_fds)) {
                char buffer[NET_BUFFER_SIZE];
                int valread = read(client_sock, buffer, NET_BUFFER_SIZE);
                
                if(valread == 0) {
                    // Client disconnected
                    getpeername(client_sock, (struct sockaddr*)&server->address, 
                                (socklen_t*)&addr_len);
                    printf("Client disconnected\n");
                    close(client_sock);
                    server->client_sockets[i] = 0;
                } else if(valread > 0) {
                    // Process message (simplified)
                    MessageHeader* header = (MessageHeader*)buffer;
                    if(header->type == MSG_OBJECT_UPDATE) {
                        // Echo back for now
                        send(client_sock, buffer, valread, 0);
                    }
                }
            }
        }
        pthread_mutex_unlock(&server->client_mutex);
        
        usleep(10000); // 10ms
    }
    
    return NULL;
}

/* === 6. MAIN MULTIVERSE INTEGRATION === */

void init_multiverse(MultiVerse* mv, const char* name) {
    strncpy(mv->name, name, MAX_NAME_LEN);
    mv->universe_count = 0;
    mv->portal_count = 0;
    mv->running = true;
    
    pthread_mutex_init(&mv->global_mutex, NULL);
    
    // Initialize global state
    mv->global_state.multiverse_version = PROTOCOL_VERSION;
    strncpy(mv->global_state.multiverse_name, name, MAX_NAME_LEN);
    mv->global_state.global_time = 0;
    mv->global_state.total_transactions = 0;
}

void add_universe(MultiVerse* mv, const char* id, const char* name) {
    if(mv->universe_count >= MAX_UNIVERSES) return;
    
    UniverseInstance* uni = &mv->universes[mv->universe_count];
    strncpy(uni->id, id, 32);
    strncpy(uni->name, name, MAX_NAME_LEN);
    
    // Initialize universe state
    uni->state.version = STATE_VERSION;
    strncpy(uni->state.universe_id, id, 32);
    strncpy(uni->state.name, name, MAX_NAME_LEN);
    uni->state.creation_timestamp = (double)time(NULL);
    uni->state.object_count = 0;
    
    // Initialize physics constants
    uni->state.physics_constants[0] = 6.67430e-11; // G
    uni->state.physics_constants[1] = 299792458.0; // c
    uni->state.physics_constants[2] = 1.0;         // time scale
    
    // Initialize chaotic systems
    uni->lorenz = (LorenzSystem){0.1, 0.1, 0.1, 10.0, 28.0, 8.0/3.0};
    uni->fractal_seed = (double)rand() / RAND_MAX;
    
    // Initialize mutex
    pthread_mutex_init(&uni->state_mutex, NULL);
    
    // Start physics thread
    uni->physics_time_step = 0.01;
    uni->physics_running = true;
    pthread_create(&uni->physics_thread, NULL, universe_physics_thread, uni);
    
    mv->universe_count++;
    strncpy(mv->global_state.universe_ids[mv->global_state.universe_count++], id, 32);
}

/* === 7. COMMAND LINE INTERFACE & MAIN FUNCTION === */

void print_usage() {
    printf("MultiVerse Simulation System\n");
    printf("Usage:\n");
    printf("  Server mode: multiverse --server --port PORT --universe-count N\n");
    printf("  Client mode: multiverse --client --host HOST --port PORT\n");
    printf("  Standalone:  multiverse --standalone --universes 3\n");
    printf("\nOptions:\n");
    printf("  --save FILE          Save state to file\n");
    printf("  --load FILE          Load state from file\n");
    printf("  --with-visualizer    Enable 3D visualization (requires SDL2)\n");
}

int main(int argc, char** argv) {
    srand(time(NULL));
    
    bool server_mode = false;
    bool client_mode = false;
    bool standalone = false;
    bool with_visualizer = false;
    int port = NET_PORT;
    int universe_count = 1;
    char* save_file = NULL;
    char* load_file = NULL;
    char* host = "127.0.0.1";
    
    // Parse arguments
    for(int i=1; i<argc; i++) {
        if(strcmp(argv[i], "--server") == 0) server_mode = true;
        else if(strcmp(argv[i], "--client") == 0) client_mode = true;
        else if(strcmp(argv[i], "--standalone") == 0) standalone = true;
        else if(strcmp(argv[i], "--port") == 0 && i+1 < argc) port = atoi(argv[++i]);
        else if(strcmp(argv[i], "--host") == 0 && i+1 < argc) host = argv[++i];
        else if(strcmp(argv[i], "--universe-count") == 0 && i+1 < argc) universe_count = atoi(argv[++i]);
        else if(strcmp(argv[i], "--save") == 0 && i+1 < argc) save_file = argv[++i];
        else if(strcmp(argv[i], "--load") == 0 && i+1 < argc) load_file = argv[++i];
        else if(strcmp(argv[i], "--with-visualizer") == 0) with_visualizer = true;
        else if(strcmp(argv[i], "--help") == 0) {
            print_usage();
            return 0;
        }
    }
    
    // Initialize MultiVerse
    MultiVerse mv;
    init_multiverse(&mv, "MyMultiVerse");
    
    // Create universes
    for(int i=0; i<universe_count; i++) {
        char id[32], name[32];
        snprintf(id, 32, "universe_%d", i);
        snprintf(name, 32, "Universe %d", i);
        add_universe(&mv, id, name);
    }
    
    // Add some test objects to first universe
    if(mv.universe_count > 0) {
        UniverseInstance* uni = &mv.universes[0];
        pthread_mutex_lock(&uni->state_mutex);
        
        for(int i=0; i<10; i++) {
            if(uni->state.object_count < MAX_OBJECTS_PER_UNIVERSE) {
                PersistentObject* obj = &uni->state.objects[uni->state.object_count++];
                snprintf(obj->id, 32, "obj_%d", i);
                snprintf(obj->name, MAX_NAME_LEN, "Object %d", i);
                obj->type = i % 3;
                obj->position = (Vec3){rand()%100-50, rand()%100, rand()%100-50};
                obj->velocity = (Vec3){rand()%10-5, rand()%10-5, rand()%10-5};
                obj->rotation = (Quat){1,0,0,0};
                obj->scale = 1.0;
                obj->energy = 100.0;
                obj->creation_time = (double)time(NULL);
                obj->last_update = obj->creation_time;
                snprintf(obj->metadata, 128, "{\"color\":\"#%06x\"}", rand()%0xFFFFFF);
            }
        }
        
        pthread_mutex_unlock(&uni->state_mutex);
    }
    
    // Save state if requested
    if(save_file && mv.universe_count > 0) {
        save_universe_state(save_file, &mv.universes[0].state);
        printf("Saved state to %s\n", save_file);
    }
    
    // Load state if requested
    if(load_file && mv.universe_count > 0) {
        if(load_universe_state(load_file, &mv.universes[0].state)) {
            printf("Loaded state from %s\n", load_file);
        } else {
            printf("Failed to load state from %s\n", load_file);
        }
    }
    
    // Run in appropriate mode
    if(server_mode) {
        printf("Running in server mode on port %d\n", port);
        
        NetworkServer server;
        memset(&server, 0, sizeof(server));
        server.running = true;
        pthread_mutex_init(&server.client_mutex, NULL);
        
        pthread_create(&server.thread, NULL, network_server_thread, &server);
        
        // Server main loop
        while(mv.running) {
            // Update global time
            mv.global_state.global_time += 0.1;
            
            // Save periodic snapshots
            static double last_save = 0;
            double now = (double)clock() / CLOCKS_PER_SEC;
            if(now - last_save > 60.0) { // Every minute
                for(int i=0; i<mv.universe_count; i++) {
                    char filename[64];
                    snprintf(filename, 64, "universe_%s_snapshot.bin", mv.universes[i].id);
                    pthread_mutex_lock(&mv.universes[i].state_mutex);
                    save_universe_state(filename, &mv.universes[i].state);
                    pthread_mutex_unlock(&mv.universes[i].state_mutex);
                }
                last_save = now;
            }
            
            usleep(100000); // 100ms
        }
        
        server.running = false;
        pthread_join(server.thread, NULL);
    }
    else if(client_mode) {
        printf("Running in client mode connecting to %s:%d\n", host, port);
        
        // Simple client implementation
        int sock = socket(AF_INET, SOCK_STREAM, 0);
        struct sockaddr_in serv_addr;
        serv_addr.sin_family = AF_INET;
        serv_addr.sin_port = htons(port);
        inet_pton(AF_INET, host, &serv_addr.sin_addr);
        
        if(connect(sock, (struct sockaddr*)&serv_addr, sizeof(serv_addr)) < 0) {
            perror("Connection failed");
            return 1;
        }
        
        printf("Connected to server\n");
        
        // Send connect message
        ConnectMessage conn_msg;
        conn_msg.header.type = MSG_CONNECT;
        conn_msg.header.size = sizeof(ConnectMessage);
        conn_msg.header.timestamp = time(NULL);
        strncpy(conn_msg.client_id, "test_client", 32);
        strncpy(conn_msg.universe_id, "universe_0", 32);
        conn_msg.header.crc = calculate_crc(&conn_msg, sizeof(ConnectMessage) - sizeof(MessageHeader));
        
        send(sock, &conn_msg, sizeof(conn_msg), 0);
        
        // Client main loop
        char buffer[NET_BUFFER_SIZE];
        while(mv.running) {
            int valread = read(sock, buffer, NET_BUFFER_SIZE);
            if(valread > 0) {
                MessageHeader* header = (MessageHeader*)buffer;
                printf("Received message type: %u, size: %u\n", header->type, header->size);
            }
            
            usleep(100000); // 100ms
        }
        
        close(sock);
    }
    else {
        printf("Running in standalone mode with %d universe(s)\n", mv.universe_count);
        printf("Press Ctrl+C to exit\n");
        
        // Standalone main loop
        while(mv.running) {
            // Display universe stats
            static int counter = 0;
            if(counter++ % 10 == 0) {
                printf("\n=== MultiVerse Status ===\n");
                printf("Global Time: %.2f\n", mv.global_state.global_time);
                printf("Total Transactions: %lu\n", mv.global_state.total_transactions);
                
                for(int i=0; i<mv.universe_count; i++) {
                    pthread_mutex_lock(&mv.universes[i].state_mutex);
                    printf("Universe %s: %d objects, last update: %.2f\n", 
                           mv.universes[i].id, 
                           mv.universes[i].state.object_count,
                           mv.universes[i].state.last_save_timestamp);
                    pthread_mutex_unlock(&mv.universes[i].state_mutex);
                }
            }
            
            // Fractal math demonstration
            if(mv.universe_count > 0) {
                UniverseInstance* uni = &mv.universes[0];
                pthread_mutex_lock(&uni->state_mutex);
                
                // Generate fractal terrain for demonstration
                double x = sin(mv.global_state.global_time * 0.1);
                double z = cos(mv.global_state.global_time * 0.1);
                double terrain_height = fractal_noise_3d(x, z, uni->fractal_seed, 4, 0.5) * 20;
                
                if(uni->state.object_count > 0) {
                    // Update first object with fractal-based position
                    uni->state.objects[0].position.y = terrain_height;
                    
                    // Calculate Mandelbrot set value for color metadata
                    Complex c = {x * 2 - 1, z * 2 - 1};
                    int mandel_val = mandelbrot_iter(c, 100);
                    snprintf(uni->state.objects[0].metadata, 128, 
                            "{\"fractal_value\":%d,\"terrain_height\":%.2f}", 
                            mandel_val, terrain_height);
                }
                
                pthread_mutex_unlock(&uni->state_mutex);
            }
            
            usleep(100000); // 100ms
            mv.global_state.global_time += 0.1;
        }
    }
    
    // Cleanup
    for(int i=0; i<mv.universe_count; i++) {
        mv.universes[i].physics_running = false;
        pthread_join(mv.universes[i].physics_thread, NULL);
        pthread_mutex_destroy(&mv.universes[i].state_mutex);
    }
    
    pthread_mutex_destroy(&mv.global_mutex);
    
    printf("MultiVerse shutdown complete\n");
    return 0;
}