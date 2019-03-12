/*
Arttu Hyvönen 3/2019
TODO:
Only updating changed regions in rendering
Sort z-axis
Barnes-Hut algorithm... Maybe in another version
file inut/output

*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <SDL2/SDL.h>
#include <ctime>
#include <algorithm>
#include <vector>
#include <random>

using std::cout;
using std::cin;
using std::endl;

/*
Units:
distance - 1e8 meters
mass     - 1e18 kilograms
time     - 1e3 seconds
*/

//Simulation constants
const double G = 6.674E-11;             // Gravitational constant real:-11
const double pi = 3.1416;
const int n = 1000;                     // Number of particles
const int d = 3;                        // Number of dimensions
const double crash = 4;                 // Min distance between particles
double dt = 1;                          // Time step in time units                     
const double scale = 1;                 // Size of pixel in distance units
const int mass_scale = 5;               // Masses range 1e(18+s)-1e(20+s)
const double start_speed = 0.2;        // Multiplier for initial speeds
const double extermination_zone = 6000; // Place where particles die
const double pos_dist_dev_z = 1.0;        // Deviation of particle position distribution
const double pos_dist_dev_xy = 1.0;        // Deviation of particle position distribution
const double vel_dist_dev = 5.0;        // Deviation of particle velocity distribution
const double rotation_bias = 0;       // Rotation in x-y plane
double system_mass = 0;                 // Total mass of the system

double mr[d] = {};                           // Center of mass
double vr[d] = {};                           // Velocity of center of mass

const int steps = 0;                    // Limit amount of steps to be taken. 0 = no limit
const int part_size = 8;                // Size of particles on screen

//SDL stuff
SDL_Window* gWindow = NULL;
SDL_Renderer* gRenderer = NULL;

bool screen = true;

const int SCREEN_WIDTH = 1000;
const int SCREEN_HEIGHT = 1000;

//Log
bool sim_log = true;
double max_mass = 0;
double cm_vel = 0;

//Line properties
bool line = false;
const int line_len = 100;
const int line_res = 20;
int line_point = 0;
int line_array[n][line_len][2];

//BH
const int d2 = 8;           // BH sub nodes
const int cases[8][3] = {{1,1,1},{1,1,-1},{1,-1,1},{1,-1,-1},
                        {-1,1,1},{-1,1,-1},{-1,-1,1},{-1,-1,-1}};
const double theta = 0.5;   // 0 = brute force
int part_ind = 0;
const bool BH = false;

/*
Globals end, code begins
*/

struct Part {
    double pos [d];
    double vel [d];
    double mass;
    bool e = true;
} particles[n];

double force(double m1, double m2, int s) {
    return G*m1*m2/(s*s);
}

double distance(double p1[d], double p2[d]) {
    double s = 0;
    for(int j = 0; j < d; j++) s += (p2[j]-p1[j])*(p2[j]-p1[j]);
    return sqrt(s);
}

bool particle_sorter(Part const& p1, Part const& p2) {
    return p1.pos[d-1] < p2.pos[d-1];
}

bool init() {
    
    bool success = true;
    //srand (time(NULL));
    std::default_random_engine generator;
    std::normal_distribution<double> pos_dist_z(0, pos_dist_dev_z);
    std::normal_distribution<double> pos_dist_xy(0, pos_dist_dev_xy);
    std::normal_distribution<double> vel1_dist(0, vel_dist_dev);
    std::normal_distribution<double> vel2_dist(rotation_bias, vel_dist_dev);
    
    //Simulation init
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < d; j++) {
            if(j <= 2) {
                particles[i].pos[j] = 100*pos_dist_xy(generator)*scale;
                for(int k = 0; k < line_len; k++) line_array[i][k][j] = particles[i].pos[j];
            } else {
                particles[i].pos[j] = 100*pos_dist_z(generator)*scale;
            }
        }
        
        double angle = atan2(particles[i].pos[1],particles[i].pos[0]) + pi*0.5;
        double v1 = start_speed*vel1_dist(generator)/100;
        double v2 = start_speed*vel2_dist(generator)/100*(1.0/(1.0+2*rotation_bias));
        
        particles[i].vel[0] = cos(angle)*v1-sin(angle)*v2;
        particles[i].vel[1] = sin(angle)*v1+cos(angle)*v2;
        
        for(int j = 2; j < d; j++) particles[i].vel[j] = (1.0/(1.0+2*rotation_bias))*start_speed*vel1_dist(generator)/100;
        
        particles[i].mass = (rand()%100)*pow(10, mass_scale)+10;
        system_mass += particles[i].mass;
    }

    if(!screen) return success;
    
    //SDL init
    if(SDL_Init(SDL_INIT_VIDEO) < 0) {
        cout << "SDL could not initialize";
        success = false;
    } else {
        
        if(!SDL_SetHint( SDL_HINT_RENDER_SCALE_QUALITY, "1")) {
			cout << "Warning: Linear texture filtering not enabled!";
		}
        
        gWindow = SDL_CreateWindow("N-Body Simulation", 
                                    SDL_WINDOWPOS_UNDEFINED, 
                                    SDL_WINDOWPOS_UNDEFINED, 
                                    SCREEN_WIDTH, 
                                    SCREEN_HEIGHT, 
                                    SDL_WINDOW_SHOWN);
        if(gWindow == NULL) {
            cout << "Can't create window";
            success = false;
        } else {
            gRenderer = SDL_CreateRenderer( gWindow, -1, SDL_RENDERER_ACCELERATED );
			if( gRenderer == NULL )
			{
				cout << "Renderer could not be created";
				success = false;
			}
			else
			{
				SDL_SetRenderDrawColor( gRenderer, 0xFF, 0xFF, 0xFF, 0xFF );
			}
        }
    }
    
    return success;
}

void close() {
    SDL_DestroyRenderer( gRenderer );
    SDL_DestroyWindow(gWindow);
    gRenderer = NULL;
    gWindow = NULL;
    
    SDL_Quit();
}

//Regular n^2 update
void update() {
    
    // Position update
    for(int i = 0; i < n; i++) {
        
        if(!particles[i].e) continue;
        
        Part p1 = particles[i];
        double f[d] = {};        
        for(int k = 0; k < n; k++) {
            if(i==k || !particles[k].e) continue;
            Part p2 = particles[k];
            double r [d];
            double s = 0;
            for(int j = 0; j < d; j++) {
                r[j] = p2.pos[j]-p1.pos[j];
                s += r[j]*r[j];
            }
            s = sqrt(s);
            double c = force(p1.mass, p2.mass, s);
            for(int j = 0; j < d; j++) f[j] += c*r[j]/abs(s);
        }
        
        for(int j = 0; j < d; j++) {
            double a = f[j]/p1.mass;
            particles[i].vel[j] += dt*a;
            particles[i].pos[j] += dt*particles[i].vel[j];
            //limiting area
            //if(particles[i].pos[j] < 0) particles[i].pos[j] += 1000;
            //particles[i].pos[j] = fmod(particles[i].pos[j], 1000*scale);
        }
    }
    
    // Crash check
    
    for(int i = 0; i < n; i++) {
        
        if(!particles[i].e) continue;
        
        Part p1 = particles[i];
        for(int k = i+1; k < n; k++) {    
            if(!particles[k].e) continue;
            Part p2 = particles[k];
            double s = 0;
            for(int j = 0; j < d; j++) s += (p2.pos[j]-p1.pos[j])*(p2.pos[j]-p1.pos[j]);
            s = sqrt(s);
            if(s < crash) {
                double pm1[d];
                double pm2[d];
                double tot_m = p1.mass+p2.mass;
                for(int j = 0; j < d; j++) {
                    pm1[j] = p1.vel[j]*p1.mass;
                    pm2[j] = p2.vel[j]*p2.mass;
                }
                double vel3[d];
                double pos3[d];
                for(int j = 0; j < d; j++) {
                    vel3[j] = (pm1[j]+pm2[j])/(tot_m);
                    pos3[j] = (p1.mass*p1.pos[j]+p2.mass*p2.pos[j])/(tot_m);
                }
                Part p3;
                for(int j = 0; j < d; j++) {
                    p3.pos[j] = pos3[j];
                    p3.vel[j] = vel3[j];
                }
                p3.mass = tot_m;
                particles[i] = p3;
                particles[k].e = false;
            }    
        }
    }
}

class Node {
        Part p;
        Node * subNodes;
        bool part = false;
        bool nodes = false;
    public:
        double f[d] = {};
        double com[d];
        double mass = 0;
        double center [d];
        double side;
        Node(double c[d], double s);
        Node() {}
        void add_particle(Part pa);
        double* force_on_particle(Part pa);
};

Node::Node(double c[d], double s) {
    for(int i = 0; i < d; i++) center[i] = c[i];
    side=s;
}

void Node::add_particle(Part pa) {
    if(!part && !nodes) {
        p = pa;
        for(int i = 0; i < d; i++) com[i] = p.pos[i];
        mass = p.mass;
        part = true;
    } else if(part && !nodes) {
        double min_dist_p = 2*side;
        double min_dist_pa = 2*side;
        int ind_p;
        int ind_pa;
        subNodes = new Node [d2];
        for(int i = 0; i < d2; i++) {
            for(int j = 0; j < d; j++) subNodes[i].center[j] = center[j]+cases[i][j]*side*0.25;
            subNodes[i].side = side*0.5;
//            subNodes[i] = Node(c, side*0.5);
            double dist_p = distance(subNodes[i].center, p.pos);
            double dist_pa = distance(subNodes[i].center, pa.pos);
            if(dist_p < min_dist_p) {
                min_dist_p = dist_p;
                ind_p = i;
            }
            if(dist_pa < min_dist_pa) {
                min_dist_pa = dist_pa;
                ind_pa = i;
            }
        }
        subNodes[ind_pa].add_particle(pa);
        subNodes[ind_p].add_particle(p);
        part = false;
        nodes = true;
        for(int i = 0; i < d2; i++) {
            mass += subNodes[i].mass;
        }
        for(int i = 0; i < d2; i++) {
            if(subNodes[i].mass == 0) continue;
            for(int j = 0; j < d; j++) com[j] += subNodes[i].com[j]*subNodes[i].mass/mass;
        }
    } else if(!part && nodes){
        double min_dist_pa = side;
        int ind_pa;
        for(int i = 0; i < d2; i++) {
            double c[d];
            for(int j = 0; j < d; j++) c[j] = subNodes[i].center[j];
            double dist_pa = distance(c, pa.pos);
            if(dist_pa < side/4) {
                ind_pa = i;
                break;
            }
            if(dist_pa < min_dist_pa) {
                min_dist_pa = dist_pa;
                ind_pa = i;
            }
        }
        subNodes[ind_pa].add_particle(pa);
    }
}

double* Node::force_on_particle(Part pa) {
    for(int i = 0; i < d; i++) f[i] = 0;
    if(mass == 0) {
        return f;
    } else if(part) {
        double r [d];
        double s = 0;
        for(int j = 0; j < d; j++) {
            r[j] = p.pos[j]-pa.pos[j];
            s += r[j]*r[j];
        }
        s = sqrt(s);
        if(s == 0) return f;
        double c = force(pa.mass, p.mass, s);
        for(int j = 0; j < d; j++) f[j] += c*r[j]/abs(s);
        return f;
    } else if(side/distance(com, pa.pos) < theta) {
        double r [d];
        double s = 0;
        for(int j = 0; j < d; j++) {
            r[j] = com[j]-pa.pos[j];
            s += r[j]*r[j];
        }
        s = sqrt(s);
        double c = force(pa.mass, mass, s);
        for(int j = 0; j < d; j++) f[j] += c*r[j]/abs(s);
        return f;
    } else {
        for(int i = 0; i < d2; i++) {
            double* fs;
            fs = subNodes[i].force_on_particle(pa);
            for(int j = 0; j < d; j++) f[j] += *(fs+j);
        }
        return f;
    }
}

//Barnes-Hut nlog(n) update
void BHupdate() {
    double c[d] = {};
    Node top (c, extermination_zone);
    for(int i = 0; i < n; i++) {
        if(!particles[i].e) continue;
        top.add_particle(particles[i]);
    }
    for(int i = 0; i < n; i++) {
        if(!particles[i].e) continue;
        double f[d] = {};
        top.force_on_particle(particles[i]);
        for(int j = 0; j < d; j++) f[j] = top.f[j];
        for(int j = 0; j < d; j++) {
            double a = f[j]/particles[i].mass;
            particles[i].vel[j] += dt*a;
            particles[i].pos[j] += dt*particles[i].vel[j];
        }
    }
    
    
    for(int i = 0; i < n; i++) {
        if(!particles[i].e) continue;
        Part p1 = particles[i];
        for(int k = i+1; k < n; k++) {    
            if(!particles[k].e) continue;
            Part p2 = particles[k];
            double s = 0;
            for(int j = 0; j < d; j++) s += (p2.pos[j]-p1.pos[j])*(p2.pos[j]-p1.pos[j]);
            s = sqrt(s);
            if(s < crash) {
                double pm1[d];
                double pm2[d];
                double tot_m = p1.mass+p2.mass;
                for(int j = 0; j < d; j++) {
                    pm1[j] = p1.vel[j]*p1.mass;
                    pm2[j] = p2.vel[j]*p2.mass;
                }
                double vel3[d];
                double pos3[d];
                for(int j = 0; j < d; j++) {
                    vel3[j] = (pm1[j]+pm2[j])/(tot_m);
                    pos3[j] = (p1.mass*p1.pos[j]+p2.mass*p2.pos[j])/(tot_m);
                }
                Part p3;
                for(int j = 0; j < d; j++) {
                    p3.pos[j] = pos3[j];
                    p3.vel[j] = vel3[j];
                }
                p3.mass = tot_m;
                particles[i] = p3;
                particles[k].e = false;
            }    
        }
    }
}

void render(int step) {
    
    //Clear screen
    SDL_SetRenderDrawColor(gRenderer, 0, 0, 0, 255);
    SDL_RenderClear(gRenderer);
    
    //System center of mass
    for(int i = 0; i < d; i++) vr[i] = mr[i];
    for(int i = 0; i < d; i++) mr[i] = 0;
    for(int i = 0; i < n; i++) {
        if(!particles[i].e) continue;
        for(int j = 0; j < d; j++) mr[j] += particles[i].mass*particles[i].pos[j]/system_mass;
    }
    
    for(int i = 0; i < d; i++) vr[i] = mr[i]-vr[i];
    
    // Delete stray particles
    for(int i = 0; i < n; i++) {
        if(!particles[i].e) continue;
        double s = 0;
        for(int j = 0; j < d; j++) {
            double r = mr[j]-particles[i].pos[j];
            s += r*r;
        }
        if (sqrt(s) > extermination_zone) { 
            particles[i].e = false;
            system_mass -= particles[i].mass;
        }
    }
    
    //Draw
    //Grid
    SDL_SetRenderDrawColor(gRenderer, 100, 100, 100, 255);
    for(int i = 0; i < 12; i++) {
        SDL_RenderDrawLine(gRenderer, i*100-fmod(mr[0],100), 0, i*100-fmod(mr[0],100), 1000);
        SDL_RenderDrawLine(gRenderer, 1000, i*100-fmod(mr[1],100), 0, i*100-fmod(mr[1],100));
    }
    
    //Line
    SDL_SetRenderDrawColor(gRenderer, 100, 100, 150, 255);
    if(line && step % line_res == 0) {
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < 2; j++) line_array[i][line_point][j] = particles[i].pos[j];
        }
        line_point = (line_point + 1)%line_len;
    }
    
    if(line) {
        for(int i = 0; i < n; i++) {
            if(!particles[i].e) continue;
            
            for(int j = (line_point+1)%line_len; (j+1)%line_len != line_point; j = (j+1)%line_len) {
                SDL_RenderDrawLine(gRenderer, line_array[i][j][0]-(mr[0]-500), 
                                              line_array[i][j][1]-(mr[1]-500), 
                                              line_array[i][(j+1)%line_len][0]-(mr[0]-500),
                                              line_array[i][(j+1)%line_len][1]-(mr[1]-500));
            }
        }
    }
        
    //Particles
    for(Part p : particles) {
        if(!p.e) continue;
        
        float col = log10(p.mass*n/system_mass);
        col = ((col-1)/(1+abs(2*(col-1))) + 0.5);
        
        float size = 2*col*part_size+4;
        if(d > 2) size += (mr[2]-p.pos[2])*0.005;
        
        SDL_SetRenderDrawColor(gRenderer, 255, int(col*255), int(col*255), 255);
        
        int pos[d];
        for(int j = 0; j < d; j++) {
            pos[j] = (int) (p.pos[j]/scale - size/2);
        }
        SDL_Rect rect = {pos[0] - (mr[0] - 500), pos[1] - (mr[1] - 500), size, size};
        SDL_RenderFillRect(gRenderer, &rect);
    }
    
    //Update
    SDL_RenderPresent(gRenderer);
}

int main() {
    if(!init()) {
        cout << "failed init";
    } else {
        cout << "init success" << endl;
        SDL_Event e;
        
        double sim_t = 0;
        clock_t t;
        clock_t t2;
        int rt;
        t = clock();
        t2 = clock();
        int i = 0;
        bool quit = false;
        bool pause = false;
        
        while((i<steps || steps == 0) and !quit){
            
input:
            //input during run
            while(SDL_PollEvent(&e) != 0) {
                if(e.type == SDL_QUIT) quit = true;
                if(e.type == SDL_KEYDOWN) {
                    switch(e.key.keysym.sym) {
                        case SDLK_ESCAPE:
                            quit = true;
                            break;
                        case SDLK_t:
                            dt *= 2;
                            cout << "dt = " << dt << endl;
                            break;
                        case SDLK_g:
                            if(dt != 0) dt *= 0.5;
                            cout << "dt = " << dt << endl;
                            break;
                        case SDLK_l:
                            line = !line;
                            break;
                        case SDLK_p:
                            pause = !pause;
                            break;
                    }
                }
            }
            if(pause) goto input; 
            
            //update particles
            if(BH && d==3 && n > 500) BHupdate();
            else update();
            
            //render particles
            
            rt += clock()-t2;
            t2 = clock();
            if(screen && rt > 14) { 
                rt = 0;
                render(i);
            }
            
            i++;
            sim_t += dt;
            
            //Simulation log
            if(i%100 == 0 && sim_log) {
                cout << endl;
                int rem = 0;
                system_mass = 0;
                for(int j = 0; j < n; j++) {
                    if(!particles[j].e) continue;
                    system_mass += particles[j].mass;
                    rem++;
                    if(particles[j].mass > max_mass) max_mass = particles[j].mass; 
                }
                cout << "Simlulation steps: \t" << i << endl;
                cout << "Simulation time: \t" << int(sim_t*1000*0.000011574) << " days" << endl;
                cout << "Real time: \t\t" << float(clock()-t)/(CLOCKS_PER_SEC*60.0) << " minutes" << endl;
                cout << "Particles remaining: \t" << rem << endl;
                cout << "Total mass: \t\t" << system_mass*1e18 << " kg" << endl;
                cout << "Average mass: \t\t" << system_mass/rem*1e18 << " kg" << endl;
                cout << "Largest mass: \t\t" << max_mass*1e18 << " kg" << endl;
                cout << "Center of mass: \t";
                for(int j = 0; j < d; j++) cout << int(mr[j]) << "\t";
                
                cout << endl << "Velocity of CM: \t";
                double v = 0;
                for(int j = 0; j < d; j++) {
                    v += vr[j]*vr[j];
                    cout << int(vr[j]/dt*1e5) << "\t";
                }
                cout << endl << "Absolute vel of CM: \t" << sqrt(v)/dt*1e2 << " km/s" << endl;
            }
        }
        t = clock() - t;
        float sec = ((float)t)/CLOCKS_PER_SEC;
        cout << endl << i << " steps in " << sec << " seconds." << endl;
        cout << "That's " << int(i/sec) << " steps per second!";
    }
    cout << endl;
    close();
    
    return 0;
}



