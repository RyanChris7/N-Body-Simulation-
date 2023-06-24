#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <tuple>
#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
#include "./headers/logger.h"


int n_body;
int n_iteration;

int n_omp_threads;
double totalTime;


void generate_data(double *m, double *x,double *y,double *vx,double *vy, int n) {
    // TODO: Generate proper initial position and mass for better visualization
    srand((unsigned)time(NULL));
    for (int i = 0; i < n; i++) {
        m[i] = rand() % max_mass + 1.0f;
        x[i] = 2000.0f + rand() % (bound_x / 4);
        y[i] = 2000.0f + rand() % (bound_y / 4);
        vx[i] = 0.0f;
        vy[i] = 0.0f;
    }
}


void update_position(double *x, double *y, double *vx, double *vy, int i) {
    size_t a = 0;
    while (a < i) {
        x[a] = x[a] + (vx[a] * dt);
        y[a] = y[a] + (vy[a] * dt);

        if (x[a] < 0) {
            vx[a] *= -1;
        } else if (x[a] >= bound_x){
            vx[a] *= -1; 
        } else if (y[a] < 0) {
            vy[a] *= -1;
        } else if (y[a] >= bound_y){
            vy[a] *= -1;
        }
        a++;
    }
}

void update_velocity(double *m, double *x, double *y, double *vx, double *vy, int i) {
    //TODO: calculate force and acceleration, update velocity
    size_t a, b;
    a = 0;
    while(a < i) {
        for (b = 0; b < i; b++) {
            if (a == b) continue;
            double ax, ay, dx, dy, distance, acel;
            dx = x[b] - x[a];
            dy = y[a] - y[b];
            distance = (dx * dx) + (dy * dy) + err;
            distance = sqrt(distance);
            if (distance <= radius2){
                distance = radius2;
            } else {
                distance = distance;
            }
            
            acel = (gravity_const * m[b]) / (distance * distance);
            if (distance == radius2 || distance <= radius2)
                std::tie(ax,ay) = std::make_tuple(0, 0);

            ax = dx/distance * acel;
            ay = dy/distance * acel;
            std::tie(ax,ay) = std::make_tuple(ax, ay);

            vx[a] = vx[a] + (dt * ax);
            vy[a] = vy[a] + (dt * ay);
        }
        a++;
    }
}

void master() {
    double* m = new double[n_body];
    double* x = new double[n_body];
    double* y = new double[n_body];
    double* vx = new double[n_body];
    double* vy = new double[n_body];                      

    generate_data(m, x, y, vx, vy, n_body);

    Logger l = Logger("sequential", n_body, bound_x, bound_y);

    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        
        //TODO: choose better threads configuration
        omp_set_num_threads(n_omp_threads);
        #pragma omp parallel for
        for (int i = 0; i < n_body; i++) {
            update_velocity(m, x, y, vx, vy, i);
        }

        omp_set_num_threads(n_omp_threads);
        #pragma omp parallel for
        for (int i = 0; i < n_body; i++) {
            update_position(x, y, vx, vy, i);
        }

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = t2 - t1;

        printf("Iteration %d, elapsed time: %.3f\n", i, time_span);
        totalTime += time_span.count();

        l.save_frame(x, y);

        #ifdef GUI
        glClear(GL_COLOR_BUFFER_BIT);
        glColor3f(1.0f, 0.0f, 0.0f);
        glPointSize(2.0f);
        glBegin(GL_POINTS);
        double xi;
        double yi;
        for (int i = 0; i < n_body; i++){
            xi = x[i];
            yi = y[i];
            glVertex2f(xi, yi);
        }
        glEnd();
        glFlush();
        glutSwapBuffers();
        #else

        #endif
    }

    delete[] m;
    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;
    
}


int main(int argc, char *argv[]){
    
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);
    n_omp_threads = atoi(argv[3]);

    #ifdef GUI
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(500, 500);
    glutCreateWindow("N Body Simulation Sequential Implementation");
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    gluOrtho2D(0, bound_x, 0, bound_y);
    #endif
    totalTime = 0.0;
    master();

    printf("Student ID: 119010507\n"); // replace it with your student id
    printf("Name: Ryan Christopher\n"); // replace it with your name
    printf("Assignment 2: N Body Simulation OpenMP Implementation\n");
    printf("Total running time: %.4f\n", totalTime);
    
    return 0;
}


