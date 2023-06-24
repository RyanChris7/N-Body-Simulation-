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


void update_position(double *x, double *y, double *vx, double *vy, int n) {
    //TODO: update position 
    size_t i = 0;
    while (i < n) {
        x[i] = x[i] + (vx[i] * dt);
        y[i] = y[i] + (vy[i] * dt);

        if (x[i] < 0) {
            vx[i] *= -1;
        } else if (x[i] >= bound_x){
            vx[i] *= -1; 
        } else if (y[i] < 0) {
            vy[i] *= -1;
        } else if (y[i] >= bound_y){
            vy[i] *= -1;
        }
        i++;
    }
}

void update_velocity(double *m, double *x, double *y, double *vx, double *vy, int n) {
    //TODO: calculate force and acceleration, update velocity
    size_t i, j;
    i = 0;
    while(i < n) {
        for (j = 0; j < n; j++) {
            if (i == j) continue;
            double ax, ay, dx, dy, distance, a;
            dx = x[j] - x[i];
            dy = y[i] - y[j];
            distance = (dx * dx) + (dy * dy) + err;
            distance = sqrt(distance);
            if (distance <= radius2){
                distance = radius2;
            } else {
                distance = distance;
            }
            
            a = (gravity_const * m[j]) / (distance * distance);
            if (distance == radius2 || distance <= radius2)
                std::tie(ax,ay) = std::make_tuple(0, 0);

            ax = dx/distance * a;
            ay = dy/distance * a;
            std::tie(ax,ay) = std::make_tuple(ax, ay);

            vx[i] = vx[i] + (dt * ax);
            vy[i] = vy[i] + (dt * ay);
        }
        i++;
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

        update_velocity(m, x, y, vx, vy, n_body);
        update_position(x, y, vx, vy, n_body);

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
    printf("Assignment 2: N Body Simulation Sequential Implementation\n");
    printf("Total running time: %.4f\n", totalTime);
    
    return 0;

}