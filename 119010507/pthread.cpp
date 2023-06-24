#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <pthread.h>
#include <tuple>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
#include "./headers/logger.h"

int n_thd; // number of threads

int n_body;
int n_iteration;
int thread_tracking;
double totalTime;

double *m, *x, *y, *vx, *vy;

pthread_mutex_t mutex;




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

void update_position(double *x, double *y, double *vx, double *vy, int a, int b) {
    //TODO: update position 
    size_t i = a;
    while (i < b) {
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

void update_velocity(double *m, double *x, double *y, double *vx, double *vy, int a, int b, int n) {
    //TODO: calculate force and acceleration, update velocity
    size_t i, j;
    i = a;
    while (i < b){
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


typedef struct {
    //TODO: specify your arguments for threads
    int a;
    int b;
    int diff;
} Args;


void* worker(void* args) {
    Args* my_arg = (Args*) args; 
    int a = my_arg->a;
    int b = my_arg->b;
    int diff = my_arg->diff;

    update_velocity(m, x, y, vx, vy, a, b, n_body);
    update_position(x, y, vx, vy, a, b);
    pthread_exit(NULL);
    // TODO END
}


void master(){
    m = new double[n_body];
    x = new double[n_body];
    y = new double[n_body];
    vx = new double[n_body];
    vy = new double[n_body];

    generate_data(m, x, y, vx, vy, n_body);

    Logger l = Logger("Pthread", n_body, bound_x, bound_y);

    pthread_t threads[n_thd];
    Args args[n_thd];
    size_t i = 0;

    int remainder = n_body % n_thd; // remaider of data
    int store = 0; 
 
    while (i < n_thd) {
        args[i].diff = i;
        args[i].a = store;
        int my_element;
        if (i < remainder){
            my_element = (n_body / n_thd) + 1;
        } else {
            my_element = n_body / n_thd;
        }
        store += my_element;
        args[i].b = store - 1;
        i++;
    }
    
    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        //TODO: assign jobs
        
        for (size_t i = 0; i < n_thd; i++) pthread_create(&threads[i], NULL, worker, &args[i]);
        for (int i = 0; i < n_thd; i++) pthread_join(threads[i], NULL);
        //TODO End

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


int main(int argc, char *argv[]) {
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);
    n_thd = atoi(argv[3]);

    #ifdef GUI
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Pthread");
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0, bound_x, 0, bound_y);
    #endif
    totalTime = 0.0;
    master();

    printf("Total running time: %.4f\n", totalTime);

	return 0;
}

