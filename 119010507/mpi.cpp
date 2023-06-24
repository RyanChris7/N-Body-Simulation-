#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <mpi.h>
#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
#include "./headers/logger.h"


int n_body;
int n_iteration;

int rank;
int world_size;

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
            distance = sqrt(r);
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


void slave(){
    // TODO: MPI routine
    double* local_m;
    double* local_x;
    double* local_y;
    double* local_vx;
    double* local_vy;
    // TODO End
}



void master() {
    double* total_m = new double[n_body];
    double* total_x = new double[n_body];
    double* total_y = new double[n_body];
    double* total_vx = new double[n_body];
    double* total_vy = new double[n_body];

    generate_data(total_m, total_x, total_y, total_vx, total_vy, n_body);

    Logger l = Logger("sequential", n_body, bound_x, bound_y);

    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        // TODO: MPI routine
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        
        // TODO End

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = t2 - t1;

        printf("Iteration %d, elapsed time: %.3f\n", i, time_span);

        l.save_frame(total_x, total_y);

        #ifdef GUI
        glClear(GL_COLOR_BUFFER_BIT);
        glColor3f(1.0f, 0.0f, 0.0f);
        glPointSize(2.0f);
        glBegin(GL_POINTS);
        double xi;
        double yi;
        for (int i = 0; i < n_body; i++){
            xi = total_x[i];
            yi = total_y[i];
            glVertex2f(xi, yi);
        }
        glEnd();
        glFlush();
        glutSwapBuffers();
        #else

        #endif
    }

    delete[] total_m;
    delete[] total_x;
    delete[] total_y;
    delete[] total_vx;
    delete[] total_vy;

}




int main(int argc, char *argv[]) {
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	if (my_rank == 0) {
		#ifdef GUI
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize(500, 500); 
		glutInitWindowPosition(0, 0);
		glutCreateWindow("N Body Simulation MPI Implementation");
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glMatrixMode(GL_PROJECTION);
		gluOrtho2D(0, bound_x, 0, bound_y);
		#endif
        master();
	} else {
        slave();
    }

	if (my_rank == 0){
		printf("Student ID: 119010507\n"); // replace it with your student id
		printf("Name: Ryan Christopher\n"); // replace it with your name
		printf("Assignment 2: N Body Simulation MPI Implementation\n");
	}

	MPI_Finalize();

	return 0;
}

