#include "Declearation.h"

float pos[864][3];
float vel[864][3];
float init_vel[864][3];
float acc[864][3];
float last_acc[864][3];
float dist[864];
float boundry_len = 34.7787;

int i =0;
int flag[864];
float sigma = 3.4;

// Have to take care about units
// it's in epsilon/Kb
float epsilon = 120;
float kbm = .0247713;

float last_accx;
float last_accy;
float last_accz;

float delta_time = .0001;
