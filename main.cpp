#include <iostream>
#include <cmath>
#include <cstdio>
#include <random>
#include "Particle.h"
#include "utility.h"
#include <ctime>
#include <omp.h>
#include "common/CycleTimer.h"

#define PI 3.1415926
#define NUM_PARTICLES 10000
#define ITERATION 10000
#define THREAD_NUM 1
Particle* swarm;
x_y target;

std::random_device rd;
std::mt19937 generator(rd());
std::uniform_real_distribution<float> uniRand(0, 1);


x_y get_end_tip(Particle* p);
double distance(x_y p,x_y target);

void init()
{
  //create swarms
  swarm = new Particle[NUM_PARTICLES];
  Particle::set_gbest(swarm,NUM_PARTICLES);
  //set boundary
  Particle::set_boundary({0,PI},{-PI/2,PI/2},{-PI/2,PI/2},{1,3});
}


int main(int argc, char** argv)
{
  init();
  /*set target*/
  target = {2,3};
  Particle::set_target(target);
  Particle::set_initial_position(0,0,PI/2,1,&swarm[0]);
  x_y end_tip;

  double start,end;
  start = CycleTimer::currentSeconds();
  omp_set_num_threads(THREAD_NUM); 
  /*ITERATE*/
  int iter;  
  for(iter = 0 ; iter != ITERATION ; iter++)
  {
    double r1,r2;
    r1 = uniRand(generator);
    r2 = uniRand(generator);
    
    #pragma omp parallel for
    for(int i = 0 ; i < NUM_PARTICLES; i++)
    {
      swarm[i].searching(r1,r2);
    }
    bool isConverge = Particle::set_gbest(swarm,NUM_PARTICLES);
    if(isConverge)
        break;
     
  }
  end = CycleTimer::currentSeconds();
  printf("\n"); 
  end_tip = get_end_tip(Particle::gbest);
  printf("gbest is at (%f,%f,%f,%f)\n arm is at (%f,%f)\n",
         Particle::gbest->a1,Particle::gbest->a2,Particle::gbest->a3,Particle::gbest->d2,end_tip.x,end_tip.y);
  printf("fitness(sol) is : %f\n",Particle::gbest->fitness);
  printf("time elapsed: %2.5f\n",end-start);
  printf("ITER = %d\n",iter);
  return 0;   
}
