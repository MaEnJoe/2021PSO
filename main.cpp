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
#define NUM_PARTICLES 16384
#define ITERATION 10000
#define printi(v) printf(#v" = %d\n",v);
int THREAD_NUM = 6;
Particle* swarm;
x_y target;

x_y get_end_tip(Particle* p);
double distance(x_y p,x_y target);
double uniRand(int lb, int ub);

void init()
{
  //create swarms
  if(swarm != NULL )
  {
    delete[] swarm;
    delete Particle::gbest;
  }
  //set boundary
  Particle::set_boundary({0,PI},{-PI/2,PI/2},{-PI/2,PI/2},{1,15});
  swarm = new Particle[NUM_PARTICLES];
  Particle::set_gbest(swarm,NUM_PARTICLES);
}


int main(int argc, char** argv)
{
 
  if(argc > 1){
    THREAD_NUM = atoi(argv[1]);
  }

  target = {6,4};
  Particle::set_target(target);
  double start,end,refTime,m1Time,m2Time;
  x_y end_tip;
 
  omp_set_num_threads(THREAD_NUM); 
  
  printf("Number of threads: %d\n",THREAD_NUM );
  printf("parallel 1 starts: \n");
  init();
  Particle::set_initial_position(0,0,PI/2,1,&swarm[0]);
  /*ITERATE 1 */
  int counter = 0;  
  bool isConverge = false;
  start = CycleTimer::currentSeconds();
  //#pragma omp parallel 
  {
    for(int iter = 0 ; iter != ITERATION ; iter++)
    {
      //#pragma omp atomic 
      counter++;
      double r1,r2;
      r1 = uniRand(0,1);
      r2 = uniRand(0,1);
      #pragma omp parallel for   
      for(int i = 0 ; i < NUM_PARTICLES; i++)
      {
        swarm[i].searching(r1,r2);
      }
      //if(isConverge)
      //    break;
      //#pragma omp barrier
      if(omp_get_thread_num() == 0)
        isConverge = Particle::set_gbest(swarm,NUM_PARTICLES);
       
    }
  }
  end = CycleTimer::currentSeconds();
  m1Time = end-start;
  printf("\n"); 
  end_tip = get_end_tip(Particle::gbest);
  /*
  printf("gbest is at (%f,%f,%f,%f)\n arm is at (%f,%f)\n",
         Particle::gbest->a1,Particle::gbest->a2,Particle::gbest->a3,Particle::gbest->d2,end_tip.x,end_tip.y);
  printf("fitness(sol) is : %f\n",Particle::gbest->fitness);
  printf("time elapsed: %2.5f\n",m1Time);
  printf("ITER = %d\n",counter);
*/
  /*ITERATE 2 */
  printf("Serial starts: \n");
  init();
  isConverge = false;
  int counter2 = 0;
  start = CycleTimer::currentSeconds();
  for(int iter = 0 ; iter != ITERATION  ; iter++)
  {
    counter2++;
    double r1,r2;
    r1 = uniRand(0,1);
    r2 = uniRand(0,1);
    for(int i = 0 ; i < NUM_PARTICLES; i++)
    {
      swarm[i].searching(r1,r2);
    }
    //if(isConverge)
    //    break;
    isConverge = Particle::set_gbest(swarm,NUM_PARTICLES);
  }
  end = CycleTimer::currentSeconds();
  refTime = end-start;

  printf("\n"); 
  end_tip = get_end_tip(Particle::gbest);
/*
  printf("gbest is at (%f,%f,%f,%f)\n arm is at (%f,%f)\n",
         Particle::gbest->a1,Particle::gbest->a2,Particle::gbest->a3,Particle::gbest->d2,end_tip.x,end_tip.y);
  printf("fitness(sol) is : %f\n",Particle::gbest->fitness);
  printf("time elapsed: %2.5f\n",refTime);
*/
  printf("ITER = %d\n",counter2);
  printf("normalized Speedup = %f\n",refTime/m1Time*counter/counter2);
  return 0;   
}
double uniRand(int lb,int ub)
{
  static thread_local std::mt19937* generator = nullptr;
  if (!generator) {
    std::random_device rd;
    generator = new std::mt19937(rd());
  }
  std::uniform_real_distribution<double> distribution(lb, ub);
  return distribution(*generator);
}
