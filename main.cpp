#include <iostream>
#include <cmath>
#include <cstdio>
#include <random>
#include "Particle.h"
#include "utility.h"
#include <ctime>

#define PI 3.1415926
#define NUM_PARTICLES 10000
#define ITERATION 10000
Particle* swarm;
x_y target;

std::random_device rd;
std::mt19937 generator(rd());
std::uniform_real_distribution<float> uniRand(0, 1);


x_y get_end_tip(Particle* p);
double distance(x_y p,x_y target);
bool isConverge(x_y,x_y);

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
  //set target
  target = {2,3};
  Particle::set_target(target);
  Particle::set_initial_position(0,0,PI/2,1,&swarm[0]);
  x_y end_tip;

  clock_t t1 = clock();  
  /*ITERATE*/
  int iter;
  for(iter = 0 ; iter != ITERATION ; iter++)
  {
    double r1,r2;
    r1 = uniRand(generator);
    r2 = uniRand(generator);

    for(int i = 0 ; i != NUM_PARTICLES; i++)
    {
      swarm[i].searching(r1,r2);
    }
    Particle::set_gbest(swarm,NUM_PARTICLES);
    end_tip = get_end_tip(Particle::gbest);
      
  }
  clock_t t2 = clock();
  printf("\n"); 
  end_tip = get_end_tip(Particle::gbest);
  printf("gbest is at (%f,%f,%f,%f)\n arm is at (%f,%f)\n",
         Particle::gbest->a1,Particle::gbest->a2,Particle::gbest->a3,Particle::gbest->d2,end_tip.x,end_tip.y);
  printf("fitness(sol) is : %f\n",Particle::gbest->fitness);
  printf("time elapsed: %1.5f\n",(t2-(double)t1)/CLOCKS_PER_SEC);
  printf("ITER = %d\n",iter);
  return 0;   
}
bool isConverge(x_y p,x_y target)
{
//  printf("manhatan distance =  %f\n",((p.x-target.x)+(p.y-target.y)));
  return (fabs(p.x-target.x)+fabs(p.y-target.y)) < 0.01;
}
