#include <iostream>
#include <cmath>
#include <cstdio>
#include <random>
#include "Particle.h"
#include "utility.h"


#define PI 3.1415926
#define NUM_PARTICLES 30
   
Particle* swarm;

std::random_device rd;
std::mt19937 generator(rd());
std::uniform_real_distribution<float> uniRand(0, 1);

static Particle ppp;

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
  //set target
  Particle::set_target({0,2});
    
  /*ITERATE*/
  for(int iter = 0 ; iter != 100000 ; iter++)
  {
    double r1,r2;
    r1 = uniRand(generator);
    r2 = uniRand(generator);

    for(int i = 0 ; i != NUM_PARTICLES; i++)
    {
      swarm[i].searching(r1,r2);
    }
    Particle::set_gbest(swarm,NUM_PARTICLES);
  }

  x_y end_tip = get_end_tip(Particle::gbest);
  printf("gbest is at (%f,%f,%f,%f)\n arm is at (%f,%f)\n",
         Particle::gbest->a1,Particle::gbest->a2,Particle::gbest->a3,Particle::gbest->d2,end_tip.x,end_tip.y);
  printf("fitness(sol) is : %f\n",Particle::gbest->fitness);

  return 0;   
}
