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
  //set boundary
  Particle::set_boundary({0,PI},{-PI/2,PI/2},{-PI/2,PI/2},{1,3});
  Particle::set_target({2,1});
  Particle::set_gbest(swarm,NUM_PARTICLES);

}

void cheating(Particle* p)
{
  p->a1 = PI/2;
  p->a2 = -PI/2;
  p->a3 = 0;
  p->d2 = 1;
  *(p->pbest) = *p;
  x_y cc = get_end_tip(p);
  double fitness = distance(cc,{2,1});
  p->fitness = fitness;
}

int main(int argc, char** argv)
{
  printf("(c1,c2,w) = (%f,%f,%f)\n",c1,c2,w);
  init();

  //----his is cheating----//
  /*cheating(&swarm[0]);
  printf("-----------------------------------------\n");
  Particle::set_gbest(swarm,NUM_PARTICLES);
  x_y cc = get_end_tip(&swarm[0]);
  printf("cheat is at (%f,%f,%f,%f)\n", swarm[0].a1,swarm[0].a2,swarm[0].a3,swarm[0].d2);  

  x_y end_tip1 = get_end_tip(Particle::gbest);
  printf("gbest is at (%f,%f,%f,%f)\n\t\tarm is at (%f,%f)\n",
         Particle::gbest->a1,Particle::gbest->a2,Particle::gbest->a3,Particle::gbest->d2,end_tip1.x,end_tip1.y);
  cc = get_end_tip(Particle::gbest);
  double fitness = distance(cc,{2,1});
  printf("~~~~%f\n",fitness );*/
 //----chearting ends----//


  for(int iter = 0 ; iter != 100000 ; iter++)
  {
    double r1,r2;
    r1 = uniRand(generator);
    r2 = uniRand(generator);

    for(int i = 0 ; i != NUM_PARTICLES; i++)
    {
      swarm[i].searching(r1,r2);
    }
    //printf("cheat = %f,%f,%f,%f\n", swarm[0].a1,swarm[0].a2,swarm[0].a3,swarm[0].d2);  

    Particle::set_gbest(swarm,NUM_PARTICLES);
  }
  //printf("cheat = %f,%f,%f,%f\n", swarm[0].a1,swarm[0].a2,swarm[0].a3,swarm[0].d2);

  x_y end_tip = get_end_tip(Particle::gbest);
  //---print outcome! ---//
  printf("gbest is at (%f,%f,%f,%f)\n arm is at (%f,%f)\n",
         Particle::gbest->a1,Particle::gbest->a2,Particle::gbest->a3,Particle::gbest->d2,end_tip.x,end_tip.y);
  printf("fitness(sol) is : %f\n",Particle::gbest->fitness);

  return 0;   
}
