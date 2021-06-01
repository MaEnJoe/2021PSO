#include "Particle.h"
#include "utility.h"
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <random>
#include <limits>
#include <cmath>

#define debugf(a) printf(#a" = %f\n",a);
#define debugi(a) printf(#a" = %d\n",a);

x_y get_end_tip(Particle* p)
{ 
  const double d1 = 1;
  const double d3 = 1;
 double x,y;
 x = d1*cos(p->a1)+p->d2*cos(p->a1+p->a2)+d3*cos(p->a1+p->a2+p->a3);
 y = d1*sin(p->a1)+p->d2*sin(p->a1+p->a2)+d3*sin(p->a1+p->a2+p->a3);  
 return {x,y};
}

double distance(x_y p,x_y target)
{
  double dist = hypot(p.x-target.x,p.y-target.y);
  return dist;
}
double evaluate(Particle* p,x_y target)
{
  double dist = distance(get_end_tip(p),target);
  double energy_cost = (fabs(p->a1-Particle::a1_0)+
                       fabs(p->a2-Particle::a2_0)+
                       fabs(p->a3-Particle::a3_0)+
                       fabs(p->d2-Particle::d2_0))*0.1;
  return dist+energy_cost;
}

static double PI = 3.1415926;
static std::random_device rd;
static std::mt19937 generator(rd());
static std::uniform_real_distribution<float> uniRand(-1, 1);


/*static variable*/
Particle* Particle::gbest = NULL;
double Particle::a1_ub = 0.0;
double Particle::a2_ub = 0.0;
double Particle::a3_ub = 0.0;
double Particle::d2_ub = 0.0;
double Particle::a1_lb = 0.0;
double Particle::a2_lb = 0.0;
double Particle::a3_lb = 0.0;
double Particle::d2_lb = 0.0;
x_y Particle::target;
double Particle::a1_0 = 999.0;
double Particle::a2_0 = 999.0;
double Particle::a3_0 = 999.0;
double Particle::d2_0 = 999.0;

double factor = 0.01;

/*initailize xij and vij randomly*/
Particle::Particle()
{
  this->a1 = uniRand(generator)*PI;
  this->a2 = uniRand(generator)*PI/2;
  this->a3 = uniRand(generator)*PI/2;
  this->d2 = uniRand(generator) + 2;

  this->a1_dot = uniRand(generator)*factor;
  this->a2_dot = uniRand(generator)*factor;
  this->a3_dot = uniRand(generator)*factor;
  this->d2_dot = uniRand(generator)*factor;
  this->fitness = evaluate(this,target);
  this->pbest = new Particle(this->a1,this->a2,this->a3,this->d2,this->fitness);
}

/*for creating pbest and gbest*/
Particle::Particle(double a1,double a2,double a3,double d2,double fitness)
{
  this->a1 = a1;
  this->a2 = a2;
  this->a3 = a3;
  this->d2 = d2;
  this->fitness = fitness;
}
Particle::Particle(double a1,double a2,double a3,double d2)
{
  this->a1 = a1;
  this->a2 = a2;
  this->a3 = a3;
  this->d2 = d2;
  this->a1_dot = uniRand(generator)*factor;
  this->a2_dot = uniRand(generator)*factor;
  this->a3_dot = uniRand(generator)*factor;
  this->d2_dot = uniRand(generator)*factor;
  this->fitness = evaluate(this,target);
  this->pbest = new Particle(this->a1,this->a2,this->a3,this->d2,this->fitness);
 }

bool Particle::set_gbest(Particle* swarm,unsigned int particle_num)
{
  /*input: the swarm of the particles,the numnber of how many particle
   output: return if converge. Convergence is determined by counter > 100,
           if counter > 100: the gbest haven't been updated for 100 iteration
  */
  
   
  double gbest_fitness;
  if(gbest == NULL)
  {
    double a1_ = uniRand(generator)*PI;
    double a2_ = uniRand(generator)*PI/2;
    double a3_ = uniRand(generator)*PI/2;
    double d2_ = uniRand(generator) + 2;
    double realbig = std::numeric_limits<double>::max();
    gbest = new Particle(a1_,a2_,a3_,d2_,realbig);
    return false;
  }
  else
  {
     gbest_fitness = Particle::gbest->fitness;
  }
  static int counter = 0;
  counter++;
  for(int ii = 0 ; ii != particle_num ; ii++)
  {
    if( swarm[ii].fitness < gbest_fitness )
    {
      counter = 0;
      //copy xij
      *Particle::gbest = swarm[ii];
      gbest_fitness = Particle::gbest->fitness;
    }
  }
  if(counter > 100)
      return true;
  else
      return false;

}


/*static */
void Particle::set_boundary(x_y a1_b,
                            x_y a2_b,
                            x_y a3_b,
                            x_y d2_b)
{
  a1_ub = a1_b.y;
  a2_ub = a2_b.y;
  a3_ub = a3_b.y;
  d2_ub = d2_b.y;
  a1_lb = a1_b.x;
  a2_lb = a2_b.x;
  a3_lb = a3_b.x;
  d2_lb = d2_b.x;
}
/*static*/
void Particle::set_initial_position(double a1_00,double a2_00,double a3_00,double d2_00,Particle* p)
{
  a1_0= a1_00;
  a2_0= a2_00;
  a3_0= a3_00;
  d2_0= d2_00;

  /*someone in the swarm must be at the initial position*/
  Particle pp(a1_0,a2_0,a3_0,d2_0);
  assert(p != NULL);
  *p = pp;
  /*make sure initila condition is not out of boundary*/  
  assert(a1_0 <= a1_ub && a1_0 >= a1_lb);
  assert(a2_0 <= a2_ub  && a2_0 >= a2_lb);
  assert(a3_0 <= a3_ub  && a3_0 >= a3_lb);
  assert(d2_0 <= d2_ub && d2_0 >= d2_lb);
}
static int count = 0;
void Particle::searching(double r1,double r2)
{
  //1. update vij
  assert(pbest != NULL);
  assert(gbest != NULL);
  assert(a1_0 != 999.0 && "initail position is not set");
  a1_dot = w*a1_dot + c1*r1*(pbest->a1 - a1) + c2*r2*(gbest->a1 - a1);
  a2_dot = w*a2_dot + c1*r1*(pbest->a2 - a2) + c2*r2*(gbest->a2 - a2);
  a3_dot = w*a3_dot + c1*r1*(pbest->a3 - a3) + c2*r2*(gbest->a3 - a3);
  d2_dot = w*d2_dot + c1*r1*(pbest->d2 - d2) + c2*r2*(gbest->d2 - d2);
  //2. update xij
  a1 = a1+a1_dot;
  a2 = a2+a2_dot;
  a3 = a3+a3_dot;
  d2 = d2+d2_dot;

  //3. check whether out of boundary
  bool mask_a1,mask_a2,mask_a3,mask_d2;
  mask_a1 = a1 > a1_ub || a1 < a1_lb;
  mask_a2 = a2 > a2_ub || a2 < a2_lb;
  mask_a3 = a3 > a3_ub || a3 < a3_lb;
  mask_d2 = d2 > d2_ub || d2 < d2_lb;

  if( mask_a1|| mask_a2 || mask_a3 || mask_d2)
  { 
    this->a1 = uniRand(generator)*PI;
    this->a2 = uniRand(generator)*PI/2;
    this->a3 = uniRand(generator)*PI/2;
    this->d2 = uniRand(generator) + 2;
  }


  //4. update fitness 
  this->fitness = evaluate(this,Particle::target);
  //5. update pbest
  if( this->fitness < this->pbest->fitness)
  { /*copy xij form this to pbest */
    *pbest = *this;
  }
}

void Particle::set_target(x_y target)
{
  Particle::target = target;
}

Particle::~Particle()
{
}
/*Particle Particle::operator=(const Particle& rval)
{
  printf("a1 = %f\n",this->a1);
  this->a1 = rval.a1; 
  this->a2 = rval.a2;
  this->a3 = rval.a3;
  this->d2 = rval.d2;
  this->fitness = rval.fitness;
}*/
