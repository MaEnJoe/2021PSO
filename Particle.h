#ifndef __PARTICLE__
#define __PARTICLE__

#include "utility.h"

const double c1 = 0.05;//c1 coefficient of PSO
const double c2 = 0.01;//c2 coefficient of PSO
const double w = 0.05;//w coefficient of PSO

class Particle
{
public:
  explicit Particle();
  explicit Particle(double a1,double a2,double a3,double d2,double fitness);
  explicit Particle(double a1,double a2,double a3,double d2);
  ~Particle();

  //position
  double a1,a2,a3;//in rasius
  double d2;
  //velocity
  double a1_dot,a2_dot,a3_dot,d2_dot;
  double fitness;

  Particle *pbest;
  static Particle *gbest;
  static x_y target;
  static double a1_0,a2_0,a3_0,d2_0;//initial position

  void searching(double r1,double r);
  void update_pbest();

  static void set_boundary(x_y a1_b,x_y a2_b,x_y a3_b,x_y d2_b);
  static void set_target(x_y target);
  static void set_initial_position(double,double,double,double,Particle*);
  static bool set_gbest(Particle* swarm,unsigned int particle_num);
  
  //Particle operator=(const Particle& rval);
   
private:

  //position upper bound
  static double a1_ub,a2_ub,a3_ub;//in radius
  static double d2_ub;//in radius
  //position lower bound
  static double a1_lb,a2_lb,a3_lb;//in radius
  static double d2_lb;


};

#endif 
