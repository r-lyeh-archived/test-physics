// notes:
// Acceleration can be calculated, without time, by recording not 2 but 3 positions of each particle.
// Acceleration = (Pnow-Pold) - (Pold-Polder)

// @todo:
// tests with step: big/small
// tests with acl: on/off
// tests with vel: on/off

#include <stdio.h>
#include <math.h>

#define lprintf if( j == 0 || j == loops -1 ) if( i == 0 || i == count - 1 ) printf

template<typename flt = double>
struct tests {

  flt finalx, finalv;

#define RESET() \
    x = xo = xoo = v = 0.f; a = 10.f;

  void measure() {
    // use simple euler with a very small step for this
    unsigned i, j, loops = 10, hz = 10000000, count = hz;
    flt x = 0.f;
    flt xo = 0.f;
    flt xoo = 0.f;
    flt v = 10.f;
    flt a = 0.f;
    flt dt = 1.f/count;

    RESET();
    for( j=0; j < loops; ++j )
    for (i=0; i < count; i++) {
      x += v*dt;
      v += a*dt;
    } // for

    finalx = x;
    finalv = v;
  }

  void compare( flt x, flt v ) {
   flt errx = (x-finalx) * (x-finalx);
   flt pctx = finalx == 0 ? 0 : abs( 100 - 100. * abs(x/finalx));
   flt errv = (v-finalv) * (v-finalv);
   flt pctv = finalv == 0 ? 0 : abs( 100 - 100. * abs(v/finalv));
   printf( "Error X: %f%% (%f vs %f)\n", pctx, x, finalx );
   printf( "Error V: %f%% (%f vs %f)\n\n", pctv, v, finalv );
  }

  tests() {

    measure();

    // measure deviation on many loops with big timestep
    unsigned i, j, loops = 10, hz = 1, count = hz;
    flt x = 0.f;
    flt xo = 0.f;
    flt xoo = 0.f;
    flt v = 10.f;
    flt a = 0.f;
    flt dt = 1.f/count;
    flt dto = 0.f;

    // The Euler method is a first order numerical procedure for solving ordinary differential equations
    // (ODEs) with a given initial value. It is the most basic explicit method for numerical integration
    // for ordinary differential equations.
    // if accel is constant...
    RESET();
    printf("Simple explicit Euler: (requires: x,v)(optional: a)\n");
    for( j=0; j < loops; ++j ) for ( i=0; i < count; i++ ) {
      x += v*dt;
      v += a*dt;
      lprintf("\tV: %f X: %f\n",v,x);
    } // for
    compare( x, v );

    // The Euler–Cromer algorithm or symplectic Euler method or Newton-Stormer-Verlet (NSV) method is a
    // modification of the Euler method for solving Hamilton's equations, a system of ordinary differential
    // equations that arises in classical mechanics. It is a symplectic integrator, which is a class of geometric
    // integrators that is especially good at simulations of systems of undamped oscillators. Due to this property,
    // it preserves energy better than the standard Euler method and so is often used in simulations of orbital
    // mechanics.  It is a first order integrator.
    RESET();
    printf("Semi implicit Euler (NSV: Newton-Stormer-Verlet) (Symplectic Euler) (Euler Cromer): (requires: x,v)(optional: a)\n");
    for( j=0; j < loops; ++j ) for ( i=0; i < count; i++ ) {
      v += a*dt;
      x += v*dt;
      lprintf("\tV: %f X: %f\n",v,x);
    } // for
    compare( x, v );

    RESET();
    // The "suvat" method takes its name from the variables in the equations of motion. This is a first order
    // method that is very similar to the Euler method. The only difference is the addition of an acceleration
    // term to the equation for the change in position.
    // See http://www.niksula.hut.fi/~hkankaan/Homepages/gravity.html
    // Often written also as:
    //  x += v*dt + .5f*a*dt*dt;
    //  v += a*dt;
    // if accel is variable...
    printf("Half semi implicit Euler: (requires: x,v)(optional: a)\n");
    for( j=0; j < loops; ++j ) for ( i=0; i < count; i++ ) {
      v += a*dt*0.5;
      x += v*dt;
      v += a*dt*0.5;
      lprintf("\tV: %f X: %f\n",v,x);
    } // for
    compare( x, v );

   // Verlet integration is a numerical integration method originally designed for calculating the trajectories of
   // particles in molecular dynamics simulations.  The velocity verlet variant directly calculates velocity.
   // It is a second order integrator. Unfortunately, this method is unsuitable for simulations where the
   // acceleration is dependent on velocities, such as with a damper. The modified verlet variant below uses
   // some tricks to try to get around this limitation, but this method will only provide first-order accuracy.
    RESET();
    printf("Velocity-less Verlet: (requires x,a)\n");
    {
    for( j=0; j < loops; ++j ) for ( i=0; i < count; i++ ) {
      v = (x - xo) / dt;
      xo = x;
      v += a*dt;
      x += v*dt;
      lprintf("\tV: %f X: %f\n",v,x);
    }} // for
    compare( x, v );

    RESET();
    printf("Velocity-less Verlet 'LEAPFROG': (requires x,a)\n");
    {
    for( j=0; j < loops; ++j ) for ( i=0; i < count; i++ ) {
      flt dtn = dt - (dt - dto); dto = dt;
      v = (x - xo) / dt;
      xo = x;
      v += a*dt;
      x += v*(dt * 0.5f + dtn * 0.5f);
      lprintf("\tV: %f X: %f\n",v,x);
    }} // for
    compare( x, v );

    // note: dt/dto and dt/dt are constants here. therefore, they can be precomputed before loop
    RESET();
    printf("Time-Corrected Velocity-less Verlet: (requires x,a)\n");
    {
    flt xo = x, dto = dt;
    for( j=0; j < loops; ++j ) for ( i=0; i < count; i++ ) {
      v = ( (x - xo) * (dt/dto) + a*dt*dt ) / dt;
      dto = dt;
      xo = x;
      x += v*dt;
      lprintf("\tV: %f X: %f\n",v,x);
    }} // for
    compare( x, v );

    RESET();
    printf("Velocity-less Euler (custom): (requires x,a)\n");
    for( j=0; j < loops; ++j ) for ( i=0; i < count; i++ ) {
      v = (x-xo) / dt;
      xo = x;
      v += a*dt;
      x += v*dt;
      lprintf("\tV: %f X: %f\n",v,x);
    } // for
    compare( x, v );

    RESET();
    printf("Velocity-less Euler (custom): (requires x,a)\n");
    for( j=0; j < loops; ++j ) for ( i=0; i < count; i++ ) {
      v = (x-xo) / dt;
      xo = x;
      v += a*dt;
      x += v*dt;
      lprintf("\tV: %f X: %f\n",v,x);
    } // for
    compare( x, v );
  }

};

int main() {
  tests<> test;
}
