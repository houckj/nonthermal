#ifndef USER_NONTHERMAL_H
#define USER_NONTHERMAL_H

typedef struct Particle_Type Particle_Type;
struct Particle_Type
{
   Particle_Type *next;
   char *method;
   int (*spectrum)(Particle_Type *, double, double *);
   double (*momentum_min)(Particle_Type *);
   double (*momentum_max)(Particle_Type *);
   double *params;
   unsigned int num_params;
   double mass;
};
#define NULL_PARTICLE_TYPE {NULL,NULL,NULL,NULL,NULL,NULL,0,0.0}
#define PARTICLE_METHOD(n,num,m,mn,mx) {NULL,(n),m,mn,mx,NULL,num,0.0}

#define NONTHERMAL_PDF_MODULE(name,p,o) \
   extern int Pdf_##name##_init(Particle_Type *p, char *o); \
          int Pdf_##name##_init(Particle_Type *p, char *o)

#endif
