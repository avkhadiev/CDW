/*
 * cdw.h
 * 
 * Class structure for MD simulation of CDWs.
 * Specifies:
 *      - simulation settings and data structures 
 *      - model parameters
 *      - initial conditions generation
 *      - dynamics parameters and time evolution 
 *      - observables and their measurement 
 *      - statistics and their computation
 *      - visualization 
 *      - results output
 */

/* MACROS */
// simulation-dependent constants
#define MAX_SITES       1000                        // max num of lattice sites
#define MIN_SITES       2                           // min num of lattice sites
#define TEMPERATURE     100                         // default temperature, K
#define INITIAL_TIME    0                           // simulation start time
#define TIME_STEP       0.001                       // default time step
#define NUM_STEPS       1000                        // number of time steps
// model-dependent constants
#define MAX_IM_STRENGTH 3                           // max impurity strength
#define DEF_IM_PHASE    0                           // default pinning phase
#define DEF_ELASTICITY  1                           // default CDW elasticity
#define AC_FIELD        0                           // ac electric field 
#define DC_FIELD        0                           // dc electric field

class CDW {
    public:
    // constructor, destructor
        CDW();                                      // constructor
        virtual ~CDW() = {} ;                       // destructor
    // initial conditions generation
    // time evolution 
    // measurement of observables
    // calculation of statistics
    private:

    // SIMULATION SETTINGS AND DATA STRUCTURES
   
        int is_setup                                // was the simulation set?
                                                    // 1 if was, 0 otherwise

        /*
         * LATTICE SITE
         *
         * struct to represent discretized lattice spacing
         *
         *      - can contain an impurity, in which case 
         *      - specifies impurity strength
         *      - specifies impurity phase (the CDW phase to pin).
         */
        typedef struct {
            int     is_impurity ;                   // 0 if doesn't, 1 if does
            double  im_strength ;                   // between 0 and 1
            double  im_phase    ;                   // between 0 and 2 Pi
            Phase * const phase ;                   // hosts a constant phase       
        } LatticeSite;
        
        /*
         * NUMBER OF LATTICE SITES
         */ 
        size_t num_sites ;                          // >= MIN_SITES
                                                    // <= MAX_SITES
        
        /*
         * DISTANCE BETWEEN IMPURITIES
         */
        size_t im_spacing ;                         // units: number of sites
        
        /*
         * CDW PHASE
         *
         * struct to represent CDW phase at each lattice site 
         *
         *      - specifies the value of the phase
         *      - specifies current rate of change of the phase
         *      - specifies the lattice on which the phase is located 
         */
        typedef struct {
            double              phase          ;     // not restricted to 2 Pi
            double              rate_of_change ;     // initially 0
            LatticeSite * const lattice_site   ;     // const position
        } Phase;

    // MODEL PARAMETERS
        
        /* 
         * TEMPERATURE 
         */
        double t = TEMPERATURE ;                    // temperature, K

        /*
         * THERMAL NOISE
         *
         * TODO implementation of thermal noise
         */
        inline const double noise( void ) ;         // thermal noise

        /*
         * ELASTICITY
         *
         * describes the elasticity of the CDW.
         * may be dependent on the distance 
         * (in units of number of lattice sites) between the coupld phases.
         *
         * implementation will differ depending on the underlying model:
         *  - infinite range 
         *  - nearest neighbor
         *  - next to nearest neighbor, etc.
         *
         */
        inline const double j( size_t numSites ) ;
        

    // DYNAMICS PARAMETERS
       
        /*
         * TIME: START, NUMBER OF STEPS, LENGTH OF EACH STEP
         */
        double time      = INITIAL_TIME ;            // keeps track of time
        double time_step = TIME_STEP    ;            // length of each step
        double num_steps = NUM_STEPS    ;            // number of steps

        /*
         *  RATE_OF_CHANGE
         *
         *  given a Phase, outputs the rate of change of the phase according 
         *  to the underlying force field
         *
         *  input:
         *      Phase phase, a phase for which the change is to be calculated
         *  output:
         *      const double, the calculated change in the phase
         */
        inline const double rate_of_change( Phase phase ) ;

        /*
         * CHANGE_PHASE
         *
         * given a Phase *, changes its value and rate of change 
         * for the next time step.
         *
         * requires a specified way to calculate the rate of change.
         *
         *
         * input:
         *      Phase *phase, a pointer to the phase for which values are to 
         *          be changed in place
         * output:
         *      int, 0 if change is successful, 1 otherwise
         */
        int change_phase( Phase * phase ) ;
 
    // observables

    // statistics 
};
