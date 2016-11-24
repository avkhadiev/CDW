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
#define DC_FIELD        0                           // dc electric field

class CDW {
    public:

    // constructor, destructor
        CDW();                                       // constructor
        virtual ~CDW() = {} ;                        // destructor

    // simulation settings
        int set_ini_time    ( double ini_time    ) ; // initial time
        int set_time_step   ( double time_step   ) ; // length of time step
        int set_num_steps   ( size_t num_steps   ) ; // number of time steps
        
    // initial conditions setting and generation

        // set parameters  
        int set_num_sites   ( size_t num_sites   ) ; // number of lattice sites
        int set_im_spacing  ( size_t im_spacing  ) ; // impurity spacing
        int set_im_strength ( double im_strength ) ; // impurity strength 
        int set_ini_phase   ( double phase       ) ; // initial phases
        int set_temperature ( double temp        ) ; // temperature
        int set_j_model     ( char *j_model      ) ; // elasticity model
        int set_j           ( double j           ) ; // elasticity coefficient
        int set_noise_model ( char *noise_model  ) ; // TODO noise model
        int set_noise       ( double noise       ) ; // TODO noise intensity
        int set_dc_field    ( double dc_field    ) ; // dc electric field

        // display parameters 
        int get_num_sites    ( void ) ;              // number of lattice sites
        int get_im_spacing   ( void ) ;              // impurity spacing
        int get_im_strength  ( void ) ;              // impurity strength 
        int get_ini_phase    ( void ) ;              // initial phases
        int get_temperature  ( void ) ;              // temperature
        int get_j_model      ( void ) ;              // elasticity model
        int get_j            ( void ) ;              // elasticity coefficient
        int get_noise_model  ( void ) ;              // TODO noise model
        int get_noise        ( void ) ;              // TODO noise intensity
        int get_dc_field     ( void ) ;              // dc electric field
        int display_settings ( void ) ;              // display all settings
        
        // generate initial conditions -- 
        // is_setup should be set to 0
        int setup            ( void ) ;              // will set is_setup to 1

        // display lattice
        int show_impurities  ( void ) ; 
        
    // time evolution 
        int run_simulation   ( void ) ;              

    // measurement of observables
    
        // observe a certain number of phases on cites with impurities,
        // and a certain number of phases on sites free of impurities
        int observe_phases    ( size_t impurities, size_t free );
    
    // calculation of statistics
    private:

    // SIMULATION SETTINGS AND DATA STRUCTURES
   
        int is_setup                                // was the simulation set?
                                                    // 1 if was, 0 otherwise
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
         * LATTICE ARRAYS
         *
         * contain pointers to lattice cites
         *      lattice     -- contains pointers to all lattice sites
         *      impurities  -- contains pointers to the sites with impurities
         */
        std::vector<LatticeSite *> lattice( num_sites );
        std::vector<LatticeSite *> impurities;
        
        /*
         * PHASE ARRAY
         *
         * contains pointers to all phases
         */
        std::vector<Phase *> phases( num_sites );

        /*
         * DISTANCE BETWEEN IMPURITIES
         */
        size_t im_spacing ;                         // units: number of sites
        
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
         * ELASTICITY MODELS
         */
        const char* = "infinite range"           ;
        const char* = "nearest neighbor"         ; 
        const char* = "next to nearest neighbor" ;

        /*
         * ELASTICITY
         *
         * describes the elasticity of the CDW.
         *
         * implementation will differ depending on the underlying model:
         *  - infinite range 
         *  - nearest neighbor
         *  - next to nearest neighbor
         *
         *  input:
         *      const Phase *phase, 
         *          a pointer to the phase for which the force is computed
         *
         *  output:
         *      const double, a force of elastic interaction of the phase 
         *          with other phases according to the underlying model
         */
        inline const double j( const Phase *phase ) ;

        /*
         * ELECTRIC FIELDS
         *
         *      DC field strength
         */
        double dc_field = DC_FIELD;
        

    // DYNAMICS PARAMETERS
        /*
         * TIME: START, NUMBER OF STEPS, LENGTH OF EACH STEP
         */
        double time      = INITIAL_TIME ;            // keeps track of time
        double time_step = TIME_STEP    ;            // length of each step
        size_t num_steps = NUM_STEPS    ;            // number of steps

        /*
         *  RATE_OF_CHANGE
         *
         *  given a Phase, outputs the rate of change of the phase according 
         *  to the underlying force field
         *
         *  input:
         *      const *Phase phase, a phase for which the change is 
         *      to be calculated
         *  output:
         *      const double, the calculated change in the phase
         */
        inline const double rate_of_change( const *Phase phase ) ;

        /*
         * UPDATE_PHASE
         *
         * given a Phase *, updates its value and rate of change 
         * for the next time step.
         *
         * requires a specified way to calculate the rate of change.
         *
         * input:
         *      Phase *phase, a pointer to the phase for which values are to 
         *          be changed in place
         * output:
         *      int, 0 if change is successful, 1 otherwise
         */
        int update_phase( Phase * phase ) ;

        /*
         * STEP
         *
         * performs a step in the MD simulation, updating all relevant 
         * variables in place.
         *
         * input:
         *      void
         *
         * output:
         *      int, 0 on sucess, 1 on failure
         */
        int step( void ) ;

        /*
         * EVOLVE
         *
         * performs the evolution of the system for the specified amount of 
         * steps
         *
         * input:
         *      size_t steps
         *
         * output:
         *      int, 0 on success, 1 on failure
         */
        int evolve( size_t steps );
 
    // observables
    /* TODO: need a convenient way to characterize a transition
     *      phases are fixed at the impurity sites; 
     *      
     *      one could track the phase value at the impurity, 
     *      and save a pair of numbers 
     *
     *      ( sin ( [ abs( phi - beta ) mod 2 pi ] - pi/2 ), time ) 
     *
     *      to then graph a time evolution of the phase as a function of time.
     */      

        /*
         * OBSERVED PHASE
         *
         * struct to represent phase behavior at an observed site
         */
        typedef struct {
            // sin ( [ abs( phi - beta ) mod 2 pi ] - pi/2 ) 
            double observed_phase ; 
            // time at which the the observed phase is recorded
            double observed_time  ;
        } ObservedPhase;

        /*
         * GET OBSERVED PHASE
         *
         * given a Phase that records a state of a phase, and a time 
         * corresponding to that state, outputs a struct representing 
         * the behavior of the phase at an observed state at a given time.
         * 
         * this struct may later be used for visualizations
         *
         * input:
         *      const *Phase phase, a phase at a given site
         *      double time, the corresponding time
         *
         * output:
         *      ObservedPhase, the struct representing phase behavior 
         *          at the corresponding time
         *      
         */
        inline ObservedPhase get_observed_phase( const *Phase phase, 
                                                        double time );

        /*
         * ARRAY OF OBSERVES PHASES
         *
         * contain pointers to lattice cites for which the behavior of phases
         * is recorded
         */
        std::vector<ObservedPhases *> observed_phases;

        /*
         * OBSERVE PHASE
         *
         * given a pointer to a Phase, add that phase in the array of 
         * phases to be observed
         *
         * input:
         *      Phase *phase, a pointer to (const) Phase
         *
         * output:
         *      int, 0 on success, 1 on failure
         */
        int observe_phase ( const *Phase phase );

    // statistics
    /*
     * TODO: make sure each simulation provides the following statistics:
     *          - mean transition time per impurity 
     *          - mean velocity of the phase during transition per impurity
     *          - mean number of transitions per unit time per impurity
     */
};
