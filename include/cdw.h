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
#define MAX_IM_STRENGTH 3                           // max impurity strength
#define DEF_INI_TIME    0                           // simulation start time
#define DEF_TIME_STEP   0.001                       // default time step
#define DEF_NUM_STEPS   1000                        // number of time steps
// model-dependent constants
#define DEF_SITES       12                          // default num of sites
#define DEF_IM_STRENGTH 1                           // default strength
#define DEF_IM_PHASE    0                           // default pinning phase
#define DEF_IM_SPACING  3                           // default impurity spacing
#define DEF_INI_PHASE   0                           // default initial phase
#define DEF_J_STRENGTH  1                           // default elasticity 
#define DEF_TEMPERATURE 100                         // default temperature, K
#define DEF_ELASTICITY  1                           // default CDW elasticity
#define DEF_DC_FIELD    0                           // dc electric field

/*
 * CDW PHASE
 *
 * struct to represent CDW phase at each lattice site 
 *
 *      - specifies the value of the phase
 *      - specifies current rate of change of the phase
 */
typedef struct {
    double  phase          ;     // not restricted to 2 Pi
    double  rate_of_change ;     // initially 0
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
struct LatticeSite;

struct LatticeSite{
    int     is_impurity      ;              // 0 if doesn't, 1 if does
    double  im_strength      ;              // between 0 and 1
    double  im_phase         ;              // between 0 and 2 Pi
    Phase   phase            ;              // hosts a phase 
    LatticeSite * const next ;              // pointer to next site
    LatticeSite * const prev ;              // pointer to prev site
};

typedef struct LatticeSite LatticeSite;

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
 * Class CDW
 */
class CDW {
    public:

    // constructor, destructor
        CDW( void );                                 // constructor
        virtual ~CDW( void );                        // destructor

    // simulation settings
        
        // set parameters
        int set_ini_time    ( double ini_time    ) ; // initial time
        int set_time_step   ( double time_step   ) ; // length of time step
        int set_num_steps   ( size_t num_steps   ) ; // number of time steps

        // get parameters
        double get_ini_time  ( void ) { return ini_time  ;}
        double get_time_step ( void ) { return time_step ;}
        size_t get_num_steps ( void ) { return num_steps ;}
        
    // initial conditions setting and generation

        // set parameters  
        int set_num_sites   ( size_t num_sites   ) ; // number of lattice sites
        int set_im_spacing  ( size_t im_spacing  ) ; // impurity spacing
        int set_im_strength ( double im_strength ) ; // impurity strength 
        int set_ini_phase   ( double ini_phase   ) ; // initial phases
        int set_temperature ( double temperature ) ; // temperature
        int set_j_model     ( std::string j_model) ; // elasticity model
        int set_j_strength  ( double j_strength  ) ; // elasticity coefficient
        int set_noise_model ( std::string noise_model ) ; // TODO noise model
        int set_noise       ( double noise       ) ; // TODO noise intensity
        int set_dc_field    ( double dc_field    ) ; // dc electric field

        // get parameters 
        size_t      get_num_sites    ( void ) { return num_sites   ; }          
        size_t      get_im_spacing   ( void ) { return im_spacing  ; }               
        double      get_im_strength  ( void ) { return im_strength ; }               
        double      get_ini_phase    ( void ) { return ini_phase   ; }               
        double      get_temperature  ( void ) { return temperature ; }               
        std::string get_j_model      ( void ) { return j_model     ; }               
        double      get_j_strength   ( void ) { return j_strength  ; }               
        std::string get_noise_model  ( void ) { return "STUB"      ; } // TODO              
        std::string get_noise        ( void ) { return "STUB"      ; } // TODO      
        double      get_dc_field     ( void ) { return dc_field    ; } 
        
        // vector values
        const std::vector<LatticeSite> get_lattice () { return lattice ; }
        
        // print out all the settings
        void display_settings ( void );               
        
        // generate initial conditions -- 
        // is_setup should be set to 0
        int setup            ( void ) ;              // will set is_setup to 1
        int uninstall        ( void ) ;              // will set is_setup to 0

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
   
        int is_setup;                               // was the simulation set?
                                                    // 1 if was, 0 otherwise
       
        /*
         * INITIAL PHASE
         *
         * initial phase assigned to phases at all sites;
         */
        double ini_phase;
                                                    
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
        std::vector<LatticeSite> lattice;
        std::vector<LatticeSite *> impurities;
        
        /*
         * DISTANCE BETWEEN IMPURITIES
         */
        size_t im_spacing ;                         // units: number of sites
        
        /*
         * IMPURITY STRENGTH
         */
        size_t im_strength ;                        // impurity strength

    // MODEL PARAMETERS
        
        /* 
         * TEMPERATURE 
         */
        static double temperature;                  // temperature, K

        /*
         * THERMAL NOISE
         *
         * TODO implementation of thermal noise
         */
        inline const double noise( void ) ;         // thermal noise

        /*
         * ELASTICITY MODELS
         */
        double       j_strength ; 
        std::string  j_model    ;
        
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
        double dc_field;

        /*
         * GENERATE LATTICE
         *
         * input:
         *      void
         *
         * output:
         *      void
         */
        void generate_lattice ( void );
        

    // DYNAMICS PARAMETERS
        /*
         * TIME: START, NUMBER OF STEPS, LENGTH OF EACH STEP
         */
        double ini_time  ;                            // specifies initial time
        double time      ;                            // keeps track of time
        double time_step ;                            // length of each step
        size_t num_steps ;                            // number of steps

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
        inline const double rate_of_change( const Phase *phase ) ;

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
        inline ObservedPhase get_observed_phase( const Phase *phase, 
                                                    double time );

        /*
         * ARRAY OF OBSERVES PHASES
         *
         * contain pointers to lattice cites for which the behavior of phases
         * is recorded
         */
        std::vector<ObservedPhase> observed_phases;

        /*
         * ADD OBSERVED PHASE
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
        int add_observed_phase ( const Phase *phase );

    // statistics
    /*
     * TODO: make sure each simulation provides the following statistics:
     *          - mean transition time per impurity 
     *          - mean velocity of the phase during transition per impurity
     *          - mean number of transitions per unit time per impurity
     */
};
