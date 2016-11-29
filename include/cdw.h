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
#define DEF_TIME_STEP   0.01                        // default time step
#define DEF_NUM_STEPS   50000                       // number of time steps
// model-dependent constants
#define DEF_SITES       12                          // default num of sites
#define DEF_IM_STRENGTH 0                           // default strength
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
    double  phase          ;     // in units of Pi
    double  rate_of_change ;     // in units of Pi, initially 0
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
    int     is_impurity      ;              // 0 if doesn't, 1 if does
    double  im_strength      ;              // between 0 and 1
    double  im_phase         ;              // between 0 and 2 Pi
    Phase   phase            ;              // hosts a phase 
} LatticeSite;


/*
 * OBSERVED PHASE
 *
 * struct to represent phase behavior at an observed site
 */
typedef struct {
    // sin ( [ abs( phi - beta ) mod 2 pi ] - pi/2 ) 
    double phase ;
    double rate_of_change ;
    // time at which the the observed phase is recorded
    double time ;
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
        std::string get_noise_model  ( void ) { return noise_model ; } // TODO              
        double      get_dc_field     ( void ) { return dc_field    ; } 
        
        // vector values
        const std::vector<LatticeSite> get_lattice () { 
            return lattice; 
        }
        const std::vector<const LatticeSite *> get_impurities() {
            return impurities; 
        }
        const std::vector<const LatticeSite *> get_observed_sites() {
            return observed_sites; 
        }
        
        // print out all the settings
        void display_settings ( void );               
        
        // generate initial conditions -- 
        // is_setup should be set to 0
        int setup            ( void ) ;              // will set is_setup to 1
        int uninstall        ( void ) ;              // will set is_setup to 0

        // display lattice
        void display_lattice  ( void ) ; 
        
    // time evolution 
        void run_simulation   ( void ) ;              

    // measurement of observables
    
        // observe all sites with impurities, 
        // FIXME: currently not customizable. 
        //        need a nice way to specify which sites to observe
        int observe_sites( void );
    
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
        std::vector<const LatticeSite *> impurities;
        
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
        double temperature;                  // temperature, K

        /*
         * THERMAL NOISE
         *
         * TODO implementation of thermal noise
         */
        inline double noise( void ) ;               // thermal noise
        std::string noise_model     ;

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
         *          double, a force of elastic interaction of the phase 
         *          with other phases according to the underlying model
         */
        inline double j( const Phase *phase ) ;

        /*
         * ELECTRIC FIELDS
         *
         *      DC field strength
         */
        double dc_field;

        /*
         * INITIALIZE SITE
         *
         * given a pointer to LatticeSite, initializes it as a non-impurity
         * site with a phase with zero velocity and value 
         * specified by ini_phase
         *
         * input: 
         *      LatticeSite *, pointer to a lattice site
         *
         * output: 
         *      int, 0 on success, 1 on failure
         */
        int initialize_site( LatticeSite *site );

        /*
         * GENERATE IMPURITY
         *
         * given a pointer to an initialized LatticeSite, makes it an
         * impurity of specified phase.
         * the strength of the phase is determined by this->im_strength;
         *
         * input: 
         *      LatticeSite *,      pointer to a lattice site
         *      double im_strength, strength of the impurity, <= MAX_IM_STRENGTH
         *      double im_phase,    phase of the impurity between 0 and 2 pi
         *
         * output: 
         *      int, 0 on success, 1 on failure
         */
        int generate_impurity( LatticeSite *site,
                                double im_strength,
                                double im_phase );

        /*
         * GENERATE LATTICE
         *
         *     - creates a vector of LatticeSite specified by num_sites
         *     - adds impurities in the cites according to 
         *          im_spacing, im_strength, and im_phase
         *
         * input:
         *      void
         *
         * output:
         *      void
         */
        void generate_lattice ( void );

        /*
         * ADD IMPURITIES 
         *
         * given an impurity spacing and  generates impurities 
         * in the lattice and fills a vector of pointers to them for reference.
         *
         * input:
         *      size_t im_spacing, impurity spacing in units of lattice sites.
         *
         * output: 
         *      int, 0 on success, 1 on failure
         */
        int add_impurities( size_t im_spacing );
        

    // DYNAMICS PARAMETERS
        /*
         * TIME: START, NUMBER OF STEPS, LENGTH OF EACH STEP
         */
        double ini_time  ;                            // specifies initial time
        double time      ;                            // keeps track of time
        double time_step ;                            // length of each step
        size_t num_steps ;                            // number of steps

        /*
         *  COMPUTE RATE
         *  
         *  input:
         *      size_t i, the index of the site that hosts the phase which 
         *          rate of change is to be computed
         *
         *  output:
         *      double, the value of the rate
         */
        double compute_rate( size_t i ) ;

        /*
         * UPDATE_SITE
         *
         * requires a specified way to calculate the rate of change.
         * once the rate of change has been computed and saved at the site,
         * and the current value of the phase is no longer needed, updates the 
         * value of the phase at the site
         *
         *  input:
         *      size_t i, the index of the site that hosts the phase which 
         *          rate of change is to be computed
         *  output:
         *      void, updates the phase in place
         */
        inline void update_site( size_t i ) ;

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
         *      void
         */
        void step( void ) ;

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
         *      void
         */
        void evolve( size_t steps );
 
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
        inline ObservedPhase get_observed_phase( const LatticeSite *site, 
                                                    double time );

        /*
         * ARRAY OF OBSERVES PHASES
         *
         * contain pointers to lattice cites for which the behavior of phases
         * is recorded
         */
        std::vector<const LatticeSite *> observed_sites;

        /*
         * ADD OBSERVED SITE
         *
         * given a pointer to a LatticeSite, add the pointer to that site in 
         * the array of sites to be observed
         *
         * input:
         *      LatticeSite *site, a pointer to a (const) LatticeSite
         *
         * output:
         *      int, 0 on success, 1 on failure
         */
        inline int add_observed_site ( const LatticeSite *site );

    // statistics
    /*
     * TODO: make sure each simulation provides the following statistics:
     *          - mean transition time per impurity 
     *          - mean velocity of the phase during transition per impurity
     *          - mean number of transitions per unit time per impurity
     */
   
    /*
     * mean momentum
     */
    double mean_momentum;

    /*
     * UPDATE MEAN MOMENTUM
     *
     * proportional to current density
     * update at each evolution step
     *
     */
    void update_mean_momentum ( void ) ;


    /*
     * mean deviation from impurity phase
     */
    double mean_dist_to_im_phase;

    /*
     * curious to see how far the phase at impurity site is from im_phase
     *
     * use avg( fmod( phi, im_phase ) ) for each site
     * then compute average across all sites
     */
    void update_mean_dist_to_im_phase ( void ) ;

    /*
     * for checking momentum conservation
     */
    double zero_momentum;

    /*
     *
     * excluding electric field, momentum should be conserved
     *
     */
    void check_zero_momentum ( void ); 

    // results output

    /*
     * WRITE_OBSERVED_SITES
     *
     * using an array of pointers to the sites being observed, 
     * records a corresponding ObservedPhase data struct for each site into 
     * a separate file.
     *
     * To be called at each evolution state in the MD simulation.
     * The result is an text file per each observed site, that contains the 
     * observed phase, the corresponding rate of change, and corresponding 
     * times.
     *
     * the files can be later plotted with Python
     *
     * input:  
     *      void
     *
     * output:
     *      int, 0 on success, 1 on failure
     */
    int write_observed_sites( void );

    /*
     * WRITE_OBSERVED_SITES
     *
     * input:  
     *      void
     *
     * output:
     *      int, 0 on success, 1 on failure
     */
    int write_statistics( void );
};
