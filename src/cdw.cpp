#include <cstdlib>
#include <iostream>
#include <cassert>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <random>

#include "../include/cdw.h"

#define UNUSED(x) (void)(x)
const double PI = 3.141592653589793;

// random number generation
// construct a trivial random generator engine from a time-based seed:
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::ranlux24 generator( seed );
std::uniform_real_distribution<double> noise_distribution(0,1);
std::uniform_real_distribution<double> phase_distribution(0,2);

// declare known elasticity models
std::string inf_range          = "infinite range"           ;
std::string near_neighbor      = "nearest neighbor"         ;
std::string next_near_neighbor = "next to nearest neighbor" ;

// declare known noise models
std::string none               = "none"                     ; 
std::string small              = "small"                    ; 
std::string medium             = "medium"                   ;
std::string large              = "large"                    ;

// constructor
CDW::CDW( void ) {
    // simulation-dependent settings
    is_setup    = 0                  ;          // was the simulation set?
    ini_time    = DEF_INI_TIME       ;          // initial time
    time_step   = DEF_TIME_STEP      ;          // length of time step
    num_steps   = DEF_NUM_STEPS      ;          // number of time steps
    // model-dependent settings
    num_sites   = DEF_SITES          ;          // number of lattice sites
    im_spacing  = DEF_IM_SPACING     ;          // impurity spacing
    im_strength = DEF_IM_STRENGTH    ;          // impurity strength
    ini_phase   = DEF_INI_PHASE      ;          // initial phases
    temperature = DEF_TEMPERATURE    ;          // temperature
    j_model     = "nearest neighbor" ;          // elasticity model
    j_strength  = DEF_J_STRENGTH     ;          // elasticity coefficient
    dc_field    = DEF_DC_FIELD       ;          // dc electric field
    // vector variables
    lattice         = std::vector<LatticeSite>( num_sites ) ;
    impurities      = std::vector<const LatticeSite *>( 0 )        ;
    observed_sites  = std::vector<const LatticeSite *>( 0 ) ;
};

// destructor
// implemented simply because of the number of parameters
CDW::~CDW( void ) {};  

// simulation settings
/*
 * all of the methods below return 0 on successful initialization of 
 * parameters, and 1 on failure.
 * 
 * the methods will fail if the simulation has already been set up -- 
 * i.e., if CDW::is_setup = 1;
 */

    /*
     * set initial time
     * 
     * input:
     *      double ini_time, intial time to initialize
     * 
     * output:
     *      int 0, on success, 1 on failure
     */
    int CDW::set_ini_time ( double ini_time ) {
        if (is_setup) {
            printf("%s\n", "The simulation has already been set up");
            return 1;
        }
        CDW::time = ini_time;
        return 0;
    };
    
    /*
     * set time step
     * 
     * input:
     *      double tine_step, time_step to initialize
     * 
     * output:
     *      int 0, on success, 1 on failure
     */
    int CDW::set_time_step ( double time_step ) {
        if (is_setup) {
            printf("%s\n", "The simulation has already been set up");
            return 1;
        }
        CDW::time_step = time_step;
        return 0;
    };

    /*
     * set number of steps
     * 
     * input:
     *      double num_steps, num_steps to initialize
     * 
     * output:
     *      int 0, on success, 1 on failure
     */
    int CDW::set_num_steps ( size_t num_steps ) {
        if (is_setup) {
            printf("%s\n", "The simulation has already been set up");
            return 1;
        }
        CDW::num_steps = num_steps;
        return 0;
    };

// initial conditions setting and generation
/*
 * all of the methods below return 0 on successful initialization of 
 * parameters, and 1 on failure.
 * 
 * the methods will fail if the simulation has already been set up -- 
 * i.e., if CDW::is_setup = 1;
 */

    /*
     * set number of sites
     *
     * specified number of sites should be between MIN_SITES and MAX_SITES
     *
     * input:
     *      size_t num_sites, number of sites to initialize
     *
     * output:
     *      int 0, on success, 1 on failure
     */
    int CDW::set_num_sites ( size_t num_sites ) {
        if (is_setup) {
            printf("%s\n", "The simulation has already been set up");
            return 1;
        }
        if ( (num_sites < MIN_SITES) || (num_sites > MAX_SITES) ) {
            printf("%s %d and %d\n", 
                    "specified number of sites should be between",
                    MIN_SITES,
                    MAX_SITES);
            return 1;
        }
        else {
            CDW::num_sites = num_sites;
            return 0;
        }
     }

    /*
     * set impurity spacing
     * 
     * input:
     *      double im_spacing, impurity spacing to initialize
     * 
     * output:
     *      int 0, on success, 1 on failure
     */
    int CDW::set_im_spacing ( size_t im_spacing ) {
        if (is_setup) {
            printf("%s\n", "The simulation has already been set up");
            return 1;
        }
        CDW::im_spacing = im_spacing;
        return 0;
    };

    /*
     * set impurity strength
     * 
     * input:
     *      double im_strength, impurity strength to initialize
     * 
     * output:
     *      int 0, on success, 1 on failure
     */
    int CDW::set_im_strength ( double im_strength ) {
        if (is_setup) {
            printf("%s\n", "The simulation has already been set up");
            return 1;
        }
        CDW::im_strength = im_strength;
        return 0;
    };


    /*
     * set initial phase
     * 
     * input:
     *      double ini_phase, initial phase to initialize at all sites
     * 
     * output:
     *      int 0, on success, 1 on failure
     */
    int CDW::set_ini_phase ( double ini_phase ) {
        if (is_setup) {
            printf("%s\n", "The simulation has already been set up");
            return 1;
        }
        CDW::ini_phase = ini_phase;
        return 0;
    };

    /*
     * set temperature
     * 
     * input:
     *      double ini_phase, initial phase to initialize at all sites
     * 
     * output:
     *      int 0, on success, 1 on failure
     */
    int CDW::set_temperature ( double temperature ) {
        if (is_setup) {
            printf("%s\n", "The simulation has already been set up");
            return 1;
        }
        CDW::temperature = temperature;
        return 0;
    };

    /*
     * set j_model
     * 
     * input:
     *      char *j_model, a model for elasticity to be used in the simulation:
     *          - "infinite range"           for the infinite range model
     *          - "nearest neighbor"         for the nearest neighbor model
     *          - "next to nearest neighbor" for the next-to-nearest neighbor 
     *                                          model 
     * output:
     *      int 0, on success, 1 on failure
     */
    int CDW::set_j_model ( std::string j_model ) {
        if (is_setup) {
            printf("%s\n", "The simulation has already been set up");
            return 1;
        }

        // match the input string to known models:
        if ( j_model == inf_range ) {
            this->j_model = j_model;
        }
        else if ( j_model == near_neighbor ) {
            this->j_model = j_model;
        }
        else if ( j_model == next_near_neighbor ) {
            this->j_model = j_model;
        }
        else {
            printf("%s '%s' %s: %s, %s, or %s\n",
                    "the input string",
                    j_model.c_str(),
                    "was not recognized as one of the following options",
                    inf_range.c_str(),
                    near_neighbor.c_str(),
                    next_near_neighbor.c_str());
            return 1;
        }
        return 0;
    };


    /*
     * set j_strength
     * 
     * input:
     *      double j_strength, the elasticity coefficient to initialize
     * 
     * output:
     *      int 0, on success, 1 on failure
     */
    int CDW::set_j_strength( double j_strength ) {
        if (is_setup) {
            printf("%s\n", "The simulation has already been set up");
            return 1;
        }
        CDW::j_strength = j_strength;
        return 0;
    };

    
    /*
     * set_noise_model
     * 
     * input:
     *      char *noise_model, a model for noise to be used in the simulation 
     * 
     * output:
     *      int 0, on success, 1 on failure
     */
    int CDW::set_noise_model ( std::string noise_model ) {
        if (is_setup) {
            printf("%s\n", "The simulation has already been set up");
            return 1;
        }

        // match the input string to known models:
        if ( noise_model == none  ) {
            this->noise_model = none;
        }
        if ( noise_model == small ) {
            this->noise_model = noise_model;
        }
        else if ( noise_model == medium ) {
            this->noise_model = noise_model;
        }
        else if ( noise_model == large ) {
            this->noise_model = noise_model;
        }
        else {
            printf("%s '%s' %s: %s, %s, or %s\n",
                    "the input string",
                    noise_model.c_str(),
                    "was not recognized as one of the following options",
                    small.c_str(),
                    medium.c_str(),
                    large.c_str());
            return 1;
        }
        return 0;
    };

    /*
     * set DC electric field
     * 
     * input:
     *      doube dc_field, dc field strength to be initialized
     * 
     * output:
     *      int 0, on success, 1 on failure
     */
    int CDW::set_dc_field ( double dc_field ) {
        if (is_setup) {
            printf("%s\n", "The simulation has already been set up");
            return 1;
        }
        CDW::dc_field = dc_field;
        return 1;
    };

    /*
     * display settings
     *
     * input:
     *      void
     *
     * output:
     *      void
     */
    void CDW::display_settings ( void ) {
        printf("%s: %.3f\n", "initial time       " , get_ini_time()        ) ;    
        printf("%s: %.3f\n", "time step          " , get_time_step()       ) ;    
        printf("%s: %zu\n" , "number of steps    " , get_num_steps()       ) ;    
        printf("%s: %zu\n" , "number of sites    " , get_num_sites()       ) ;    
        printf("%s: %zu\n" , "impurity spacing   " , get_im_spacing()      ) ;    
        printf("%s: %.3f\n", "impurity strength  " , get_im_strength()     ) ;    
        printf("%s: %.3f\n", "initial phase      " , get_ini_phase()       ) ;    
        printf("%s: %.3f\n", "temperature        " , get_temperature()     ) ;    
        printf("%s: %s\n"  , "elasticity model   " , get_j_model().c_str() ) ;    
        printf("%s: %.3f\n", "elasticity strength" , get_j_strength()      ) ;    
        printf("%s: %s\n"  , "noise model        " , get_noise_model().c_str() ) ;  
        printf("%s: %.3f\n", "DC electric field  " , get_dc_field()        ) ;    
        return;
    }


    /*
     * INITIALIZE SITE
     *
     * given a pointer to LatticeSite, initializes it as a non-impurity
     * site with a phase with zero velocity and value 
     * specified by ini_phase
     *
     * input: 
     *      LatticeSite *
     *
     * output: 
     *      int, 0 on success, 1 on failure
     */
    int CDW::initialize_site( LatticeSite *site ) {
       
        if ( site == 0 ) {
            printf("%s\n", "invalid pointer to LatticeSite");
            return 1;
        }

        // initial phase has 0 velocity
        Phase ini_phase = { get_ini_phase(), 0.0 };

        site->is_impurity = 0         ;                 // not an impurity 
        site->im_strength = 0         ;                 // impurity strength
        site->im_phase    = 0         ;                 // impurity phase
        site->phase       = ini_phase ;                 // initial phase at site

        return 0;
    }

    /*
     * GENERATE IMPURITY
     *
     * given a pointer to an initialized LatticeSite, makes it an
     * impurity of specified phase and strength.
     *
     * input: 
     *      LatticeSite *, pointer to a lattice cite
     *      double im_strength, strength of the impurity, <= MAX_IM_STRENGTH
     *      double im_phase,    phase of the impurity, in units of pi
     *
     * output: 
     *      int, 0 on success, 1 on failure
     */
    int CDW::generate_impurity( LatticeSite *site, 
                            double im_strength, 
                            double im_phase  ) {
    
        if ( site == 0 ) {
            printf("%s\n", "invalid pointer to LatticeSite");
            return 1;
        }
    
        if ( (im_phase < 0) || (im_phase > 2) ) {
            printf("phase value %f %s\n", 
                        im_phase, 
                        "is not between 0 and 2 (specified in units of Pi)");
            return 1;
        }
        
        if ( im_strength > MAX_IM_STRENGTH ) {
            printf("strength value %f %s %d",
                        im_strength,
                        "is greater than the maximum",
                        MAX_IM_STRENGTH);
            return 1;
        }

        // generate an impurity:
        site->is_impurity = 1           ;               // is an impurity
        site->im_strength = im_strength ;               // specifed strength
        site->im_phase    = im_phase    ;               // specified phase 

        return 0;
    }
    
    
    /*
     * generate lattice
     *      - creates an vector of LatticeSite * specified by num_sites
     *      - adds impurities in the sites according to 
     *          im_spacing, im_strength, and im_phase.
     */
    void CDW::generate_lattice( void ) {

        // get required lattice dimensions 
        const size_t num_sites = (const size_t)get_num_sites();    

        // resize the lattice vector
        this->lattice.resize( num_sites );

        for(size_t i = 0; i < num_sites; ++i) {
            initialize_site( &(lattice.at(i)) );
        }
    }

    /*
     * ADD IMPURITIES 
     * 
     * TODO any random assignment of impurity strengths and phases can be 
     *      implemented in this function
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
    int CDW::add_impurities( size_t im_spacing ) {
        
        // get required lattice dimensions 
        const size_t num_sites = (const size_t)get_num_sites();

        if ( im_spacing > num_sites ){
            printf("specified impurity spacing %zu > number of sites %zu\n",
                    im_spacing,
                    num_sites);
            return 1;
        }

        if ( im_spacing == 0 ){
            printf("specified impurity spacing %zu is 0",
                    im_spacing);
            return 1;
        }

        //double im_phase     = DEF_IM_PHASE;     // CHANGED TO UNIFORM RANDOM
        double im_strength  = get_im_strength();
    
        this->impurities.clear();

        LatticeSite *required_site;              
        for( size_t i = (im_spacing - 1); i < num_sites; i += im_spacing) {
            required_site = &( this->lattice.at(i) );
            generate_impurity( required_site, 
                                im_strength, 
                                phase_distribution( generator ) );
            this->impurities.push_back( (const LatticeSite *)required_site );
        }
        return 0;
    }

    /*
     * setup
     *
     * generates all the initial conditions:
     *      - generates lattice with generate_lattice();
     *      - adds impurities with 
     *      - sets is_setup to 1
     */
    int CDW::setup ( void ) {
        if (is_setup) {
            printf("%s\n", "The simulation has already been set up");
            return 1;
        }
        
        generate_lattice();                     // generate lattice
        add_impurities( get_im_spacing() );     // add impurities
        is_setup = 1;                           // simulation has been setup

        return 1;
    }

    /*
     * uninstall
     *
     * clears lattice site and impurity vectors
     * sets is_setup to 0;
     *
     * input:
     *      void
     *
     * output:
     *      int, 0, on success, 1 on failure
     */
    int CDW::uninstall ( void ) {
        if ( !(is_setup) ) {
            printf("%s\n", "The simulation has not been yet set up");
            return 1;
        }

        this->lattice.clear();
        this->impurities.clear();
        this->observed_sites.clear();

        is_setup = 0;

        return 0;
    }

    /* display lattice 
     *
     * a simple function to display current lattice / impurities to the user
     *
     * input:
     *      void
     * 
     * output:
     *      void
     */
    void CDW::display_lattice ( void ) {
       const std::vector<LatticeSite> my_lattice = get_lattice();
       
       // get required lattice dimensions 
       const size_t num_sites = (const size_t)get_num_sites();

       printf("\n");
       printf("representation of the lattice:\n");
       for( size_t i = 0; i < num_sites; ++i) {
            my_lattice.at(i).is_impurity? printf(" X %.2f ", 
                                            my_lattice.at(i).im_phase) 
                                        : printf(" _ ");
       }
       printf("\n");

       return;
    }

// Observables
/*
 * functions to observe specific sites on the lattice
 */

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
    inline int CDW::add_observed_site ( const LatticeSite *site ) {
        
        if ( site == 0  ) {
            printf("%s\n", "invalid pointer to LatticeSite");
            return 1;
        }

        this->observed_sites.push_back( site );

        return 0;
    }


    /* observe all sites with impurities, 
     *
     * input: 
     *      void
     *
     * output:
     *      int, 0 on success, 1 on failure
     */
    int CDW::observe_sites( void ){
        if ( !(is_setup) ) {
            printf("%s\n", "The simulation has not been yet set up");
            return 1;
        }
        
        // add all sites in the list of obsered sites
        for(size_t i = 0; i < num_sites; ++i ){
            add_observed_site( &( lattice.at(i)) );
        }

        
        return 0;
    }
    
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
    inline ObservedPhase CDW::get_observed_phase( const LatticeSite *site, 
                                                    double time ){
        ObservedPhase observed;
        double cdw_phase = site->phase.phase ;
        double pin_phase = site->im_phase    ;
        
        // [ abs( phi - beta ) mod 2 pi ] - pi/2 in units of pi
        double normalized_phase = sin( (fmod( fabs( cdw_phase - pin_phase ), 
                                                    2.0 ) 
                                        - 0.5) * PI );

        observed.phase          = normalized_phase;  
        observed.rate_of_change = site->phase.rate_of_change;
        observed.time           = time;

        return observed;
    }

// dynamics
/*
 * functions for updating the state of the system
 */

    /*
     *  NOISE
     *
     *  input:
     *      void:
     *
     *   output: 
     *      void:
     */
     inline double CDW::noise( void ) {
        
        double random = noise_distribution( generator );
        double noise_force;

        if (random >= 0.5) {
            noise_force = 1;
        }
        if (random < 0.5) {
            noise_force = -1;
        }

        if( get_noise_model() == none ) {
            return 0;
        }
        else if( get_noise_model() == small ){
            return 0.2 * noise_force; 
        }
        else if( get_noise_model() == medium ){
            return 0.5 * noise_force; 
        }
        else if( get_noise_model() == large ){
            return 1.0 * noise_force;
        }
        else{
            return 0;
        }
     }

    /*
     *  UPDATE_RATE
     *
     *  input:
     *      size_t i, the index of the site that hosts the phase which 
     *          rate of change is to be computed
     *
     *  output:
     *      double rate, the corresponding rate
     */
     double CDW::compute_rate( size_t i ) {
       if (j_model == near_neighbor){

            LatticeSite *site = &( lattice.at(i) );
            
            double cdw_phase    = site->phase.phase;
            double pin_phase    = site->im_phase;
            double pin_strength = site->im_strength;
            
            double pinning_force = (-1) * pin_strength * sin( (cdw_phase - pin_phase) * PI ) ;
            
            double elastic_force = 0;
            // respect periodic boundary conditions
            size_t right_index = (i + 1) % num_sites;
            size_t left_index  = (i + num_sites - 1) % num_sites;
            double right_phase = lattice.at( right_index ).phase.phase;     
            double left_phase  = lattice.at( left_index  ).phase.phase;
            // whether a phase is to the right or to the left does not matter.
            // i.e., when a neighbor phase is greater, it's always pushing 
            // in the positive direction, regardless of whether it's a right 
            // neighbor or a left neighbor.
            elastic_force += j_strength * ( right_phase - cdw_phase );
            elastic_force += j_strength * ( left_phase  - cdw_phase );

            double rate = elastic_force + pinning_force + noise() + dc_field;
            site->phase.rate_of_change = rate;

            return rate;
       }
       else if (j_model == next_near_neighbor){
           fprintf(stderr, "j_model %s not yet implemented", 
                   next_near_neighbor.c_str() );
           exit( EXIT_FAILURE );
       }
       else {
          fprintf(stderr, "j_model '%s' not recognized by rate_of_change", 
                    j_model.c_str());
          exit( EXIT_FAILURE );
       }
    }


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
    inline void CDW::update_site( size_t i ) {
        
        LatticeSite *site = &( lattice.at(i) );
        site->phase.phase += site->phase.rate_of_change * time_step;

        return;
    }


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
    void CDW::step( void ) {
        std::vector<LatticeSite> new_lattice = get_lattice();

        // compute all rates of change, don't update sites yet (need old values 
        // to compute all the rates, and the new phases)
        for (size_t i = 0; i < num_sites; ++i) {
            new_lattice.at(i).phase.rate_of_change = compute_rate(i); 
        }

        // update all the phases, and rates
        for (size_t i = 0; i < num_sites; ++i) {
            update_site(i);
            lattice.at(i).phase.rate_of_change = new_lattice.at(i).phase.rate_of_change;
        }

        // update time
        time += time_step;

        return;
    }

    
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
    void CDW::evolve( size_t steps ){
        for(size_t i = 0; i < steps; ++i){
            step();
        }
        return;
    }


    /*
     * WRITE_OBSERVED_SITES
     *
     * using an array of pointers to the sites being observed, 
     * records a corresponding ObservedPhase data struct for each site into 
     * a separate file.
     *
     * To be called at each evolution step in the MD simulation.
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
    int CDW::write_observed_sites( void ){

        size_t observations = observed_sites.size();
        ObservedPhase observed;

        // prepare for writeout
        std::ofstream writeout;
        std::string   base_name_phase    = "output/phase";
        std::string   base_name_momentum = "output/momentum";
        std::string   extension          = ".txt";
        std::string   impurity_string;
        std::string   name_phase;
        std::string   name_momentum;

        // loop over observed sites vector
        for( size_t i = 0; i < observations; ++i ){
            
            // get observed phase
            observed = get_observed_phase( observed_sites.at(i), time );

            int is_impurity = observed_sites.at(i)->is_impurity;
            if (is_impurity) {
                impurity_string = "impure";
            }
            else {
                impurity_string = "free";
            }

            // print out (phase, time)
            name_phase = base_name_phase
                            + "_"
                            + impurity_string
                            + "_"
                            + std::to_string( i ) 
                            + extension;
            writeout.open( name_phase, std::ios::out | std::ios::app );

            if( writeout.is_open() ) {
                writeout << std::to_string( observed.phase );
                writeout << "     ";
                writeout << std::to_string( observed.time  );
                writeout << "\n";
                writeout.close();
            }
            else{
                fprintf(stderr, 
                            "Unable to open file %s\n", 
                            name_momentum.c_str());
                perror( "open:" );
                return 1;
            }

            // print out (momentum, time)

            name_momentum = base_name_momentum 
                                + "_"
                                + impurity_string
                                + "_"
                                + std::to_string( i ) 
                                + extension;
            writeout.open( name_momentum, std::ios::out | std::ios::app );

            if( writeout.is_open() ) {
                writeout << std::to_string( observed.rate_of_change );
                writeout << "     ";
                writeout << std::to_string( observed.time  );
                writeout << "\n";
                writeout.close();
            }
            else{
                fprintf(stderr, 
                            "Unable to open file %s\n", 
                            name_momentum.c_str());
                perror( "open:" );
                return 1;
            }
        }

        return 0;
    }


    /*
     * RUN SIMULATION
     */
    void CDW::run_simulation ( void ) {
        if ( !(is_setup) ) {
            printf("%s\n", "The simulation has not yet been set up");
            return;
        }

        // calculate final time:
        double final_time = ini_time + num_steps * time_step;

        write_observed_sites();         // write out the initial values
        while ( time < final_time ) {
            evolve( 1 );                // it is possible to not write out
                                        // at every step by changing this
                                        // argument
            write_observed_sites();
        }
       
        return;
    }
