#include <cstdlib>
#include <iostream>
#include <cassert>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include "../include/cdw.h"

#define UNUSED(x) (void)(x)
const double PI = 3.141592653589793;
const double TwoPI = 2 * PI; 

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
    lattice         = std::vector<LatticeSite  >( num_sites ) ;
    impurities      = std::vector<LatticeSite *>(     0     ) ;
    observed_phases = std::vector<ObservedPhase>(     0     ) ;
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

        // declare known elasticity models
        std::string inf_range          = "infinite range"           ;
        std::string near_neighbor      = "nearest neighbor"         ;
        std::string next_near_neighbor = "next to nearest neighbor" ;

        // match the input string to known models:
        if ( j_model == inf_range ) {
        }
        else if ( j_model == near_neighbor ) {
        }
        else if ( j_model == next_near_neighbor ) {
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
        // FIXME STUB
        UNUSED( noise_model );
        return 1;
    };

    /*
     * set_noise
     * 
     * input:
     *      doube noise, noise intensity to be initialized
     * 
     * output:
     *      int 0, on success, 1 on failure
     */
    int CDW::set_noise ( double noise ) {
        // FIXME STUB
        UNUSED( noise );
        return 1;
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
        printf("%s: %s\n"  , "noise model        " , "STUB"                ) ;   //TODO
        printf("%s: %s\n"  , "noise intensity    " , "STUB"                ) ;   //TODO
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

        site->is_impurity = 0         ;             // not an impurity 
        site->im_strength = 0         ;             // impurity strength
        site->im_phase    = 0         ;             // impurity phase
        site->phase       = ini_phase ;             // initial phase at site

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
     *      double im_phase,    phase of the impurity between 0 and 2 pi
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
    
        if ( (im_phase < 0) || (im_phase > TwoPI) ) {
            printf("phase value %f %s\n", 
                        im_phase, 
                        "is not between 0 and 2 pi");
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
        // FIXME
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

        double im_phase     = DEF_IM_PHASE;
        double im_strength  = get_im_strength();
    
        this->impurities.clear();

        LatticeSite *required_site;              
        for( size_t i = (im_spacing - 1); i < num_sites; i += im_spacing) {
            required_site = &( this->lattice.at(i) );
            generate_impurity( required_site, 
                                im_strength, 
                                im_phase );
            this->impurities.push_back( required_site );
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
            my_lattice.at(i).is_impurity? printf("X ") : printf("_ ");
       }
       printf("\n");

       return;
    }
    
