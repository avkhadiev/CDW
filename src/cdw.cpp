#include <cstdlib>
#include <iostream>
#include <cassert>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include "include/cdw.h"

#define UNUSED(x) (void)(x)

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
        printf("%s: %.3f\n", "initial time"       , get_ini_time()        ) ;    
        printf("%s: %.3f\n", "time step"          , get_time_step()       ) ;    
        printf("%s: %zu\n" , "number of steps"    , get_num_steps()       ) ;    
        printf("%s: %zu\n" , "number of sites"    , get_num_sites()       ) ;    
        printf("%s: %zu\n" , "impurity spacing"   , get_im_spacing()      ) ;    
        printf("%s: %.3f\n", "impurity strength"  , get_im_strength()     ) ;    
        printf("%s: %.3f\n", "initial phase"      , get_ini_phase()       ) ;    
        printf("%s: %.3f\n", "temperature"        , get_temperature()     ) ;    
        printf("%s: %s\n"  , "elasticity model"   , get_j_model().c_str() ) ;    
        printf("%s: %.3f\n", "elasticity strength", get_j_strength()      ) ;    
        printf("%s: %s\n"  , "noise model"        , "STUB"                ) ;   //TODO
        printf("%s: %s\n"  , "noise intensity"    , "STUB"                ) ;   //TODO
        printf("%s: %.3f\n", "DC electric field"  , get_dc_field()        ) ;    
        return;
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

        // clear and resize the lattice vector
        lattice.clear();
        lattice.resize( num_sites );

        // initial phase has 0 velocity
        Phase initial_phase = { get_ini_phase(), 0.0 };

        // initialize first and last element of 
        lattice.at(1).phase = initial_phase;
        
        for(size_t i = 0; i < num_sites; ++i) {
            lattice.at(i).phase = initial_phase;
        }
    }

    /*
     * setup
     *
     * generates all the initial conditions:
     *      - generates lattice with generate_lattice();
     *      - sets is_setup to 1
     */
    int CDW::setup (  ) {
        if (is_setup) {
            printf("%s\n", "The simulation has already been set up");
            return 1;
        }
        // get required lattice dimensions 
        const size_t num_sites = (const size_t)get_num_sites();    
        
        // clear the vectors
        lattice.clear();
        impurities.clear();

        lattice.resize( num_sites );

        // initial phase has 0 velocity
        Phase initial_phase = { get_ini_phase(), 0.0 };
        
        // iterate over lattice sites and phases
        for(size_t i = 0; i < num_sites; ++i) {
        }
        return 1;
    }
