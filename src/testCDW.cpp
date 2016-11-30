#include "cdw.cpp"

int test_lattice ( void ) { 
    std::cout<<"**************************************************"<<std::endl;
    std::cout<<"WILL TEST:                                        "<<std::endl;
    std::cout<<"    - LATTICE GENERATION:                         "<<std::endl;
    std::cout<<"**************************************************"<<std::endl;
    
    CDW simulation;                                       // declare simulation
    simulation.display_settings();                        // display settings
    simulation.setup();                                   // setup

    simulation.display_lattice();                         // display lattice
    
    printf("done testing lattice generation!\n");
    return 0;
}


int test_default ( void ) { 
    std::cout<<"**************************************************"<<std::endl;
    std::cout<<"WILL TEST:                                        "<<std::endl;
    std::cout<<"    - RUN DEFAULT SIMULATION                      "<<std::endl;
    std::cout<<"**************************************************"<<std::endl;
    
    CDW simulation;                                       // declare simulation
    
    simulation.set_dc_field( 0.0 );                       // no field                
    simulation.set_noise_amp( 0.0 );                      // no noise
    simulation.display_settings();                        // display settings
    
    simulation.setup();                                   // setup
    simulation.observe_sites();                           // observe some sites

    simulation.run_simulation();                          // run simulation
    
    printf("done testing simulation run!\n");
    return 0;
}

int test_field_small ( void ) { 
    std::cout<<"**************************************************"<<std::endl;
    std::cout<<"WILL TEST:                                        "<<std::endl;
    std::cout<<"    - RUN SMALL FIELD SIMULATION                  "<<std::endl;
    std::cout<<"**************************************************"<<std::endl;
    
    CDW simulation;                                       // declare simulation
    
    simulation.set_im_strength( 1.0 );                     // with impurities
    simulation.set_im_spacing( 4 ) ;                       // with impurities
    simulation.set_dc_field( 0.01 );                      // small field                
    simulation.set_noise_amp( 0.0);                 // no noise
    simulation.display_settings();                        // display settings
    
    simulation.setup();                                   // setup
    simulation.observe_sites();                           // observe some sites

    simulation.run_simulation();                          // run simulation
    
    printf("done testing simulation run!\n");
    return 0;
}

int test_field_large ( void ) { 
    std::cout<<"**************************************************"<<std::endl;
    std::cout<<"WILL TEST:                                        "<<std::endl;
    std::cout<<"    - RUN LARGE FIELD SIMULATION                  "<<std::endl;
    std::cout<<"**************************************************"<<std::endl;
    
    CDW simulation;                                       // declare simulation
    
    simulation.set_im_strength( 1.0 );                     // with impurities
    simulation.set_im_spacing( 4 ) ;                       // with impurities
    simulation.set_dc_field( 1.0 );                       // large field                
    simulation.set_noise_amp( 0.0 );                 // no noise
    simulation.display_settings();                        // display settings
    
    simulation.setup();                                   // setup
    simulation.observe_sites();                           // observe some sites

    simulation.run_simulation();                          // run simulation
    
    printf("done testing simulation run!\n");
    return 0;
}

int test_noise_medium ( void ) { 
    std::cout<<"**************************************************"<<std::endl;
    std::cout<<"WILL TEST:                                        "<<std::endl;
    std::cout<<"    - RUN LARGE FIELD SIMULATION                  "<<std::endl;
    std::cout<<"**************************************************"<<std::endl;
    
    CDW simulation;                                       // declare simulation
    
    simulation.set_dc_field( 0.0 );                       // no field                
    simulation.set_im_strength( 1.0 );                     // with impurities
    simulation.set_im_spacing( 4 ) ;                       // with impurities
    simulation.set_noise_amp( 0.5 );                      // medium noise
    simulation.display_settings();                        // display settings
    
    simulation.setup();                                   // setup
    simulation.observe_sites();                           // observe some sites

    simulation.run_simulation();                          // run simulation
    
    printf("done testing simulation run!\n");
    return 0;
}

int test_noise_large ( void ) { 
    std::cout<<"**************************************************"<<std::endl;
    std::cout<<"WILL TEST:                                        "<<std::endl;
    std::cout<<"    - RUN LARGE FIELD SIMULATION                  "<<std::endl;
    std::cout<<"**************************************************"<<std::endl;
    
    CDW simulation;                                       // declare simulation
    
    simulation.set_im_strength( 1.0 );                // with impurities
    simulation.set_im_spacing( 4 );                      // with impurities
    simulation.set_dc_field( 0.0 );                       // no field                
    simulation.set_noise_amp( 1.0 );                      // large noise
    simulation.display_settings();                        // display settings
    
    simulation.setup();                                   // setup
    simulation.observe_sites();                           // observe some sites

    simulation.run_simulation();                          // run simulation
    
    printf("done testing simulation run!\n");
    return 0;
}

int main(int argc, const char** argv) {
    std::cout<<"**************************************************"<<std::endl;
    std::cout<<"                    CDW TESTER                    "<<std::endl;
    std::cout<<"**************************************************"<<std::endl;

    std::vector<std::string> args(0);

    if (argc == 1) {
        // HAVE TO SPECIFY AT LEAST ONE TEST SUITE
        fprintf(stderr, "no test suites to run!\n");
        return 1;
    }
    else {
        for( int i = 0; i < (argc - 1); ++i ){
            args.push_back( argv[i + 1] );
            std::string current_arg = args.at(i);

            // PARSE THROUGH ARGUMENTS AND MATCH THEM TO TEST SUITES
            if ( current_arg == "lattice" ) {
                test_lattice();
            }
            else if ( current_arg == "default"     ) {
                test_default();
            }
            else if ( current_arg == "field_small"     ) {
                test_field_small();
            }
            else if ( current_arg == "field_large"     ) {
                test_field_large();
            }
            else if ( current_arg == "noise_medium"    ) {
                test_noise_medium();
            }
            else if ( current_arg == "noise_large"     ) {
                test_noise_large();
            }
            else {
                fprintf(stderr, 
                        "argument '%s' is not recognized\n", 
                        current_arg.c_str());
            }
        }
        fprintf(stdout, "done testing!\n");
    }

    return 0;
}
