#include "cdw.cpp"

int free( void ) {
    std::cout<<"**************************************************"<<std::endl;
    std::cout<<"CONTROL FOR:                                      "<<std::endl;
    std::cout<<"    - RUN LARGE FIELD SIMULATION                  "<<std::endl;
    std::cout<<"    - WITHOUT IMPURITIES                          "<<std::endl;
    std::cout<<"**************************************************"<<std::endl;
    
    CDW simulation;                                       // declare simulation
    
    simulation.set_num_sites( 100 );
    simulation.set_im_spacing( 4 );
    simulation.set_noise_amp(   0.5 );                    // no noise
    simulation.set_im_strength( 0.0 );                    // no impurities
    
    simulation.display_settings();                        // display settings
    
    for( double i = 0.0; i <= 1.0; i = i + 0.02 ){
        simulation.set_dc_field( i );                     // set new field
        simulation.setup();                               // setup
        simulation.display_settings();                    // display settings
        simulation.run_simulation();                      // run simulation
        simulation.uninstall();                           // uninstall
    }
    
    
    printf("done running control experiment!\n");
    return 0;
}

int impure( void ) {
    std::cout<<"**************************************************"<<std::endl;
    std::cout<<"CONTROL FOR:                                      "<<std::endl;
    std::cout<<"    - RUN LARGE FIELD SIMULATION                  "<<std::endl;
    std::cout<<"    - WITH IMPURITIES                             "<<std::endl;
    std::cout<<"**************************************************"<<std::endl;
    
    CDW simulation;                                       // declare simulation
    
    simulation.set_num_sites( 100 );
    simulation.set_im_spacing( 4 );
    simulation.set_noise_amp(  0.8 );                     // some noise?
    simulation.set_im_strength( 1.0 );                    // with impurities
    
    simulation.display_settings();                        // display settings
   
    for( double i = 0.0; i <= 0.2; i = i + 0.002 ){
        simulation.set_dc_field( i );                     // set new field
        simulation.setup();                               // setup
        simulation.display_settings();                    // display settings
        simulation.run_simulation();                      // run simulation
        simulation.uninstall();                           // uninstall
    }

    for( double i = 0.2; i <= 1.0; i = i + 0.02 ){
        simulation.set_dc_field( i );                     // set new field
        simulation.setup();                               // setup
        simulation.display_settings();                    // display settings
        simulation.run_simulation();                      // run simulation
        simulation.uninstall();                           // uninstall
    }
    
    printf("done running control experiment!\n");
    return 0;
}


int main(int argc, const char** argv) {
    std::cout<<"**************************************************"<<std::endl;
    std::cout<<"                    CONTROL FIELD                 "<<std::endl;
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
            if ( current_arg == "free" ) {
                free();
            }
            else if ( current_arg == "impure"     ) {
                impure();
            }
            else {
                fprintf(stderr, 
                        "argument '%s' is not recognized\n", 
                        current_arg.c_str());
            }
        }
        fprintf(stdout, "done running experiments!\n");
    }

    return 0;
}
