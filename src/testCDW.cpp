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


int run_simulation ( void ) { 
    std::cout<<"**************************************************"<<std::endl;
    std::cout<<"WILL TEST:                                        "<<std::endl;
    std::cout<<"    - RUN DEFAULT SIMULATION                      "<<std::endl;
    std::cout<<"**************************************************"<<std::endl;
    
    CDW simulation;                                       // declare simulation
    
    simulation.set_dc_field( 0.0 );                      
    simulation.set_noise_model( "none" );                 // noise
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
            else if ( current_arg == "run"     ) {
                run_simulation();
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
