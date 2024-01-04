#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include <mutex>
#include <bits/stdc++.h>
#include "utimer.hpp"

// Global verobsity variables.
bool print_logs = false;
bool print_times = false;
bool subtask_time_analysis = false;

int num_workers = 1;

class City {     // The class
  public:       // Access specifier
    int id;    // Index in the .tsp file
    float x;  // Latitude
    float y; // Longitude
    
    City(int index = 0, float longitude = 0, float latitude = 0) { // Constructor with parameters
      id = index;
      x = latitude;
      y = longitude;
    }

    float distance(City dest) { 
        return pow(pow(x - dest.x, 2) + pow(y - dest.y, 2), 0.5);
    }
};

class Route { 
  public: 
    int N;
    std::vector<City> *cities;
    std::vector<int> ordering;
    float distance = 0;

    Route() = default;

    Route(int n, std::vector<City> *cities_list){
        N = n;
        cities = cities_list;

        for (int i = 0; i< N; i++){
            ordering.push_back(i);
        } 

        std::random_shuffle (ordering.begin(), ordering.end());
    }

    Route(int n, std::vector<City> *cities_list, std::vector<int>::iterator ordering_cities){
        N = n;
        cities = cities_list;
        ordering.assign(ordering_cities, ordering_cities + N);
    }


    float compute_total_distance(){ 
        distance = 0;
        if(print_times && subtask_time_analysis){
            utimer full_route_distance("Full route distance: ");
            for(int i=0; i < N; i++){
                float delta_d = cities->at(ordering[i]).distance(cities->at(ordering[(i+1)%N]));
                distance = distance + delta_d;
            }
        } else { 
            for(int i=0; i < N; i++){
                float delta_d = cities->at(ordering[i]).distance(cities->at(ordering[(i+1)%N]));
                distance = distance + delta_d;
            }
        }
        return distance;
        //@TODO save values in an upper triangular (half of a symmetric) matrix??? 
    }

    Route PMX(Route mating_route, int crossing_point){
        if(print_times && subtask_time_analysis){
            utimer PMX("Breeding (PMX): ");
        }
        if (crossing_point == -1){
            crossing_point = rand()%N;
        }
        
        std::vector<int> child_ordering;
        child_ordering.reserve(N);

        std::vector<int>::iterator partner_begin = mating_route.ordering.begin();
        std::vector<int>::iterator child_begin = child_ordering.begin();
        std::vector<int>::iterator chid_end = child_begin + N;

        for(int i = 0; i < N; i++){
            child_ordering[i] = ordering[i];
        }
    
        for (int i = 0; i < crossing_point; i++){
            int mate_gene = *(partner_begin + i);
         
            // If mutating gene is already present in 
            auto correspondant_idx = std::find(child_begin + i, chid_end, mate_gene);
            if (correspondant_idx != chid_end){
                child_ordering[(correspondant_idx - child_begin)] = *(child_begin + i);
            } 

            child_ordering[i] = mate_gene;
        }

        Route r = Route(N, cities, child_ordering.begin());
        return r;
    }      

    void mutate(){
        int first = rand () % N;
        int second = rand () % (N-1);

        if (second == first){
            ++second;
        }  

        int copy = ordering[first];

        ordering[first] = ordering[second];
        ordering[second] = copy;
    } 

    void print_cities_ids () { 
        for(unsigned int i = 0; i < N; i++){
            std::cout << ordering[i] << '\n';
        }
    }
};

class Nation {       // The class
  private:          // Access specifier
    std::vector<City> cities;

  public:
    int N;

  public:
    Nation(std::string filename) { // Constructor with parameters
        std::ifstream myfile;
        myfile.open("./NATIONS/" + filename);
        std::string myline;

        int point_pos = filename.find(".tsp");
        N = std::stoi(filename.replace(point_pos, 4, "").replace(0, 2, ""));

        if ( myfile.is_open() ) {
            int skipline = 7;
            while ( myfile ) { // equivalent to myfile.good()
                std::getline (myfile, myline);
                skipline --;
                if ( skipline < 0 && skipline > - (N + 1)) {

                    std::string s;
                    std::vector<std::string> raw_line;
                    std::stringstream ss(myline);

                    int j = 0;
                    while (std::getline(ss, s, ' ')) {
                        // store token string in the vector
                        raw_line.push_back(s);
                        j ++;
                    
                        if (j % 3 == 0){
                            addCity(City(std::stoi(raw_line[0]), std::stof(raw_line[1]), std::stof(raw_line[2])));
                            raw_line.clear();
                        }
                    }
                }
            }
        } else {
            std::cout << "Couldn't open file\n";
        }
    }

    void addCity(City city) { 
        cities.push_back(city);
    }

    City get_city(int id){ 
        return cities[id];
    }

    void print_cities_ids () { 
        for(unsigned int i = 0; i < N; i++){
            std::cout << cities[i].id << '\n';
        }
    }

    Route random_route() { 
        if(print_times && subtask_time_analysis){
            utimer random_route_generation("Random route generation: ");
            return Route(N, &cities);
        }

        return Route(N, &cities);
    }
};

class Population { 
    private:
        bool isOrdered = false;
        std::vector<float> probabilities;
        std::discrete_distribution<> picker;
        std::random_device rd;
        std::mt19937 gen = std::mt19937(rd());
        int N;

    public: 
        int K; // Number of genes, a.k.a. routes 
        Route best_route;
        std::vector<Route> routes;
    
        Population(int k, Nation *nation){ 
            K = k;
            N = nation->N;

            std::mutex mutex;
            std::vector<std::thread> mapThreads;

            std::vector<Route>* routed_pointer = &routes;
            std::vector<float>* probabilities_pointer = &probabilities;
            
            for (int i = 0; i< num_workers; i++){
                // Use threads to generate the K initial routes contemporarely. 
                // Assume numThreads == K.
                mapThreads.emplace_back([&nation, routed_pointer, k, &mutex]() {
                    std::vector<Route> random_routes;
                    for (int j = 0; j < k/num_workers; j++){
                        random_routes.push_back((*nation).random_route());
                    } 
                    std::lock_guard<std::mutex> lock(mutex);
                    for (int j = 0; j < k/num_workers; j++){
                        (*routed_pointer).push_back(random_routes[j]);
                    }
                });
            }
            probabilities = std::vector<float>(K, 1/K);
                // Wait for all Map threads to finish
            for (auto& thread : mapThreads) {
                thread.join();
            }

            rank_all();
            best_route = routes[0]; 
        }

        // The distance computation step should be parallelized within the evolution phase and not in the ranking/reduce pahse. 
        void rank_all() { 
            float total_inverse_distance = 0;
            for (int i = 0; i< K; i++){
                total_inverse_distance += 1/routes[i].compute_total_distance(); // This also updates the distance value for the route object. 
            }
            
            std::sort(routes.begin(), routes.begin() + K, [](Route a, Route b){ return a.distance > b.distance; });  //bizarre cpp lambda function syntax

            // This part can be "parallelized" being moved at the beginning of the evolution step and made individual for each worker, 
            // considering K different distributions to avoid self picking. 
            for (int i = 0; i< K; i++){
                probabilities[i] = (1/routes[i].distance)/total_inverse_distance;
            }

            picker = std::discrete_distribution<>(probabilities.begin(), probabilities.end());
            return;
        }

        Route pick_route(bool allow_self_pick = true){ 
            //@TODO implement prevention from a route picking itself as partner
            if (print_times && subtask_time_analysis){ 
                utimer pick_route("Partner pickup: ");
                return routes[picker(gen)]; 
            } else { 
                return routes[picker(gen)]; 
            }
           
        }

        void evolve_once(float mutation_rate = 0){ 
            std::vector<std::thread> mapThreads;
            std::mutex mutex;
            for (int i = 0; i < num_workers; i++){
                mapThreads.emplace_back([=, &mutex](){
                    for (int j = 0; j < K/num_workers; j++) { 
                        int crossing_point = rand()% N;//rand()%2 ==1 ? rand()% N : -1;
                        // Randomly choose which PMX ordering to choose
                        int route_idx = i*K/num_workers + j;

                        mutex.lock();
                        Route current_route = routes[route_idx];
                        Route mate_route = pick_route();
                        mutex.unlock();

                        Route new_route;
                        
                        if(rand() % 2 == 1){
                            new_route = current_route.PMX(mate_route, crossing_point);
                        } else { 
                            new_route = mate_route.PMX(current_route, crossing_point);
                        }

                        // Randomly check if the mutation should happen or not. 
                        if ( mutation_rate != 0 && (rand() % int(1/mutation_rate)) == 1){ 
                            if(print_times && subtask_time_analysis){
                                utimer mutationTime("Mutation time:");
                                new_route.mutate();
                            } else { 
                                new_route.mutate();
                            }
                        }

                        routes[route_idx] = new_route;
                    }
                });
            }
            for (auto& thread : mapThreads) {
                thread.join();
            }
        }

        void evolve(int epochs, float mutation_rate){ 
            for(int i = 0; i < epochs; i++ ){ 
                if (print_times){
                    utimer ranking("Evolution (Breed + Mutate): ");
                    evolve_once(mutation_rate);
                } else { 
                    evolve_once(mutation_rate);
                }
                if (print_times){
                    utimer ranking("Ranking (Sort + Normalize): ");
                    rank_all();
                } else { 
                    rank_all();
                }

                // Hold a copy of the best route ever found, in case evoltuion cause a distance increasement.
                if (routes[0].compute_total_distance() < best_route.compute_total_distance()){ 
                    best_route.distance = routes[0].distance;
                    best_route.ordering = routes[0].ordering;
                } 

                if(print_logs){ 
                    std::cout << "Epoch." << i << ")  Best Distance:";
                    std::cout << best_route.compute_total_distance() << '\n';
                }
            }
        }
};


int main (int argc, char* argv[]){
    // This code allows to read the .tsp files skipping the fixed-lenght headers (7 lines).
    // The remaining lines just corrresponds to the coordinates of the cities in the nation.
    std::string filename = argc > 1 ? argv[1] : "wi29.tsp";
    std::cout << filename << '\n';

    // User can pass all the algorithm parameters. 
    int epochs =  argc > 2 ? atoi(argv[2]) : 500;
    std::cout << "Iterations: " << epochs << "\n";
    int genes_num = argc > 3 ? atoi(argv[3]) : 1000;
    std::cout << "Genes number: " << genes_num << "\n";
    float mutation_rate = argc > 4 ? atof(argv[4]) : 0.0001;
    std::cout << "Mutaton Rate: " << mutation_rate <<"\n";

    num_workers = argc > 5 ? atoi(argv[5]) : 10;
    std::cout << "Number of worker: " << num_workers <<"\n";
    std::cout << "Maximum number of workers possible: "<< std::thread::hardware_concurrency() << "\n"; 
    // Costumizable verbosity: 
    // 2 ( > 1) for logs only
    // 1 For logs and timers
    // 0 for no prints at all
    // -1 for timers only
    int verbosity = argc > 6 ? atoi(argv[6]) : 0;
    print_logs = verbosity > 0; 
    print_times = verbosity == -1 || verbosity == 1;

    // Last param to decide wether to save logs on a file or just use the terminal. 
    bool redirect_logs = argc > 7; 
    if (redirect_logs && atoi(argv[7]) == -1){ 
        subtask_time_analysis = true;
    }

    if (redirect_logs){
        std::time_t ct = std::time(0);
        //TODO include a log folder and a t= 
        std::stringstream sstm;
        sstm << "LOGS/PAR/E"<< epochs << "G"<< genes_num << "_" << ct << ".txt";
        std::string logfile_path = sstm.str(); 
        std::cout << "Trying redirecting the logs to :" << logfile_path << "\n";
        std::freopen(logfile_path.c_str(),"w",stdout);
    }

    Nation nation = Nation(filename); 

    Population* genome;

    if(print_times){
        utimer initialization("Random routes initialization");
        genome = new Population(genes_num, &nation);
    } else {
        genome = new Population(genes_num, &nation);
    }    

    if (print_times){
        utimer overall("Overall execution time: ");
        (*genome).evolve(epochs, mutation_rate);
    } else { 
        (*genome).evolve(epochs, mutation_rate);
    }

    if(print_logs){ 
        (*genome).best_route.print_cities_ids();
    }
}

void lazyRouteTest(Nation nation){ 
    nation.print_cities_ids();

    std::cout << nation.N << '\n';
    std::cout << nation.get_city(1).distance(nation.get_city(9))  << '\n';
    std::cout << nation.random_route().compute_total_distance() << '\n';

    Route route1 =  nation.random_route();
    Route route2 =  nation.random_route();
    
    route1.print_cities_ids();
    std::cout << "++++++++++++++++++++++++++++++++++++++++" << '\n';
    route2.print_cities_ids();
    std::cout << "========================================" << '\n';

    Route child = (route2.PMX(route1, 8));
    child.print_cities_ids();
    std::cout << child.ordering[0] <<'\n';
}

void time_tests(){ 
    if (print_times ){ 
        utimer evolution("Single evolution loop (including ranking):");
        // TODO evakyate the single tasks here.
    }
}