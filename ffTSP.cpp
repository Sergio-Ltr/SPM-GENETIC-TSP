#include <iostream>
#include <fstream>
#include <string>
#include <bits/stdc++.h>
#include "utimer.hpp"
#include <functional>

#include <ff/ff.hpp>
#include <ff/parallel_for.hpp>

using namespace ff;

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

        std::random_shuffle(ordering.begin(), ordering.end());
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

//The emitter of the farm, it computes the merge phase
class merge_em : public ff_monode{
    private: 
        int epochs ; //number of iteratons
        int workers;
        int counter = 0; //counter used to understand if all worker completed their iteration
        std::function<void()> callback;
        bool first_gen = true;
    
    public:
        merge_em(int epochs, int workers, std::function<void()> callback){
            this->epochs = epochs;
            this->workers = workers;
            this->callback = callback;
        }

        void* svc(void* inp){
            if(epochs != 0){ //still iteration to do
                if(first_gen){ //no need to do anything, just let the workers start
                    broadcast_task(GO_ON);
                    first_gen = false;
                    return GO_ON;
                }
                if(counter < (workers-1)){ //understand if all workers have completed their termination
                    counter++;
                    return GO_ON;
                }
                callback();
                counter = 0; //reset counter for next iteration
                epochs--; //decrease the number of iteration
                broadcast_task(GO_ON);
                return GO_ON;
            }
            else{
                return EOS;
            }
        }
};

//The worker of the farm, it computes Selection, Crossover, Mutation, Fitness
class worker_SCMF : public ff_node{
    private: 
        int worker_idx;
        float mutation_rate;
        std::function<void(int, float)> work;

    public: 
        worker_SCMF(int worker_idx, float mutation_rate ,std::function<void(int, float)> work){ 
            this->worker_idx = worker_idx;
            this->mutation_rate = mutation_rate;
            this->work = work;
        }

        void* svc(void* inp){
            work(worker_idx, mutation_rate);
            ff_send_out(GO_ON);
            return GO_ON;
        }
};

class Population { 
    private:
        std::vector<float> probabilities;
        std::discrete_distribution<> picker;
        std::random_device rd;
        std::mt19937 gen = std::mt19937(rd());
        float total_inverse_distance = 0;
        int N;

    public: 
        int K; // Number of genes, a.k.a. routes 
        Route best_route;
        std::vector<Route> current_gen_routes;
        std::vector<Route> new_gen_routes;
        
        Population(int k, Nation *nation){ 
            K = k;
            N = nation->N;
            
            K = k;
            N = nation->N;

            current_gen_routes = std::vector<Route>(K);
            new_gen_routes= std::vector<Route>(K);
            probabilities = std::vector<float>(K, 1/K);

            parallel_for(0, K, [=](unsigned int i){current_gen_routes[i] = (*nation).random_route();}, num_workers);

            // Assign to the best route just the first random route, it will be updated very soon.
            best_route.distance = current_gen_routes[0].compute_total_distance();
            best_route.ordering = current_gen_routes[0].ordering;

            picker = std::discrete_distribution<>(probabilities.begin(), probabilities.end());

            std::cout<< "Initial best distance: " << best_route.distance << "\n";
        }


        Route pick_route(bool allow_self_pick = true){ 
            //@TODO implement prevention from a route picking itself as partner
            if (print_times && subtask_time_analysis){ 
                utimer pick_route("Partner pickup: ");

                return current_gen_routes[picker(gen)]; 
            } else { 
                return current_gen_routes[picker(gen)]; 
            }
           
        }

        // Here we should do the more of the parallelization part. 
        void evolve_with_threads(int worker_idx, float mutation_rate){ 
            int first_target = worker_idx * int(K/num_workers) + (worker_idx < K % num_workers ? worker_idx : K % num_workers );
            int last_target = first_target + int(K/num_workers) + (worker_idx < K % num_workers? 1 : 0);

            //first = worker_idx * int(K/num_workers) + (worker_idx if worker_idx < K % num_workers else K % num_workers )
            //last = first + int(K/num_workers) + ( 1 if worker_idx < K%num_workers else 0 )
            for (int j = first_target; j < last_target; j++) { 
                Route mate_route = pick_route();
                int crossing_point = rand()% N;// rand()%2 ==1 ? rand()% N : -1;

                // int route_idx = worker_idx*K/num_workers + j;
                Route current_route = current_gen_routes[j];
                //std::cout << "Rotta letta" << "\n";

                Route new_route;
                
                // Randomly choose which PMX ordering to choose
                if(rand() % 2 == 1){
                    new_route = current_route.PMX(mate_route, crossing_point);
                } else { 
                    new_route = mate_route.PMX(current_route, crossing_point);
                }
                //std::cout << "PMX fatto" << "\n";

                // Randomly check if the mutation should happen or not. 
                if ( mutation_rate != 0 && (rand() % int(1/mutation_rate)) == 1){ 
                    if(print_times && subtask_time_analysis){
                        utimer mutationTime("Mutation time:");
                        new_route.mutate();
                    } else { 
                        new_route.mutate();
                    }
                }

                //std::cout << "MuTASSAO" << std::size(new_gen_routes) << "\n";
                float distance = new_route.compute_total_distance();

                new_gen_routes[j] = new_route;
                total_inverse_distance += 1/distance;
            }
            
        }

        void rank_all(){
            if(print_logs){ 
                std::cout << "Epoch." << 0 << ")  Best Distance:";
                std::cout << best_route.compute_total_distance() << '\n';
            }
            
            for (int j = 0; j < K; j++){  
                probabilities[j] = (1/probabilities[j])/total_inverse_distance;
                current_gen_routes[j].ordering = new_gen_routes[j].ordering;
                current_gen_routes[j].distance = new_gen_routes[j].distance;

                if (current_gen_routes[j].distance < best_route.distance){ 
                    best_route.distance = current_gen_routes[j].distance;
                    best_route.ordering = current_gen_routes[j].ordering;
                }
            }

            picker = std::discrete_distribution<>(probabilities.begin(), probabilities.end());
            total_inverse_distance = 0;
        }

        void evolve(int epochs, float mutation_rate){ 
            merge_em e(epochs, num_workers, [=](){rank_all();});
            std::vector<std::unique_ptr<ff_node>> workers_pool;

            for(int i = 0; i < num_workers; i++){
                workers_pool.push_back(make_unique<worker_SCMF>(i, mutation_rate, [=](int i, float mutation_rate){evolve_with_threads(i, mutation_rate);}));
            }

            ff_Farm<> farm(std::move(workers_pool));
            
            farm.add_emitter(e); //Same as a barrier with a callback.
            
            farm.remove_collector();
            farm.wrap_around();
            farm.run_and_wait_end();

            rank_all();
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
        sstm << "LOGS/FF/E"<< epochs << "G"<< genes_num << "W" << num_workers << "_" << ct << ".txt";
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
    std::cout << "Best solution length: "<<(*genome).best_route.distance << "\n";
}

void time_tests(){ 
    if (print_times ){ 
        utimer evolution("Single evolution loop (including ranking):");
        // TODO evakyate the single tasks here.
    }
}