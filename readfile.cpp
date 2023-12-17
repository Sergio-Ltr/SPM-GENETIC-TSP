#include <iostream>
#include <fstream>
#include <string>

#include <bits/stdc++.h>

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

    int distance(City dest) { 
        return pow(pow(x - dest.x, 2) + pow(y - dest.y, 2), 0.5);
    }
};

class Route { 
  public: 
    int N;
    std::vector<City> *cities;
    std::vector<int> ordering;
    float distance = 0;

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

    float totalDistance(){ //@TODO save values in a symmetric matrix??? 
        distance = 0;
        for(int i=0; i < N; i++){
            float delta_d = cities->at(ordering[i]).distance(cities->at(ordering[(i+1)%N]));
            distance = distance + delta_d;

            std::cout << (cities->at(ordering[i])).id ;
            std::cout << " --> ";
            std::cout << (cities->at(ordering[(i+1)%N])).id ;
            std::cout << " == " <<  delta_d << '\n';
        }
        return distance;
    }

    Route PMX(Route mating_route, int crossing_point){
        std::vector<int> child_ordering;
        child_ordering.reserve(N);

        std::vector<int>::iterator current_begin = ordering.begin();
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

    void printCitiesId () { 
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
        // std::cout << N << '\n';

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

    City getCity(int id){ 
        return cities[id];
    }

    void printCitiesId () { 
        for(unsigned int i = 0; i < N; i++){
            std::cout << cities[i].id << '\n';
        }
    }

    Route randomRoute() { 
        return Route(N, &cities);
    }
};


int main (int argc, char* argv[]){
    // This code allows to read the .tsp files skipping the fixed-lenght headers (7 lines).
    // The remaining lines just corrresponds to the coordinates of the cities in the nation.
    std::string filename = argc > 1 ? argv[1] : "wi29.tsp";
    std::cout << filename << '\n';
    
    Nation nation = Nation(filename);    

    //nation.printCitiesId();
    //std::cout << nation.N << '\n';
    //std::cout << nation.getCity(1).distance(nation.getCity(9))  << '\n';
    std::cout << nation.randomRoute().totalDistance() << '\n';
    //std::cout << (c->at(0)).x << '\n';
    Route route1 =  nation.randomRoute();
    Route route2 =  nation.randomRoute();
    

    route1.printCitiesId();
    std::cout << "++++++++++++++++++++++++++++++++++++++++" << '\n';
    route2.printCitiesId();
    std::cout << "========================================" << '\n';
    Route child = (route2.PMX(route1, 8));
    child.printCitiesId();
    //std::cout << child.ordering[0] <<'\n';
    return 0;      
}