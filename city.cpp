class City {        // The class
  public:          // Access specifier
    int id;    // Index in the .tsp file
    float x;  // Latitude
    float y;  // Longitude
    
    City(int index, float longitude, float latitude) { // Constructor with parameters
      id = index;
      x = latitude;
      y = longitude;

    }
};
