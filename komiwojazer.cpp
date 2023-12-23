#include <iostream>
#include <stdlib.h>
#include <utility>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>
#include <math.h>

std::vector<std::pair<int, int>> ChooseInput(int argc, char* argv[]);
std::vector<std::pair<int, int>> GenerateRandomGraph(int size);
std::vector<std::pair<int, int>> GetGraphFromFile(std::string file);
std::vector<std::vector<double>> Matrix(std::vector<std::pair<int, int>> coords, int size);
std::pair<double, std::vector<int>> greedyAlgorithmTSP(int size, std::vector<std::vector<double>> matrix);
std::pair<double, std::vector<int>> AntsAlgorithm(int ants, int iterations, float alpha, float beta, float p, float Q, float c, std::vector<std::vector<double>> &Matrix);
std::pair<std::vector<float>, std::vector<int>> DataForAnts();
void OutputToTerminal(double allDistance, std::vector<int> path, int size);
void SaveOutputGreedy(double allDistance, std::vector<int> path);
void SaveOutputAnts(double allDistance, std::vector<int> path, std::pair<std::vector<float>, std::vector<int>> dataAnts);
void probabilityValue(int current, float alpha, float beta, std::vector<std::vector<double>> &Matrix, std::vector<long double> &probabilityAnts, std::vector<std::vector<float>> &trailIntensity, std::vector<int> &allowAnt);
int next(int size, std::vector<std::vector<double>> &Matrix, std::vector<long double> &probability, std::vector<int> &allowAnt);
int ChooseAnAlgorithm(int choice);
double GetDistance(std::vector<std::pair<int, int>> coords, int i, int j);
long double randomLongDouble();

int main(int argc, char *argv[]){
    std::srand(time(nullptr));
    std::vector<std::pair<int, int>> coords;
    std::vector<std::vector<double>> matrix;
    std::pair<double, std::vector<int>> output;
    std::pair<std::vector<float>, std::vector<int>> dataAnts;
    std::vector<float> FloatDataForAnts;
    std::vector<int> IntDataForAnts;
    std::vector<int> path;
    double allDistance;
    int choice;
    
    
    coords = ChooseInput(argc, argv);
    int size = coords.size();
    matrix = Matrix(coords, size);

    choice = ChooseAnAlgorithm(choice);

    if(choice){
        dataAnts = DataForAnts();
        FloatDataForAnts = dataAnts.first;
        IntDataForAnts = dataAnts.second;
        output = AntsAlgorithm(IntDataForAnts.front(), IntDataForAnts.back(), FloatDataForAnts.front(), FloatDataForAnts[1], FloatDataForAnts[2], FloatDataForAnts[3], FloatDataForAnts.back(), matrix);
        SaveOutputAnts(output.first, output.second, dataAnts);
    } else {
        output = greedyAlgorithmTSP(size, matrix);
        SaveOutputGreedy(output.first, output.second);
    }
    

    OutputToTerminal(output.first, output.second, size);
    return 0;
}



int ChooseAnAlgorithm(int choice){
    std::cout<<"Choose a method:"<<std::endl;
    std::cout<<"(0) Greedy algorithm"<<std::endl;
    std::cout<<"(1) Ant-colony algorithm"<<std::endl;
    std::cin>>choice;
    return choice;
}

std::pair<std::vector<float>, std::vector<int>> DataForAnts() {
    std::vector<float> FloatData;
    std::vector<int> IntData;
    float alpha, beta, p, Q, c;
    int ants, iterations, choice;

    std::cout<<"(0) Get default data for ants"<<std::endl;
    std::cout<<"(1) Enter your data for ants"<<std::endl;
    std::cin>>choice;

    if((int)choice){
        std::cout<<"Enter Alpha: "<<std::endl;
        std::cin>>alpha;
        FloatData.push_back(alpha);
        std::cout<<"Enter Beta:"<<std::endl;
        std::cin>>beta;
        FloatData.push_back(beta);
        std::cout<<"Enter p: "<<std::endl;
        std::cin>>p;
        FloatData.push_back(p);
        std::cout<<"Enter Q: "<<std::endl;
        std::cin>>Q;
        FloatData.push_back(Q);
        std::cout<<"Enter c: "<<std::endl;
        std::cin>>c;
        FloatData.push_back(c);

        std::cout<<"Enter number of ants: "<<std::endl;
        std::cin>>ants;
        IntData.push_back(ants);
        std::cout<<"Enter number of iterations: "<<std::endl;
        std::cin>>iterations;
        IntData.push_back(iterations);
    } else {
        FloatData.push_back(1);
        FloatData.push_back(5);
        FloatData.push_back(0.5);
        FloatData.push_back(100);
        FloatData.push_back(1);

        IntData.push_back(100);
        IntData.push_back(100);
    }
    return std::pair<std::vector<float>, std::vector<int>> (FloatData, IntData);
}

std::vector<std::pair<int, int>> ChooseInput(int argc, char* argv[]){
    std::string file;
    std::vector<std::pair<int, int>> coords;
    int size, choice;

    if (argc==2){
        file = argv[1];
        coords = GetGraphFromFile(file);
    } else if (argc == 1) {
        std::cout<<"Input graph size: "<<std::endl;
        std::cin>>size;
        coords = GenerateRandomGraph(size);
    }

    return coords;
}

std::vector<std::pair<int, int>> GenerateRandomGraph(int size) {
    std::ofstream outputFile("random.txt");
    std::vector<std::pair<int, int>> coords;
    std::pair<int, int> newPair;
    int coord1, coord2;


    for (int i = 0; i < size; i++) { // setting coordinates for points
        coord1 = rand() % 10000;
        coord2 = rand() % 10000;
        newPair.first = coord1;
        newPair.second = coord2;
        std::cout << newPair.first << " ";
        std::cout << newPair.second << "\n";
        coords.push_back(newPair);
    }
    if (outputFile.is_open()) {
        outputFile << size << std::endl;
        for(int i = 0; i<size; i++){
            outputFile << i+1 << " " << coords[i].first <<" "<< coords[i].second << std::endl;
        }
    }
    return coords;
}
std::vector<std::pair<int, int>> GetGraphFromFile(std::string file) {
    std::ifstream inputFile(file);
    std::vector<std::pair<int, int>> coords;
    std::string coord1, coord2;
    int size, nr; 
    if (inputFile.is_open()) {
        inputFile >> size;
        while(inputFile >> nr >> coord1 >> coord2){
            coords.push_back(std::pair<int,int>(std::stoi(coord1), std::stoi(coord2)));
        }
        return coords;
    }
    return coords; 

}
double GetDistance(std::vector<std::pair<int, int>> coords, int i, int j){
    return sqrt((coords[j].first - coords[i].first)*(coords[j].first - coords[i].first) + (coords[j].second - coords[i].second)*(coords[j].second - coords[i].second));
}
std::vector<std::vector<double>> Matrix(std::vector<std::pair<int, int>> coords, int size){
    std::vector<std::vector<double>> matrix;
    for ( int i = 0; i < size; i++) {
        matrix.push_back(std::vector<double>(size, -1));
        for (int j = 0; j < size; j++){
            if (i != j){
                matrix[i][j] = GetDistance(coords, i, j);
            }
        }
    }
    return matrix;
}
std::pair<double, std::vector<int>> greedyAlgorithmTSP(int size, std::vector<std::vector<double>> matrix){
    std::vector<int> path;
    std::vector<int> visited(size, 0);
    double allDistance = 0.0;
    int currentVertex, destination;
    double currentValue;

    visited[0] = 1;
    path.push_back(0);
    
    while(path.size()<size){
        currentValue = 100000000000000000;
        currentVertex = path.back();
        for (int i = 0; i < size; i++){
            if(currentVertex != i && !visited[i] && matrix[currentVertex][i] < currentValue) {
                currentValue = matrix[currentVertex][i];
                destination = i;
            }
        }
        path.push_back(destination);
        visited[destination] = 1;
        allDistance += currentValue;
    }
    path.push_back(0);
    allDistance += matrix[0][destination];

    return std::pair<double, std::vector<int>>(allDistance, path);
}
void SaveOutputAnts(double allDistance, std::vector<int> path, std::pair<std::vector<float>, std::vector<int>> dataAnts) {
    std::ofstream output("output.txt");
    std::vector<float> FloatDataForAnts;
    std::vector<int> IntDataForAnts;
    FloatDataForAnts = dataAnts.first;
    IntDataForAnts = dataAnts.second;
    if (output.is_open()) {
        output << "------"<<std::endl;
        output << "ants: "<< IntDataForAnts.front() << ", iterations: "<< IntDataForAnts.back() <<", Alpha: "<< FloatDataForAnts.front() << ", Beta: "<<FloatDataForAnts[1]<<", p: "<< FloatDataForAnts[2]<< ", Q: " << FloatDataForAnts[3] <<", c: "<<FloatDataForAnts.back() ;
        output << std::endl;
        output << "Distance: "<< allDistance <<std::endl;
        for (int i = 0; i < path.size(); i++){
            output << path[i] + 1<< ", ";
        }
        output << std::endl;
    } else {
        std::cout<<"Could not open result file"<<std::endl;
    }
}
void SaveOutputGreedy(double allDistance, std::vector<int> path) {
    std::ofstream output("output.txt");
    if (output.is_open()) {
        output << "------"<<std::endl;
        output << "Distance: "<< allDistance <<std::endl;
        for (int i = 0; i < path.size(); i++){
            output << path[i] + 1<< ", ";
        }
        output << std::endl;
    } else {
        std::cout<<"Could not open result file"<<std::endl;
    }
}
void OutputToTerminal(double allDistance, std::vector<int> path, int size) {
    std::cout << "Length of path: " << allDistance << std::endl << "Path: ";
    for ( int i = 0; i <= size; i++ ){
        std::cout << path[i] + 1 << " ";
    }
}

//Ants
long double randomLongDouble() {
    return (long double)(rand()) / (long double) (RAND_MAX);
}

int next(int size, std::vector<std::vector<double>> &Matrix, std::vector<long double> &probabilityAnts, std::vector<int> &allowAnt){
    long double random = randomLongDouble();
    for (int r = 0; r < size; r++) {
        if (random < probabilityAnts[r]) {
            return r;
        } else if ( random >= probabilityAnts[r]) {
            random -= probabilityAnts[r];
        }
    }
    for (int j = 0; j < size; j++) {
        if(allowAnt[j] == 1) {
            return j;
        }
    }
    return -1;
}
void probabilityValue(int size, int current, float alpha, float beta, std::vector<std::vector<double>> &Matrix, std::vector<long double> &probabilityAnts, std::vector<std::vector<float>> &trailIntensity, std::vector<int> &allowAnt){
    long double valueEnd, nValue;
    long double denominator = 0.0;

    for (int i = 0; i < size; i++){
        if(allowAnt[i] == 0 || i == current) {
            if(i == current){
                allowAnt[i] = 0;
            }
        }
        else {
            nValue = 1.0/Matrix[current][i];
            if(pow(trailIntensity[current][i], alpha) * pow(nValue, beta)>=0){
                denominator+= pow(trailIntensity[current][i], alpha) * pow(nValue, beta);
            }
        }
    }
    for (int j = 0; j < size; j++){
        if(allowAnt[j] == 0){
            probabilityAnts[j] = 0;
        } else if (denominator <= 0) {
            probabilityAnts[j] = 0;
        } else {
            probabilityAnts[j] = pow(trailIntensity[current][j], alpha)/pow(Matrix[current][j], beta)/denominator;
        }
    }
}
std::pair<double, std::vector<int>> AntsAlgorithm(int ants, int iterate, float alpha, float beta, float p, float Q, float c, std::vector<std::vector<double>> &Matrix){
    int size = Matrix.size();
    int point;
    float divide;
    float dist = 0.0;
    double allDist = 100000000000000000;
    std::vector<std::vector<float>> trailIntensity (size, std::vector<float>(size, c));
    std::vector<std::vector<float>> trailIntensity2 (size, std::vector<float>(size, 0));

    std::vector<std::vector<int>> allowAnt (ants, std::vector<int>(size, 1));
    std::vector<std::vector<int>> antsTrail (ants, std::vector<int>());

    std::vector<std::vector<long double>> probability (ants, std::vector<long double> (size, 0));
    std::vector<int> path;

    for (int i = 0; i < iterate; i++){

        for (int x = 0; x < ants; x++){
            antsTrail[x].push_back(rand() % size);
            allowAnt[x][antsTrail[x][0]] = 0;
        }


        for(int k = 0; k < size -1; k++) {
            for (int x = 0; x < ants; x++){
                probabilityValue(size, antsTrail[x].back(), alpha, beta, Matrix, probability[x], trailIntensity, allowAnt[x]);
                antsTrail[x].push_back(point = next(size, Matrix, probability[x], allowAnt[x]));
                allowAnt[x][point] = 0;
            }
        }



        for (int y = 0; y < ants; y++) {
            dist = 0.0;
            antsTrail[y].push_back(antsTrail[y][0]);
            for (int x = 0; x < size; x++) {
                dist += Matrix[antsTrail[y][x]][antsTrail[y][x + 1]];
            }
            for(int x = 0; x < size - 1; x++) {
                divide = Q/dist;
                trailIntensity2[antsTrail[y][x]][antsTrail[y][x+1]] += divide;
                trailIntensity2[antsTrail[y][x+1]][antsTrail[y][x]] += divide;
            }
            if (dist < allDist){
                allDist = dist;
                path = antsTrail[y];
            }
        }

        for (int i = 0; i < ants; i++) {
            antsTrail[i].clear(); 
            allowAnt[i] = std::vector<int> (size, 1);    
        }

        for (int x = 0; x < size; x++) {
            for (int y = x+1; y < size; y++) {
                trailIntensity2[y][x] = 0.0;

                trailIntensity[x][y] = p*trailIntensity[x][y] + trailIntensity2[x][y];
                trailIntensity2[x][y] = 0.0;
                trailIntensity[y][x] = trailIntensity[x][y];
            }
        }
    }
    return std::pair<double, std::vector<int>>(allDist, path);
}

