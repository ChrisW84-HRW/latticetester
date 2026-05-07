#include <cstdlib>
#include <iostream>

int main() {
    std::cout << "##########################\nTesting Coordinates.h \n########################## \n \n";
    int result = system("./testCoordinates");    
    // std::cout << "##########################\nTesting Util.h \n########################## \n \n";
    // result = system("./testUtil");    
    std::cout << "##########################\nTesting BasisConstruction.h \n########################## \n \n";
    result = system("./testBasisConstruction");
    std::cout << "##########################\nTesting Rank1Lattice.h \n########################## \n \n";
    result = system("./testRank1Lattice");        
    std::cout << "##########################\nTesting IntLattice.h \n########################## \n \n";
    result = system("./testIntLattice");    
    std::cout << "##########################\nTesting ReducerStatic.h \n########################## \n \n";
    result = system("./testReducerStatic");    
    std::cout << "##########################\nTesting ReducerBB.h \n########################## \n \n";
    result = system("./testReducerBB");    
    std::cout << "##########################\nTesting FigureOfMerit.h and FigureOfMeritDual.h \n########################## \n \n";
    result = system("./testFigureOfMerit");
    return result;
}
