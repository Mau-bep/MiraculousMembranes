#include <iostream>
#include <omp.h>

int main() {
    std::cout<<"Welcome to the parallel version with "<<omp_get_num_threads()<< " threads\n\n";
    const size_t N = 1000000000000;
    int sum = 0;

    #pragma omp parallel for reduction(+:sum)
    for (size_t i = 0; i < N; ++i) {
        sum += i;
    }

    std::cout << "Sum: " << sum << std::endl;

    return 0;
}



