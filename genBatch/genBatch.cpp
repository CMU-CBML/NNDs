#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib> // For std::atoi

void generate_job_file(int job_number) {
    std::ofstream file;
    std::string filename = "Single" + std::to_string(job_number) + ".job";
    file.open(filename);

    file << "#!/bin/sh\n\n";
    file << "#SBATCH -p RM-shared\n";
    file << "#SBATCH --ntasks-per-node=16\n";
    file << "#SBATCH --job-name=S" << job_number << "\n";
    file << "#SBATCH --output=S" << job_number << ".out\n";
    file << "#SBATCH --time=48:00:00\n";
    // file << "#SBATCH --mail-type=ALL\n";
    // file << "#SBATCH --mail-user=your_email@domain.com\n\n";
    file << "mkdir ../io2D_single" << job_number << "/\n";
    file << "mkdir ../io2D_single" << job_number << "/outputs/\n";
    file << "cd ../2D_solver_coarse_smallGC_04182024/\n";
    file << "/ocean/projects/eng170006p/ussqww/petsc/arch-linux-c-opt/bin/mpiexec -np 16 ./2DNG 1 350000 ../io2D_single" << job_number << "/\n\n";
    file << "#-- exit\n#\n";

    file.close();
    std::cout << "Generated file: " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " [end number]" << std::endl;
        return 1;
    }

    int end = std::atoi(argv[1]);  // Convert command-line input to integer

    for (int i = 0; i <= end; i++) {
        generate_job_file(i);
    }

    return 0;
}