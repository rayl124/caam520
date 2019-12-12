#include <iostream>
#include <mpi.h>

class MPISession
{
public:
  MPISession(int *argc, char ***argv)
  {
    MPI_Init(&argc, &argv);
  }

  ~MPISession()
  {
    MPI_Finalize();
  }
};

int main(int argc, char **argv)
{
  MPISession session(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  std::cout << "hello, world from rank " << rank << "/" << size << std::endl;

  return 0;
}
