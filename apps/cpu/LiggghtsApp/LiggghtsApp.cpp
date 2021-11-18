#include <iostream>
#include <string>
#include <memory>

#include "VirtualFluids.h"

//#include "lammps.h"
//#include "input.h"
//#include "atom.h"
//#include "modify.h"
//#include "fix_lb_coupling_onetoone.h"

#include "LiggghtsCouplingCoProcessor.h"

#include "LiggghtsCouplingWrapper.h"


using namespace std;


int main(int argc, char *argv[])
{
    SPtr<Communicator> comm = MPICommunicator::getInstance();
    LiggghtsCouplingWrapper wrapper(argv, (MPI_Comm)comm->getNativeCommunicator());

	//LAMMPS_NS::LAMMPS *lmp;
 //   // custom argument vector for LAMMPS library
 //   const char *lmpargv[] {"liblammps", "-log", "none"};
 //   int lmpargc = sizeof(lmpargv)/sizeof(const char *);

 //   // explicitly initialize MPI
 //   MPI_Init(&argc, &argv);

 //   // create LAMMPS instance
 //   lmp = new LAMMPS_NS::LAMMPS(lmpargc, (char **)lmpargv, MPI_COMM_WORLD);
 //   lmp->input->file("in.lbdem");
 //   //lmp->input->one("run 1");
 //   
 //   //# Try extracting a global value
 //   //    print("")
 //   //    print("Attempting to get the number of atoms in simulation")
 //   //    numAtoms = lmp.extract_global("natoms", 0)
 //   //    print("natoms =", numAtoms)

 //   //    # Try extracting atom's positions
 //   //    print("")
 //   //    print("Attempting to get the atom's positions")
 //   //    pos = lmp.extract_atom("x",3)
 //   //    for k in range(0,numAtoms):
 //   //        print("Pos[%i] = [%f, %f, %f]" % (k, pos[k][0], pos[k][1], pos[k][2]))

 //   LAMMPS_NS::FixLbCouplingOnetoone 
 //       *couplingFix 
 //       = dynamic_cast<LAMMPS_NS::FixLbCouplingOnetoone*>
 //       (lmp->modify->find_fix_style("couple/lb/onetoone",0));

 //   cout << "test1\n";
 //   
 //   //double **t_liggghts = couplingFix->get_torque_ptr();
 //   cout << "test2\n";

 //   lmp->input->one("run 9 upto");

 //   for (int step = 0; step < 10; step++)
 //   {
 //       

 //       int numAtoms = lmp->atom->natoms;

 //       //double** pos = (double**)lmp->atom->extract("x");
 //       double** pos = lmp->atom->x;
 //       
 //       //double* forceX = lmp->atom->fx;

 //       for (int i = 0; i < numAtoms; i++)
 //       {
 //           double **f_liggghts = couplingFix->get_force_ptr();
 //           double** force = lmp->atom->f;
 //           cout << "Pos[" << i << "] = [" << pos[i][0] << ", " << pos[i][1] << ", " << pos[i][2] << "]\n";
 //           cout << "Force1[" << i << "] = [" << f_liggghts[i][0] << ", " << f_liggghts[i][1] << ", " << f_liggghts[i][2] << "]\n";
 //           f_liggghts[i][0] += 0;
 //           f_liggghts[i][1] += 0;
 //           f_liggghts[i][2] += 500;
 //           cout << "Force2[" << i << "] = [" << force[i][0] << ", " << force[i][1] << ", " << force[i][2] << "]\n";
 //       }

 //       couplingFix->comm_force_torque();

 //       lmp->input->one("run 10000");
 //      
 //   }

 //   // delete LAMMPS instance
 //   delete lmp;

 //   // stop MPI environment
 //   MPI_Finalize();
    return 0;
}
