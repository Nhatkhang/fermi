/* frees all the memory */

#include "fermi.h"

int ferfini(void)
{

    // eliminamos las estructuras de comunicadores
    // e intercomunicadores
//    commdom_delete(); 

    SlepcFinalize();
    
    return 0;
}
